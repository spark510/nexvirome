#!/usr/bin/env python3
"""
Profile the virome classification pipeline to identify bottlenecks.

This script runs the full pipeline with timing and memory profiling
to identify performance bottlenecks.
"""
from __future__ import annotations

import time
import sys
import argparse
from pathlib import Path
from typing import Optional, Callable, Any
import pandas as pd

from ..core import FilterCriteria, log_info, set_verbose
from ..taxonomy import TaxonomyDB
from ..alignment import AlignmentParser, MaskingFilter


class Timer:
    """Context manager for timing code blocks."""

    def __init__(self, name: str, verbose: bool = True):
        self.name = name
        self.verbose = verbose
        self.start_time = None
        self.elapsed = None

    def __enter__(self):
        if self.verbose:
            log_info(f"\n{'='*70}")
            log_info(f"⏱️  Starting: {self.name}")
            log_info(f"{'='*70}")
        self.start_time = time.time()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.elapsed = time.time() - self.start_time
        if self.verbose:
            log_info(f"✅ Completed: {self.name}")
            log_info(f"⏱️  Time: {self.elapsed:.2f} seconds ({self.elapsed/60:.2f} minutes)")
            log_info(f"{'='*70}")


def get_memory_usage() -> float:
    """Get current memory usage in MB."""
    try:
        import psutil
        process = psutil.Process()
        return process.memory_info().rss / (1024 ** 2)
    except ImportError:
        return 0.0


def profile_step(name: str, func: Callable, *args, **kwargs) -> tuple[Any, float, float]:
    """
    Profile a single pipeline step.

    Returns:
        (result, elapsed_time, memory_used)
    """
    mem_before = get_memory_usage()
    start_time = time.time()

    log_info(f"\n{'='*70}")
    log_info(f"⏱️  {name}")
    log_info(f"{'='*70}")

    result = func(*args, **kwargs)

    elapsed = time.time() - start_time
    mem_after = get_memory_usage()
    mem_used = mem_after - mem_before

    log_info(f"✅ Completed: {name}")
    log_info(f"   Time: {elapsed:.2f}s ({elapsed/60:.2f}m)")
    log_info(f"   Memory: {mem_after:.1f} MB (Δ {mem_used:+.1f} MB)")

    return result, elapsed, mem_used


def profile_pipeline(
    r1_file: str,
    r2_file: str,
    taxonomy_db: str,
    mask_file: str,
    min_identity: float = 0.8,
    min_length: int = 50,
    min_coverage: float = 0.1,
    verbose: bool = True,
):
    """
    Run full pipeline with profiling.
    """
    if verbose:
        set_verbose(True)

    timings = {}
    memory_usage = {}

    print("\n" + "=" * 70)
    print("PIPELINE PROFILING")
    print("=" * 70)
    print(f"R1: {r1_file}")
    print(f"R2: {r2_file}")
    print(f"Taxonomy: {taxonomy_db}")
    print(f"Mask: {mask_file}")
    print("=" * 70)

    # Step 1: Load Taxonomy
    def load_tax():
        return TaxonomyDB.from_sqlite(taxonomy_db, root_taxid=10239)

    tax, t, m = profile_step("1. Load Taxonomy Database", load_tax)
    timings["1_load_taxonomy"] = t
    memory_usage["1_load_taxonomy"] = m

    # Step 2: Parse R1
    def parse_r1():
        parser = AlignmentParser(normalize_headers=True)
        df = parser.parse(Path(r1_file))
        df['strand'] = '+'
        return df

    r1_hits, t, m = profile_step("2. Parse R1 Alignments", parse_r1)
    timings["2_parse_r1"] = t
    memory_usage["2_parse_r1"] = m
    log_info(f"   R1 hits: {len(r1_hits):,}")

    # Step 3: Parse R2
    def parse_r2():
        parser = AlignmentParser(normalize_headers=True)
        df = parser.parse(Path(r2_file))
        df['strand'] = '-'
        return df

    r2_hits, t, m = profile_step("3. Parse R2 Alignments", parse_r2)
    timings["3_parse_r2"] = t
    memory_usage["3_parse_r2"] = m
    log_info(f"   R2 hits: {len(r2_hits):,}")

    # Step 4: Combine
    def combine():
        return pd.concat([r1_hits, r2_hits], ignore_index=True)

    combined_hits, t, m = profile_step("4. Combine R1 + R2", combine)
    timings["4_combine"] = t
    memory_usage["4_combine"] = m
    log_info(f"   Combined hits: {len(combined_hits):,}")

    # Step 5: Quality Filter
    def quality_filter():
        parser = AlignmentParser()
        criteria = FilterCriteria(
            min_identity=min_identity,
            min_alignment_length=min_length,
            max_evalue=1e-3,
            min_query_coverage=0.5,
        )
        return parser.filter(combined_hits, criteria)

    filtered_hits, t, m = profile_step("5. Quality Filtering", quality_filter)
    timings["5_quality_filter"] = t
    memory_usage["5_quality_filter"] = m
    log_info(f"   Filtered hits: {len(filtered_hits):,} ({len(filtered_hits)/len(combined_hits)*100:.1f}%)")

    # Step 6: Load Mask
    def load_mask():
        return MaskingFilter.from_bed_file(mask_file)

    masking_filter, t, m = profile_step("6. Load Masking Filter", load_mask)
    timings["6_load_mask"] = t
    memory_usage["6_load_mask"] = m

    # Step 7: Masking Filter
    def masking():
        return masking_filter.filter_by_unmasked_coverage(filtered_hits, min_coverage=min_coverage)

    result, t, m = profile_step("7. Apply Masking Filter", masking)
    timings["7_masking_filter"] = t
    memory_usage["7_masking_filter"] = m
    log_info(f"   Passed: {len(result.passed):,} hits ({result.n_passed_targets} targets)")

    # Step 8: LCA Classification
    def lca_classify():
        query_groups = result.passed.groupby('query')
        lca_results = []

        for query, group in query_groups:
            taxids = group['taxid'].unique().tolist()
            lca_taxid = tax.compute_lca(taxids)

            if lca_taxid and lca_taxid > 0:
                qlen = int(group.iloc[0]['qlen']) if 'qlen' in group.columns else 100
                lca_results.append({
                    'query': query,
                    'lca_taxid': lca_taxid,
                    'qlen': qlen,
                    'read_count': 1,
                })

        return pd.DataFrame(lca_results)

    lca_df, t, m = profile_step("8. LCA Classification", lca_classify)
    timings["8_lca_classification"] = t
    memory_usage["8_lca_classification"] = m
    log_info(f"   Classified: {len(lca_df):,} queries")

    # Get LCA cache statistics
    cache_stats = tax.get_lca_cache_stats()
    log_info(f"   LCA Cache: {cache_stats['hits']:,} hits, {cache_stats['misses']:,} misses")
    log_info(f"   LCA Cache Hit Rate: {cache_stats['hit_rate_percent']:.1f}%")
    log_info(f"   LCA Cache Size: {cache_stats['cache_size']:,} unique taxid combinations")

    # Print Summary
    print("\n" + "=" * 70)
    print("PROFILING SUMMARY")
    print("=" * 70)

    total_time = sum(timings.values())
    total_memory = sum(memory_usage.values())

    # Sort by time
    sorted_timings = sorted(timings.items(), key=lambda x: x[1], reverse=True)

    print("\n⏱️  TIME BREAKDOWN:")
    print(f"{'Step':<30} {'Time (s)':<12} {'%':<8} {'Time (m)':<10}")
    print("-" * 70)
    for step, elapsed in sorted_timings:
        percentage = (elapsed / total_time * 100) if total_time > 0 else 0
        minutes = elapsed / 60
        print(f"{step:<30} {elapsed:>10.2f}s {percentage:>6.1f}% {minutes:>8.2f}m")
    print("-" * 70)
    print(f"{'TOTAL':<30} {total_time:>10.2f}s {'100.0%':>7} {total_time/60:>8.2f}m")

    print("\n💾 MEMORY USAGE:")
    print(f"{'Step':<30} {'Memory (MB)':<15}")
    print("-" * 70)
    for step, mem in sorted(memory_usage.items(), key=lambda x: x[1], reverse=True):
        if mem > 0:
            print(f"{step:<30} {mem:>13.1f} MB")
    print("-" * 70)
    print(f"{'TOTAL DELTA':<30} {total_memory:>13.1f} MB")
    print(f"{'CURRENT':<30} {get_memory_usage():>13.1f} MB")

    # Identify bottlenecks
    print("\n🔍 BOTTLENECK ANALYSIS:")
    bottlenecks = [step for step, elapsed in sorted_timings if elapsed / total_time > 0.15]
    if bottlenecks:
        print("Steps taking >15% of total time:")
        for step in bottlenecks:
            percentage = (timings[step] / total_time * 100)
            print(f"  ⚠️  {step}: {timings[step]:.2f}s ({percentage:.1f}%)")
    else:
        print("  ✅ No major bottlenecks detected (all steps <15% of total time)")

    # Data size analysis
    print("\n📊 DATA SIZE ANALYSIS:")
    print(f"  Input R1: {len(r1_hits):,} hits")
    print(f"  Input R2: {len(r2_hits):,} hits")
    print(f"  Combined: {len(combined_hits):,} hits")
    print(f"  After quality filter: {len(filtered_hits):,} hits ({len(filtered_hits)/len(combined_hits)*100:.1f}%)")
    print(f"  After masking: {len(result.passed):,} hits ({len(result.passed)/len(filtered_hits)*100:.1f}%)")
    print(f"  LCA classified: {len(lca_df):,} queries")

    # Optimization suggestions
    print("\n💡 OPTIMIZATION SUGGESTIONS:")

    # Check parsing bottleneck
    parse_time = timings.get("2_parse_r1", 0) + timings.get("3_parse_r2", 0)
    if parse_time / total_time > 0.3:
        print("  • Parsing is slow (>30% of time):")
        print("    - Consider using pre-filtered alignments")
        print("    - Use faster storage (SSD)")
        print("    - Check file format (binary formats may be faster)")

    # Check filtering bottleneck
    filter_time = timings.get("5_quality_filter", 0)
    if filter_time / total_time > 0.2:
        print("  • Quality filtering is slow (>20% of time):")
        print("    - Consider stricter thresholds to reduce data size early")
        print("    - Profile pandas operations")

    # Check masking bottleneck
    mask_time = timings.get("7_masking_filter", 0)
    if mask_time / total_time > 0.2:
        print("  • Masking filter is slow (>20% of time):")
        print("    - Coverage calculation may be expensive")
        print("    - Consider pre-filtering targets")
        print("    - Check mask file size")

    # Check LCA bottleneck
    lca_time = timings.get("8_lca_classification", 0)
    if lca_time / total_time > 0.2:
        print("  • LCA classification is slow (>20% of time):")

        # Check cache effectiveness
        if cache_stats['hit_rate_percent'] < 50:
            print(f"    ⚠️  Cache hit rate is low ({cache_stats['hit_rate_percent']:.1f}%)")
            print("    - Most queries have unique taxid combinations")
            print("    - Consider grouping queries by taxid sets before LCA")
        elif cache_stats['hit_rate_percent'] > 80:
            print(f"    ✅ Cache hit rate is good ({cache_stats['hit_rate_percent']:.1f}%)")
            print("    - Many queries share common taxid combinations")
            print("    - Remaining time is likely in first-time calculations")
        else:
            print(f"    ℹ️  Cache hit rate is moderate ({cache_stats['hit_rate_percent']:.1f}%)")

        print(f"    - {cache_stats['cache_size']:,} unique taxid combinations cached")
        print("    - Profile taxonomy tree traversal for further optimization")
        print("    - Check if queries have many unique taxids per query")

    # Memory suggestions
    if total_memory > 8000:  # > 8 GB
        print("  • High memory usage (>8 GB):")
        print("    - Consider processing in batches")
        print("    - Use data type optimization (int32 instead of int64)")
        print("    - Free intermediate DataFrames explicitly")

    print("\n" + "=" * 70)
    print("PROFILING COMPLETE")
    print("=" * 70)

    return {
        "timings": timings,
        "memory": memory_usage,
        "data_sizes": {
            "r1_hits": len(r1_hits),
            "r2_hits": len(r2_hits),
            "combined": len(combined_hits),
            "filtered": len(filtered_hits),
            "passed": len(result.passed),
            "lca": len(lca_df),
        }
    }


def main():
    parser = argparse.ArgumentParser(description="Profile virome classification pipeline")

    parser.add_argument("--r1", required=True, help="R1 alignment file")
    parser.add_argument("--r2", required=True, help="R2 alignment file")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy database")
    parser.add_argument("--mask", required=True, help="Mask BED file")
    parser.add_argument("--min-identity", type=float, default=0.8)
    parser.add_argument("--min-length", type=int, default=50)
    parser.add_argument("--min-coverage", type=float, default=0.1)
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    try:
        profile_pipeline(
            r1_file=args.r1,
            r2_file=args.r2,
            taxonomy_db=args.taxonomy,
            mask_file=args.mask,
            min_identity=args.min_identity,
            min_length=args.min_length,
            min_coverage=args.min_coverage,
            verbose=args.verbose,
        )
        return 0
    except Exception as e:
        print(f"\n❌ Error: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
