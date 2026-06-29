#!/usr/bin/env python3
"""
OTU Table Merger - Merge multiple LCA results into OTU tables

This script merges LCA classification results from multiple samples
into unified OTU tables (sample × taxon matrices) at various taxonomic ranks.

Usage:
    python -m virome_classifier.cli.merge_otu \\
        --input-dir results/ \\
        --pattern "*_read_classification.csv" \\
        --taxonomy taxonomy.db \\
        --output otu_results/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict

from ..core import log_info, set_verbose
from ..taxonomy import TaxonomyDB
from ..reporting import create_otu_pipeline
from ..reporting.otu_table import build_abundance_table, export_otu_table


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Merge LCA results into OTU tables",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect LCA files in directory
  %(prog)s --input-dir results/ --taxonomy tax.db --output otu_tables/

  # Specify pattern for LCA files
  %(prog)s --input-dir results/ --pattern "*_lca.csv" \\
      --taxonomy tax.db --output otu_tables/

  # Generate tables at specific ranks
  %(prog)s --input-dir results/ --taxonomy tax.db \\
      --ranks genus species family --output otu_tables/

  # With filtering
  %(prog)s --input-dir results/ --taxonomy tax.db \\
      --min-count 10 --min-samples 2 --output otu_tables/
        """,
    )

    # Input options
    input_group = parser.add_argument_group("Input options")
    input_group.add_argument(
        "--input-dir",
        "-i",
        type=str,
        required=True,
        help="Directory containing LCA classification CSV files",
    )
    input_group.add_argument(
        "--pattern",
        "-p",
        type=str,
        default="*_read_classification.csv",
        help="Glob pattern for LCA files (default: *_read_classification.csv)",
    )
    input_group.add_argument(
        "--taxonomy",
        "-t",
        type=str,
        required=True,
        help="Taxonomy database (SQLite)",
    )
    input_group.add_argument(
        "--abundance-pattern",
        type=str,
        default="*_coverage_abundance.tsv",
        help="Glob pattern for per-sample coverage_abundance.tsv files. When "
             "found, ALSO merge their read_count and TPM into sample × taxon "
             "matrices at each rank (otu_table_<rank>_readcount.csv / _tpm.csv). "
             "Set to '' to skip abundance merging. (default: *_coverage_abundance.tsv)",
    )

    # Output options
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output directory for OTU tables",
    )
    output_group.add_argument(
        "--ranks",
        "-r",
        nargs="+",
        default=["genus", "species"],
        help="Taxonomic ranks for OTU tables (default: genus species)",
    )
    output_group.add_argument(
        "--phage-host",
        action="store_true",
        help="ALSO emit otu_table_phage_host.csv: phage with a known host are "
             "rolled up to their bacterial/archaeal host genus (phages_of_<Host>) "
             "via the same post-processor used by classify --phage-host-rollup; "
             "non-phage and host-unknown phage stay at species. This is an ADDED "
             "table — the per-rank tables above are unaffected.",
    )

    # Filtering options
    filter_group = parser.add_argument_group("Filtering options")
    filter_group.add_argument(
        "--min-count",
        type=int,
        default=10,
        help="Minimum total count across samples (default: 10)",
    )
    filter_group.add_argument(
        "--min-samples",
        type=int,
        default=1,
        help="Minimum number of samples where taxon must be present (default: 1)",
    )
    filter_group.add_argument(
        "--normalize",
        action="store_true",
        help="Normalize counts to relative abundances",
    )

    # Taxonomy options
    tax_group = parser.add_argument_group("Taxonomy options")
    tax_group.add_argument(
        "--root-taxid",
        type=int,
        default=10239,
        help="Root taxonomy ID (default: 10239 for Viruses)",
    )

    # Other options
    other_group = parser.add_argument_group("Other options")
    other_group.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )

    return parser.parse_args()


def find_lca_files(input_dir: str, pattern: str) -> Dict[str, Path]:
    """
    Find LCA classification files in directory.

    Args:
        input_dir: Input directory path
        pattern: Glob pattern for files

    Returns:
        Dictionary mapping sample names to file paths
    """
    input_path = Path(input_dir)

    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    # Find matching files
    lca_files = {}
    for file_path in input_path.glob(pattern):
        # Extract sample name from filename by stripping the known per-sample
        # output suffixes (longest first so '_read_classification' wins over
        # '_classification', and the abundance suffix is handled too).
        sample_name = file_path.stem
        for suffix in ["_read_classification", "_coverage_abundance",
                       "_coverage_summary", "_lca_classification",
                       "_classification", "_abundance", "_lca"]:
            if sample_name.endswith(suffix):
                sample_name = sample_name[: -len(suffix)]
                break

        lca_files[sample_name] = file_path

    return lca_files


def main():
    """Main entry point."""
    args = parse_args()

    # Set verbosity
    if args.verbose:
        set_verbose(True)

    # Print header
    print("=" * 70)
    print("OTU TABLE MERGER")
    print("=" * 70)

    try:
        # 1. Find LCA files
        log_info(f"\n🔍 Searching for LCA files in {args.input_dir}...")
        log_info(f"   Pattern: {args.pattern}")

        lca_files = find_lca_files(args.input_dir, args.pattern)

        if not lca_files:
            log_info(f"❌ No LCA files found matching pattern: {args.pattern}")
            return 1

        log_info(f"✅ Found {len(lca_files)} LCA files:")
        for sample_name, file_path in lca_files.items():
            log_info(f"   {sample_name}: {file_path.name}")

        # 2. Load taxonomy
        log_info(f"\n📚 Loading taxonomy from {args.taxonomy}...")
        tax = TaxonomyDB.from_sqlite(args.taxonomy, root_taxid=args.root_taxid)
        log_info(f"✅ Loaded {len(tax):,} taxa")

        # 3. Generate OTU tables
        log_info(f"\n📊 Generating OTU tables at ranks: {', '.join(args.ranks)}")
        log_info(f"   Min count: {args.min_count}")
        log_info(f"   Min samples: {args.min_samples}")
        log_info(f"   Normalize: {args.normalize}")

        otu_tables = create_otu_pipeline(
            lca_files={name: str(path) for name, path in lca_files.items()},
            tax=tax,
            output_dir=args.output,
            ranks=args.ranks,
            min_count=args.min_count,
            min_samples=args.min_samples,
            normalize=args.normalize,
        )

        # 3b. Optional ADDED table: phage rolled up to host genus.
        # Re-uses classify's post-processor (apply_phage_host_rollup) on the same
        # per-query CSVs, then runs them through the normal species-rank OTU
        # builder. The synthetic host-genus rows (negative pseudo-taxids labelled
        # "phages_of_<Host>") survive as their own leaves; everything else is
        # unchanged, so this is purely additive to the per-rank tables above.
        if getattr(args, "phage_host", False):
            import tempfile
            import pandas as pd
            from ..classification.phage_host_rollup import (
                build_phage_host_map, apply_phage_host_rollup)

            log_info("\n🦠 Building phage→host-genus rolled-up OTU table...")
            tax2host, phage_set = build_phage_host_map(args.taxonomy)

            with tempfile.TemporaryDirectory() as tmp:
                rolled_files = {}
                for name, path in lca_files.items():
                    df = pd.read_csv(path)
                    if "lca_taxid" in df.columns and len(df):
                        df = apply_phage_host_rollup(df, tax2host, phage_set)
                    rp = Path(tmp) / f"{name}_phagehost.csv"
                    df.to_csv(rp, index=False)
                    rolled_files[name] = str(rp)

                # Build the species-level table over the rolled-up CSVs INTO the
                # temp dir, so it can never clobber the real otu_table_species.csv
                # already written in step 3. The host-genus pseudo-leaves come
                # through as species-rank rows.
                ph_tables = create_otu_pipeline(
                    lca_files=rolled_files,
                    tax=tax,
                    output_dir=tmp,
                    ranks=["species"],
                    min_count=args.min_count,
                    min_samples=args.min_samples,
                    normalize=args.normalize,
                )
                tmp_sp = Path(tmp) / "otu_table_species.csv"
                ph_out = Path(args.output) / "otu_table_phage_host.csv"
                if tmp_sp.exists():
                    import shutil
                    shutil.copyfile(tmp_sp, ph_out)
            if "species" in ph_tables:
                otu_tables["phage_host"] = ph_tables["species"]
            log_info(f"✅ phage-host OTU table: {ph_out}")

        # 3c. Abundance matrices (read_count + TPM) from per-sample
        # coverage_abundance.tsv. These metrics are precomputed per taxon by the
        # coverage classifier, so we only assemble them across samples (taxon-rolled
        # to each rank). read assignment files (_read_classification.csv, .kraken)
        # are intentionally NOT merged — only taxon-level summaries are.
        if getattr(args, "abundance_pattern", ""):
            abund_files = find_lca_files(args.input_dir, args.abundance_pattern)
            if abund_files:
                log_info(f"\n📈 Merging abundance metrics (read_count, TPM) from "
                         f"{len(abund_files)} coverage_abundance.tsv files...")
                # taxid-level (raw) + each requested rank
                abund_ranks = [None] + list(args.ranks)
                for metric in ("read_count", "TPM"):
                    for rk in abund_ranks:
                        tab = build_abundance_table(
                            {n: str(p) for n, p in abund_files.items()},
                            tax=tax, value_column=metric, rank=rk,
                        )
                        suffix = "raw" if rk is None else rk
                        mlabel = "readcount" if metric == "read_count" else "tpm"
                        out_f = Path(args.output) / f"otu_table_{suffix}_{mlabel}.csv"
                        export_otu_table(tab, out_f, transpose=True)
                        otu_tables[f"{suffix}_{mlabel}"] = tab
                log_info(f"✅ Abundance matrices written (read_count + TPM per rank).")
            else:
                log_info(f"\n(no coverage_abundance.tsv found for pattern "
                         f"'{args.abundance_pattern}' — skipping TPM merge)")

        # Summary
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        print(f"Input samples: {len(lca_files)}")
        print(f"Generated OTU tables: {len(otu_tables)}")
        for rank, otu_table in otu_tables.items():
            print(f"  {rank}: {len(otu_table)} samples × {len(otu_table.columns)} taxa")
        print(f"\n✅ OTU tables saved to: {args.output}")

        return 0

    except Exception as e:
        log_info(f"\n❌ Error: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
