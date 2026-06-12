#!/usr/bin/env python3
"""
Virome Classifier CLI - Main classification pipeline

Complete workflow for viral metagenomics classification:
1. Parse paired-end alignment results (R1 + R2)
2. Quality filtering
3. Classification (LCA, Coverage-based, or EM)
4. Kraken report generation

Modes:
  lca       - Simple LCA classification (fast, default)
  coverage  - Coverage-based classification with segment awareness,
              genome-size correction, and real/fake species judgment
  em        - EM algorithm for multi-mapping read resolution
              (iterative abundance estimation)

Usage:
    python -m virome_classifier.cli.classify \\
        --r1 alignments_R1.result \\
        --r2 alignments_R2.result \\
        --taxonomy taxonomy.db \\
        --mask masked_regions.bed \\
        --output results/ \\
        --sample sample_name \\
        --mode coverage
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..core import FilterCriteria, log_info, set_verbose
from ..taxonomy import TaxonomyDB
from ..classification import LCAClassifier
from ..alignment import AlignmentParser, MaskingFilter
from ..reporting import write_all_outputs
from ..coverage_based_classifier3 import CoverageBasedClassifier3, HitQualityFilter, CoverageThresholds


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Virome Classifier - Viral metagenomics classification pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Paired-end mode
  %(prog)s --r1 sample_R1.result --r2 sample_R2.result \\
      --taxonomy tax.db --mask masked.bed --output results/ --sample sample1

  # Single-end mode
  %(prog)s --input sample.result \\
      --taxonomy tax.db --mask masked.bed --output results/ --sample sample1

  # With custom thresholds
  %(prog)s --r1 R1.result --r2 R2.result --taxonomy tax.db \\
      --min-identity 0.9 --min-length 100 --min-coverage 0.25
        """,
    )

    # Input files
    input_group = parser.add_argument_group("Input files")
    input_mode = input_group.add_mutually_exclusive_group(required=True)
    input_mode.add_argument(
        "--input",
        type=str,
        help="Single alignment file (single-end mode)",
    )
    input_mode.add_argument(
        "--r1",
        type=str,
        help="R1 alignment file (paired-end mode, use with --r2)",
    )
    input_group.add_argument(
        "--r2",
        type=str,
        help="R2 alignment file (paired-end mode, requires --r1)",
    )
    input_group.add_argument(
        "--taxonomy",
        "-t",
        type=str,
        required=True,
        help="Taxonomy database (SQLite)",
    )
    input_group.add_argument(
        "--mask",
        "-m",
        type=str,
        required=True,
        help="Masked regions BED file",
    )

    # Output options
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "--output",
        "-o",
        type=str,
        required=True,
        help="Output directory",
    )
    output_group.add_argument(
        "--sample",
        "-s",
        type=str,
        required=True,
        help="Sample name (for output files)",
    )

    # Filtering parameters
    filter_group = parser.add_argument_group("Filtering parameters")
    filter_group.add_argument(
        "--min-identity",
        type=float,
        default=0.85,
        help="Minimum alignment identity (default: 0.85, per ALGORITHM_SELECTION.md sweep)",
    )
    filter_group.add_argument(
        "--min-length",
        type=int,
        default=60,
        help="Minimum alignment length (default: 60, per ALGORITHM_SELECTION.md sweep)",
    )
    filter_group.add_argument(
        "--max-evalue",
        type=float,
        default=1e-3,
        help="Maximum E-value (default: 1e-3)",
    )
    filter_group.add_argument(
        "--min-query-coverage",
        type=float,
        default=0.5,
        help="Minimum query coverage (default: 0.5)",
    )
    filter_group.add_argument(
        "--min-unmasked-coverage",
        type=float,
        default=0.05,
        help="Minimum unmasked genome coverage (default: 0.05, locked 2026-05-28. "
             "Recovers commensal-phage signal that 0.08 truncates while keeping "
             "GT precision; FP control is delegated to --min-rel-abundance "
             "(0.0005) rather than the breadth gate. Sweep history: 0.08 was the "
             "ALGORITHM_SELECTION.md choice; 0.05 is now production.)",
    )
    filter_group.add_argument(
        "--drop-all-masked",
        action="store_true",
        help="Drop reads whose every hit falls (>=50%%) inside a mask region "
             "(all_masked reads = host-mimic/vector/rRNA/MAG contamination), and "
             "remove masked hits of surviving reads, BEFORE LCA. Validated to cut "
             "FP with no KIT TP loss (project_allmasked_drop_validated). Off by "
             "default pending production lock.",
    )
    filter_group.add_argument(
        "--no-breadth-gate",
        action="store_true",
        help="Skip the unmasked-coverage breadth gate entirely (lca mode). Lets "
             "masking act alone so its effect can be measured independently of "
             "breadth. All targets pass; LCA/rel-abundance still apply.",
    )

    # Taxonomy options
    tax_group = parser.add_argument_group("Taxonomy options")
    tax_group.add_argument(
        "--root-taxid",
        type=int,
        default=10239,
        help="Root taxonomy ID (default: 10239 for Viruses)",
    )
    tax_group.add_argument(
        "--virus-root",
        action="store_true",
        default=True,
        help="Use Viruses (10239) as root in reports (default: True)",
    )

    # Classification mode
    mode_group = parser.add_argument_group("Classification mode")
    mode_group.add_argument(
        "--mode",
        type=str,
        choices=["lca", "lca_conditional", "coverage", "em", "ml_filter"],
        default="lca",
        help="Classification mode: lca (LCA with the mode-agnostic realness "
             "post-filter; recommended default per ALGORITHM_SELECTION.md — best "
             "detection F1/FP and zero false-positive abundance under the "
             "expected-genome-length breadth denominator), coverage (coverage-based "
             "with segment awareness), em (EM algorithm for multi-mapping). "
             "Default: lca",
    )
    mode_group.add_argument(
        "--classification-rank",
        type=str,
        choices=["strain", "species", "genus", "family"],
        default="strain",
        help="Taxonomic rank for coverage/EM classification (default: strain)",
    )
    mode_group.add_argument(
        "--read-assign",
        type=str,
        choices=["best_hit", "lca"],
        default="best_hit",
        help="STAGE-0 read-level taxonomy in lca mode: 'best_hit' (default, user "
             "core / Strategy B) assigns each read its single max-bitscore hit's "
             "taxid — one taxon per read, lowest FP; 'lca' assigns the common "
             "ancestor of all the read's hits (cross-map reads rise to a higher "
             "rank). Both feed the SAME downstream rank roll-up (--lca-fix-rank), "
             "phage-host roll-up, and FP post-filter. Default: best_hit",
    )
    mode_group.add_argument(
        "--lca-fix-rank",
        type=str,
        choices=["none", "species", "genus"],
        default="none",
        help="lca mode: lift each read's natural LCA UP to this rank "
             "(get_taxid_at_rank). 'none' (default) keeps the natural common "
             "ancestor (varies per read; ~species for most). 'species' ≈ natural. "
             "'genus' merges same-genus relatives into one taxon (loses species "
             "resolution and dissolves the genus-competition gate). Reads whose "
             "natural LCA is already above the target rank are left unchanged.",
    )
    mode_group.add_argument(
        "--multi-mapping-mode",
        type=str,
        choices=["all", "best_hit", "local_depth"],
        default="best_hit",
        help="Multi-mapping read handling in coverage mode. best_hit gives the "
             "lowest false-positive rate (ALGORITHM_SELECTION.md); 'all' is the "
             "former default. Default: best_hit",
    )
    mode_group.add_argument(
        "--breadth-denominator",
        type=str,
        choices=["expected", "aligned"],
        default="expected",
        help="Genome-coverage breadth denominator (coverage mode). 'expected' "
             "divides by the taxon's whole-genome length so reads on only the "
             "short segment of a multipartite virus do not inflate breadth "
             "(removes partial/short-reference FPs; needs DB column "
             "expected_genome_length). 'aligned' is the former behaviour "
             "(hit segments only). Default: expected",
    )
    mode_group.add_argument(
        "--use-depth-entropy",
        action="store_true",
        default=False,
        help="Use depth entropy for masking-aware coverage (coverage mode)",
    )
    # --- FP-leak fix toggles (coverage mode; off by default) ---
    mode_group.add_argument(
        "--min-species-confidence",
        type=float,
        default=0.5,
        help="lca_conditional mode: species with unique-read fraction below this "
             "are retreated to genus rank (hierarchical LCA). Default 0.5.",
    )
    mode_group.add_argument(
        "--min-rel-abundance",
        type=float,
        default=0.0005,
        help="DEFAULT false-positive control. Drop any taxon assigned < this "
             "fraction of all classified reads (default 0.0005 = 0.05%%, "
             "ALGORITHM_SELECTION.md sweep). Removes same-genus cross-mapping "
             "shadows (true species sit >=0.06%% vs shadows ~0.001%%) and the "
             "low-abundance contaminant tail, depth-robustly; 0 = off.",
    )
    mode_group.add_argument(
        "--min-unique-fraction",
        type=float,
        default=0.0,
        help="Require unique/total reads >= this for a taxon to be real "
             "(0 = off, default). Optional alternative gate; the rel-abundance "
             "cut is the default FP control.",
    )
    mode_group.add_argument(
        "--min-read-count",
        type=int,
        default=0,
        help="Absolute per-taxon read floor in the mode-agnostic FP post-filter: "
             "drop any taxon supported by fewer than this many classified reads "
             "(applies to lca/em/coverage alike). Depth-DEPENDENT, so 0 = off by "
             "default; the LOCKED policy value is 10 (pass --min-read-count 10). "
             "Complements the depth-robust --min-rel-abundance gate.",
    )
    mode_group.add_argument(
        "--phage-host-rollup",
        action="store_true",
        help="After the FP post-filter, roll phage detections up to their "
             "bacterial/archaeal HOST GENUS (phages_of_<Host>), collapsing "
             "phage cross-map dispersion into an interpretable per-host read-out. "
             "Only phage with a KNOWN host are merged; host-unknown phage stay "
             "per-species and NON-phage (herpes/human/...) stay at species. Host "
             "is from refseq_metadata.host (+ title fill), VMR-independent. Off by "
             "default. See docs/phage_host_processing.md.",
    )
    mode_group.add_argument(
        "--ml-model",
        type=str,
        default=None,
        help="Trained ML filter model (.joblib) for --mode ml_filter.",
    )
    mode_group.add_argument(
        "--ml-threshold",
        type=float,
        default=0.5,
        help="TP probability cutoff for --mode ml_filter (default: 0.5).",
    )
    mode_group.add_argument(
        "--em-iterations",
        type=int,
        default=20,
        help="Number of EM iterations (em mode, default: 20)",
    )
    mode_group.add_argument(
        "--em-convergence",
        type=float,
        default=1e-6,
        help="EM convergence threshold (default: 1e-6)",
    )

    # Other options
    other_group = parser.add_argument_group("Other options")
    other_group.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    other_group.add_argument(
        "--no-kraken",
        action="store_true",
        help="Skip Kraken report generation",
    )

    args = parser.parse_args()

    # Validate paired-end mode
    if args.r1 and not args.r2:
        parser.error("--r1 requires --r2 (paired-end mode)")
    if args.r2 and not args.r1:
        parser.error("--r2 requires --r1 (paired-end mode)")

    return args


def load_taxonomy(db_path: str, root_taxid: int, verbose: bool) -> TaxonomyDB:
    """Load taxonomy database."""
    log_info(f"📚 Loading taxonomy from {db_path} (root={root_taxid})...")
    tax = TaxonomyDB.from_sqlite(db_path, root_taxid=root_taxid)
    log_info(f"✅ Loaded {len(tax):,} taxa")
    return tax


def parse_alignments(
    r1_file: Optional[str],
    r2_file: Optional[str],
    input_file: Optional[str],
    verbose: bool,
) -> pd.DataFrame:
    """Parse alignment files (paired-end or single-end)."""
    parser = AlignmentParser(normalize_headers=True)

    if r1_file and r2_file:
        # Paired-end mode
        log_info("\n📄 Parsing paired-end alignment results...")
        r1_hits = parser.parse(Path(r1_file))
        r1_hits["strand"] = "+"
        log_info(f"  R1: {len(r1_hits):,} hits")

        r2_hits = parser.parse(Path(r2_file))
        r2_hits["strand"] = "-"
        log_info(f"  R2: {len(r2_hits):,} hits")

        combined_hits = pd.concat([r1_hits, r2_hits], ignore_index=True)
        log_info(f"  Combined: {len(combined_hits):,} hits")
        return combined_hits
    else:
        # Single-end mode
        log_info(f"\n📄 Parsing alignment results from {input_file}...")
        hits = parser.parse(Path(input_file))
        log_info(f"  Loaded: {len(hits):,} hits")
        return hits


def apply_quality_filter(
    hits: pd.DataFrame,
    min_identity: float,
    min_length: int,
    max_evalue: float,
    min_query_coverage: float,
    verbose: bool,
) -> pd.DataFrame:
    """Apply quality filtering."""
    log_info("\n🔍 Applying quality filters...")
    log_info(f"  Identity ≥ {min_identity}")
    log_info(f"  Length ≥ {min_length}")
    log_info(f"  E-value ≤ {max_evalue}")
    log_info(f"  Query coverage ≥ {min_query_coverage}")

    parser = AlignmentParser()
    criteria = FilterCriteria(
        min_identity=min_identity,
        min_alignment_length=min_length,
        max_evalue=max_evalue,
        min_query_coverage=min_query_coverage,
    )

    filtered = parser.filter(hits, criteria)
    retention = len(filtered) / len(hits) * 100 if len(hits) > 0 else 0
    log_info(f"✅ Filtered: {len(filtered):,}/{len(hits):,} hits ({retention:.1f}% retention)")

    return filtered


def build_masking_filter(mask_file: str) -> "MaskingFilter":
    """Construct a MaskingFilter from a BED file (thin factory wrapper)."""
    return MaskingFilter.from_bed_file(mask_file)


def apply_breadth_gate(
    masking_filter: "MaskingFilter",
    hits: pd.DataFrame,
    min_coverage: float,
):
    """Run the unmasked-coverage breadth gate (thin wrapper + logging)."""
    log_info(f"\n🎭 Applying masking filter (min coverage: {min_coverage})...")
    result = masking_filter.filter_by_unmasked_coverage(hits, min_coverage=min_coverage)
    log_info(result.summary())
    log_info(f"✅ Passed: {len(result.passed):,} hits from {result.n_passed_targets} targets")
    return result


def perform_lca_classification(
    hits: pd.DataFrame,
    tax: TaxonomyDB,
    verbose: bool,
    fix_rank: str = "none",
) -> pd.DataFrame:
    """Thin wrapper around LCAClassifier.classify (logic lives in the class)."""
    return LCAClassifier(tax, verbose).classify(hits, fix_rank=fix_rank)


def perform_conditional_lca_classification(
    hits: pd.DataFrame,
    tax: TaxonomyDB,
    verbose: bool,
    min_species_confidence: float = 0.5,
) -> pd.DataFrame:
    """Thin wrapper around LCAClassifier.classify_conditional."""
    return LCAClassifier(tax, verbose).classify_conditional(
        hits, min_species_confidence=min_species_confidence)


def export_results(
    lca_df: pd.DataFrame,
    final_hits: pd.DataFrame,
    stats_df: pd.DataFrame,
    output_dir: Path,
    sample_name: str,
    verbose: bool,
) -> None:
    """Export classification results to CSV files."""
    log_info("\n💾 Exporting classification results...")

    output_dir.mkdir(parents=True, exist_ok=True)

    # Save LCA results
    lca_file = output_dir / f"{sample_name}_read_classification.csv"
    lca_df.to_csv(lca_file, index=False)
    log_info(f"  LCA results: {lca_file}")

    # Save filtered hits
    hits_file = output_dir / f"{sample_name}_final_hits.csv"
    final_hits.to_csv(hits_file, index=False)
    log_info(f"  Final hits: {hits_file}")

    # Save coverage statistics
    stats_file = output_dir / f"{sample_name}_coverage_stats.csv"
    stats_df.to_csv(stats_file)
    log_info(f"  Coverage stats: {stats_file}")


def generate_kraken_reports(
    lca_df: pd.DataFrame,
    tax: TaxonomyDB,
    output_dir: Path,
    sample_name: str,
    virus_root: bool,
    verbose: bool,
) -> None:
    """Generate Kraken-format reports."""
    log_info("\n📊 Generating Kraken reports...")

    kraken_files = write_all_outputs(
        results_df=lca_df,
        tax=tax,
        output_dir=str(output_dir),
        sample_name=sample_name,
        virus_root=virus_root,
    )

    log_info("✅ Kraken outputs generated:")
    for file_type, file_path in kraken_files.items():
        log_info(f"  {file_type}: {file_path}")


def run_coverage_mode(args, tax, combined_hits, output_dir):
    """Run coverage-based classification pipeline."""
    from ..coverage_based_classifier3 import CoverageBasedClassifier3, HitQualityFilter, CoverageThresholds

    log_info(f"\n🧬 Running COVERAGE-BASED classification (rank={args.classification_rank})...")

    # Load segment info from DB
    segment_info = tax.get_segment_info()

    classifier = CoverageBasedClassifier3(
        taxonomy_db=tax,
        segment_info=segment_info,
        hit_filter=HitQualityFilter(
            min_identity=args.min_identity,
            min_aligned_length=args.min_length,
            min_query_coverage=args.min_query_coverage,
        ),
        thresholds=CoverageThresholds(
            min_unmasked_coverage=args.min_unmasked_coverage,
            min_unique_fraction=getattr(args, "min_unique_fraction", 0.0),
        ),
        mask_bed_file=args.mask,
        multi_mapping_mode=args.multi_mapping_mode,
        classification_rank=args.classification_rank,
        use_depth_entropy=args.use_depth_entropy,
        breadth_denominator=getattr(args, "breadth_denominator", "expected"),
        verbose=args.verbose,
    )

    assignment_df, abundance_df, taxon_cov = classifier.classify(combined_hits)

    # Save coverage-specific outputs
    output_dir.mkdir(parents=True, exist_ok=True)
    abundance_df.to_csv(output_dir / f"{args.sample}_coverage_abundance.tsv", sep='\t', index=False)
    log_info(f"  Coverage abundance: {output_dir / f'{args.sample}_coverage_abundance.tsv'}")

    # Save summary (+ vectorized paired-end support feature: both_frac)
    summary_df = classifier.get_summary_dataframe(taxon_cov)
    try:
        from ..alignment.pairing import species_pair_support
        if "strand" in combined_hits.columns and "taxon_taxid" in summary_df.columns:
            ps = species_pair_support(combined_hits, tax)  # species, reads, both_frac
            summary_df = summary_df.merge(
                ps[["species", "both_frac"]].rename(columns={"species": "taxon_taxid"}),
                on="taxon_taxid", how="left",
            )
            summary_df["both_frac"] = summary_df["both_frac"].fillna(0.0)
    except Exception as e:
        log_info(f"  (pairing feature skipped: {e})")
    summary_df.to_csv(output_dir / f"{args.sample}_coverage_summary.tsv", sep='\t', index=False)
    log_info(f"  Coverage summary: {output_dir / f'{args.sample}_coverage_summary.tsv'}")
    # Stash for ml_filter mode to reuse in-memory (avoids re-reading the TSV).
    run_coverage_mode._last_summary = summary_df

    # Convert to LCA-compatible format for Kraken reporting
    lca_df = _coverage_to_lca_format(assignment_df, tax)

    return lca_df, assignment_df, abundance_df, taxon_cov


def run_ml_filter_mode(args, tax, combined_hits, output_dir):
    """Coverage detection, then a trained ML model re-judges each detected species
    TP vs FP (CellTypist-style). Reads assigned to ML-rejected species are dropped."""
    from ..classification.ml_filter import score

    log_info("\n🧬 Running ML_FILTER (coverage detection + learned TP/FP refinement)...")
    lca_df, assignment_df, abundance_df, taxon_cov = run_coverage_mode(args, tax, combined_hits, output_dir)

    model_path = getattr(args, "ml_model", None)
    if not model_path:
        log_info("  ⚠️  --ml-model not given; returning coverage result unchanged.")
        return lca_df, assignment_df, abundance_df, taxon_cov

    # Build the feature summary in-memory from taxon_cov (no disk round-trip).
    summary_df = run_coverage_mode._last_summary if hasattr(run_coverage_mode, "_last_summary") else None
    if summary_df is None or summary_df.empty:
        log_info("  ⚠️  no coverage summary to score; returning coverage result.")
        return lca_df, assignment_df, abundance_df, taxon_cov

    scored = score(summary_df, model_path, threshold=getattr(args, "ml_threshold", 0.5))
    scored.to_csv(output_dir / f"{args.sample}_ml_scored.tsv", sep='\t', index=False)
    keep = set(scored.loc[scored["ml_real"], "taxon_taxid"].astype(int))
    rejected = set(scored.loc[~scored["ml_real"], "taxon_taxid"].astype(int))
    log_info(f"  ML kept {len(keep)} species, rejected {len(rejected)} (FP) "
             f"@ threshold {getattr(args, 'ml_threshold', 0.5)}")

    # Filter assignment_df to ML-kept species, rebuild lca_df
    if assignment_df is not None and not assignment_df.empty:
        assignment_df = assignment_df[assignment_df["working_taxid"].astype(int).isin(keep)].copy()
    lca_df = _coverage_to_lca_format(assignment_df, tax)
    if abundance_df is not None and "taxon_taxid" in abundance_df.columns:
        abundance_df = abundance_df[abundance_df["taxon_taxid"].astype(int).isin(keep)].copy()
    return lca_df, assignment_df, abundance_df, taxon_cov


def _coverage_to_lca_format(assignment_df, tax):
    """Convert coverage assignment to LCA-compatible DataFrame for Kraken reporting."""
    if assignment_df is None or assignment_df.empty:
        return pd.DataFrame()

    rows = []
    for _, row in assignment_df.iterrows():
        taxid = int(row['working_taxid'])
        rows.append({
            'query': row['query'],
            'lca_taxid': taxid,
            'lca_name': tax.get_name(taxid) or f"Unknown ({taxid})",
            'lca_rank': tax.get_rank(taxid) or 'no rank',
            'qlen': 100,
            'read_count': 1,
            'n_hits': 1,
            'n_unique_taxids': 1,
            'all_taxids': str(taxid),
        })

    return pd.DataFrame(rows)


def run_em_mode(args, tax, combined_hits, output_dir):
    """Run EM algorithm classification pipeline."""
    from ..classification.em_classifier import EMClassifier

    log_info(f"\n🧬 Running EM classification (iterations={args.em_iterations})...")

    # Load segment info
    segment_info = tax.get_segment_info()

    classifier = EMClassifier(
        taxonomy_db=tax,
        segment_info=segment_info,
        min_identity=args.min_identity,
        min_length=args.min_length,
        max_evalue=args.max_evalue,
        min_query_coverage=args.min_query_coverage,
        mask_bed_file=args.mask,
        max_iterations=args.em_iterations,
        convergence_threshold=args.em_convergence,
        verbose=args.verbose,
    )

    lca_df, abundance_df = classifier.classify(combined_hits)

    # Save EM-specific outputs
    output_dir.mkdir(parents=True, exist_ok=True)
    abundance_df.to_csv(output_dir / f"{args.sample}_em_abundance.tsv", sep='\t', index=False)
    log_info(f"  EM abundance: {output_dir / f'{args.sample}_em_abundance.tsv'}")

    return lca_df, abundance_df


def run_lca_mode(args, tax, combined_hits, output_dir):
    """LCA pipeline: quality filter -> [masking: classify -> drop_all_masked]
    -> breadth gate -> LCA -> export.

    Masking (classification + optional all_masked drop) is now separated from
    the breadth gate so masking can act independently. With --drop-all-masked
    OFF and the breadth gate ON (defaults), the path is byte-identical to the
    previous behaviour: classify_hits adds a label column that is dropped before
    the gate, and filter_by_unmasked_coverage runs unchanged on the same hits.
    """
    filtered_hits = apply_quality_filter(
        combined_hits, args.min_identity, args.min_length,
        args.max_evalue, args.min_query_coverage, args.verbose,
    )

    mf = build_masking_filter(args.mask) if args.mask else None

    # ── masking stage (independent of breadth gate) ──
    hits = filtered_hits
    if mf is not None and getattr(args, "drop_all_masked", False):
        n_in, reads_in = len(hits), hits["query"].nunique()
        hits = mf.drop_all_masked_reads(hits)
        log_info(f"🧽 all_masked drop: {n_in:,}→{len(hits):,} hits, "
                 f"{reads_in:,}→{hits['query'].nunique():,} reads")

    # ── breadth gate (optional) ──
    if mf is not None and not getattr(args, "no_breadth_gate", False):
        result = apply_breadth_gate(mf, hits, args.min_unmasked_coverage)
        passed_hits, stats_df = result.passed, result.stats
    else:
        # breadth gate skipped (or no mask): masking-only / unfiltered passthrough
        passed_hits = hits
        stats_df = mf.calculate_stats(hits) if mf is not None else pd.DataFrame()

    if getattr(args, "mode", "lca") == "lca_conditional":
        lca_df = perform_conditional_lca_classification(
            passed_hits, tax, args.verbose,
            min_species_confidence=getattr(args, "min_species_confidence", 0.5),
        )
    elif getattr(args, "read_assign", "best_hit") == "best_hit":
        # STAGE-0: best-hit (user core / Strategy B) — one taxon per read
        from ..classification.besthit_classifier import BestHitClassifier
        lca_df = BestHitClassifier(tax, args.verbose).classify(
            passed_hits, fix_rank=getattr(args, "lca_fix_rank", "none"))
    else:
        lca_df = perform_lca_classification(
            passed_hits, tax, args.verbose,
            fix_rank=getattr(args, "lca_fix_rank", "none"),
        )
    export_results(lca_df, passed_hits, stats_df, output_dir, args.sample, args.verbose)
    return lca_df


# Mode registry: name -> callable(args, tax, combined_hits, output_dir) -> lca_df.
# The coverage/em/ml_filter runners return richer tuples (lca_df first); the
# adapters keep the dispatch uniform. Register a new mode here — no if/elif edits.
MODE_DISPATCH = {
    "lca":             run_lca_mode,
    "lca_conditional": run_lca_mode,   # same pipeline; mode flag triggers retreat
    "coverage":        lambda a, t, h, o: run_coverage_mode(a, t, h, o)[0],
    "em":              lambda a, t, h, o: run_em_mode(a, t, h, o)[0],
    "ml_filter":       lambda a, t, h, o: run_ml_filter_mode(a, t, h, o)[0],
}


def main():
    """Main entry point."""
    args = parse_args()

    # Set verbosity
    if args.verbose:
        set_verbose(True)

    # Print header
    print("=" * 70)
    print(f"VIROME CLASSIFIER - {args.mode.upper()} MODE")
    print("=" * 70)

    try:
        # 1. Load taxonomy
        tax = load_taxonomy(args.taxonomy, args.root_taxid, args.verbose)

        # 2. Parse alignments
        combined_hits = parse_alignments(
            args.r1,
            args.r2,
            args.input,
            args.verbose,
        )

        output_dir = Path(args.output)

        # ====== MODE DISPATCH (registry) ======
        # Each mode is a function (args, tax, combined_hits, output_dir) -> lca_df.
        # Adding a mode = register it in MODE_DISPATCH; no new if/elif needed.
        if args.mode not in MODE_DISPATCH:
            raise ValueError(f"Unknown mode: {args.mode}")
        lca_df = MODE_DISPATCH[args.mode](args, tax, combined_hits, output_dir)

        # ====== COMMON: mode-agnostic FP post-filter ======
        # Relative-abundance cut (default FP control) + optional unique-fraction,
        # applied uniformly to lca / em / coverage so the three modes are compared
        # under the SAME FP control (only removes taxa; harmless if a mode already
        # filtered). ml_filter uses its own trained model, so skip it.
        if args.mode in ("lca", "em", "coverage") and len(lca_df) > 0 and (
            getattr(args, "min_rel_abundance", 0.0) > 0
            or getattr(args, "min_unique_fraction", 0.0) > 0
            or getattr(args, "min_read_count", 0) > 0
        ):
            from ..classification.fp_postfilter import apply_fp_postfilter
            lca_df = apply_fp_postfilter(
                lca_df, tax,
                min_rel_abundance=getattr(args, "min_rel_abundance", 0.0),
                min_unique_fraction=getattr(args, "min_unique_fraction", 0.0),
                min_read_count=getattr(args, "min_read_count", 0),
            )

        # ====== COMMON: optional phage -> host-genus roll-up ======
        # Collapses phage cross-map dispersion into per-host taxa; non-phage and
        # host-unknown phage are left untouched. Applied to every counting mode.
        if getattr(args, "phage_host_rollup", False) and len(lca_df) > 0:
            from ..classification.phage_host_rollup import (
                build_phage_host_map, apply_phage_host_rollup)
            tax2host, phage_set = build_phage_host_map(args.taxonomy)
            lca_df = apply_phage_host_rollup(lca_df, tax2host, phage_set)

        # ====== COMMON: per-query classification CSV ======
        # LCA mode writes this via export_results(); coverage/em modes write it
        # here so every mode emits a consistent {sample}_read_classification.csv
        # (consumed by the Nextflow module and OTU merge).
        if args.mode in ("coverage", "em") and len(lca_df) > 0:
            output_dir.mkdir(parents=True, exist_ok=True)
            lca_df.to_csv(output_dir / f"{args.sample}_read_classification.csv", index=False)
            log_info(f"  Classification: {output_dir / f'{args.sample}_read_classification.csv'}")

        # ====== COMMON: Kraken reports ======
        # Always emit a .kreport, even for a zero-read / fully-filtered sample, so
        # a downstream per-sample merge does not fail on a missing file. The
        # writer (write_kraken_output) handles an empty lca_df by writing an
        # empty report rather than crashing.
        if not args.no_kraken:
            generate_kraken_reports(
                lca_df,
                tax,
                output_dir,
                args.sample,
                args.virus_root,
                args.verbose,
            )

        # Summary
        print("\n" + "=" * 70)
        print("PIPELINE SUMMARY")
        print("=" * 70)
        print(f"Mode: {args.mode}")
        print(f"Input hits: {len(combined_hits):,}")
        # (per-mode quality/masking counts are logged inside each runner via
        # log_info; they are not in main()'s scope, so don't reference them here.)
        print(f"Classified: {len(lca_df):,} queries")
        if len(lca_df) > 0:
            print(f"Unique taxa: {lca_df['lca_taxid'].nunique()}")
        print(f"\n✅ Pipeline complete! Results saved to: {output_dir}")

        return 0

    except Exception as e:
        log_info(f"\n❌ Error: {e}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())
