"""
OTU table generation from multiple sample LCA results.

This module provides functions to merge multiple LCA classification results
into a unified OTU/ASV table (sample × taxon matrix).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional, Union
import pandas as pd

from ..core import log_info, log_verbose
from ..taxonomy import TaxonomyDB


def partition_nonempty_samples(
    lca_files: Dict[str, Union[str, Path]],
    count_column: str = "read_count",
) -> tuple[Dict[str, Union[str, Path]], List[str]]:
    """Split sample files into (kept, dropped) by whether they have any
    classified reads.

    A sample is treated as EMPTY (0 identified reads) when its per-query
    classification CSV has no data rows, or its read_count sum is 0, or it lacks
    the 'lca_taxid' column entirely (e.g. a header-only file written for a
    fully-filtered sample). Empty samples are excluded from the OTU merge so a
    single 0-read sample cannot break the whole table; the dropped names are
    returned so the caller can log them.
    """
    kept: Dict[str, Union[str, Path]] = {}
    dropped: List[str] = []
    for sample_name, file_path in lca_files.items():
        try:
            df = pd.read_csv(file_path)
        except Exception:
            dropped.append(sample_name)
            continue
        if "lca_taxid" not in df.columns or len(df) == 0:
            dropped.append(sample_name)
            continue
        if count_column in df.columns and pd.to_numeric(
            df[count_column], errors="coerce").fillna(0).sum() <= 0:
            dropped.append(sample_name)
            continue
        kept[sample_name] = file_path
    return kept, dropped


def load_lca_results(file_path: Union[str, Path]) -> pd.DataFrame:
    """
    Load LCA classification results from CSV file.

    Args:
        file_path: Path to LCA classification CSV file

    Returns:
        DataFrame with LCA results
    """
    df = pd.read_csv(file_path)
    required_cols = ["query", "lca_taxid"]

    if not all(col in df.columns for col in required_cols):
        raise ValueError(f"LCA file must contain columns: {required_cols}")

    return df


def build_otu_table(
    lca_files: Dict[str, Union[str, Path]],
    tax: Optional[TaxonomyDB] = None,
    count_column: str = "read_count",
    normalize: bool = False,
) -> pd.DataFrame:
    """
    Build OTU table from multiple LCA classification files.

    Creates a sample × taxon matrix where:
    - Rows: Samples
    - Columns: Taxa (taxids)
    - Values: Read counts (or normalized abundances)

    Args:
        lca_files: Dictionary mapping sample names to LCA CSV file paths
        tax: Optional TaxonomyDB for adding taxonomy metadata
        count_column: Column name containing read counts (default: "read_count")
        normalize: If True, normalize counts to relative abundances per sample

    Returns:
        OTU table DataFrame (samples as rows, taxa as columns)

    Example:
        >>> lca_files = {
        ...     "sample1": "sample1_read_classification.csv",
        ...     "sample2": "sample2_read_classification.csv",
        ... }
        >>> otu = build_otu_table(lca_files)
    """
    log_info(f"📊 Building OTU table from {len(lca_files)} samples...")

    # Load all LCA results
    all_data = []
    for sample_name, file_path in lca_files.items():
        log_verbose(f"  Loading {sample_name}...")
        df = load_lca_results(file_path)

        # Count reads per taxon
        if count_column in df.columns:
            taxon_counts = df.groupby("lca_taxid")[count_column].sum()
        else:
            # Count queries if no read_count column
            taxon_counts = df["lca_taxid"].value_counts()

        # Add sample name
        taxon_counts.name = sample_name
        all_data.append(taxon_counts)

    # Merge into OTU table
    otu_table = pd.DataFrame(all_data).fillna(0).astype(int)
    otu_table.index.name = "sample"

    log_info(f"✅ OTU table: {len(otu_table)} samples × {len(otu_table.columns)} taxa")

    # Normalize if requested
    if normalize:
        log_verbose("  Normalizing to relative abundances...")
        otu_table = otu_table.div(otu_table.sum(axis=1), axis=0)

    # Add taxonomy metadata if provided
    if tax is not None:
        log_verbose("  Adding taxonomy metadata...")
        otu_table = add_taxonomy_metadata(otu_table, tax)

    return otu_table


def add_taxonomy_metadata(
    otu_table: pd.DataFrame,
    tax: TaxonomyDB,
) -> pd.DataFrame:
    """
    Add taxonomy metadata columns to OTU table.

    Converts columns from taxids to multi-level columns with:
    - taxid: Taxonomy ID
    - name: Taxon name
    - rank: Taxonomic rank

    Args:
        otu_table: OTU table with taxids as columns
        tax: TaxonomyDB instance

    Returns:
        OTU table with taxonomy metadata in column names
    """
    # Create metadata for each taxid
    metadata = []
    for taxid in otu_table.columns:
        metadata.append(
            {
                "taxid": int(taxid),
                "name": tax.get_name(int(taxid)) or f"taxid_{taxid}",
                "rank": tax.get_rank(int(taxid)) or "no rank",
            }
        )

    # Create multi-index columns
    metadata_df = pd.DataFrame(metadata)
    otu_table.columns = pd.MultiIndex.from_frame(metadata_df[["taxid", "name", "rank"]])

    return otu_table


def build_otu_table_at_rank(
    lca_files: Dict[str, Union[str, Path]],
    tax: TaxonomyDB,
    rank: str = "genus",
    count_column: str = "read_count",
    normalize: bool = False,
) -> pd.DataFrame:
    """
    Build OTU table collapsed at specific taxonomic rank.

    Instead of using LCA taxids directly, aggregate counts at a specific rank
    (e.g., genus, species) by traversing lineages.

    Args:
        lca_files: Dictionary mapping sample names to LCA CSV file paths
        tax: TaxonomyDB instance
        rank: Taxonomic rank to collapse to (e.g., "genus", "species", "family")
        count_column: Column name containing read counts
        normalize: If True, normalize to relative abundances

    Returns:
        OTU table at specified rank

    Example:
        >>> # Genus-level OTU table
        >>> otu_genus = build_otu_table_at_rank(
        ...     lca_files, tax, rank="genus"
        ... )
    """
    log_info(f"📊 Building {rank}-level OTU table from {len(lca_files)} samples...")

    all_data = []
    for sample_name, file_path in lca_files.items():
        log_verbose(f"  Processing {sample_name}...")
        df = load_lca_results(file_path)

        # Get counts
        if count_column in df.columns:
            counts = df.set_index("lca_taxid")[count_column]
        else:
            counts = df["lca_taxid"].value_counts()

        # Map taxids to target rank
        rank_counts = {}
        for taxid, count in counts.items():
            # Get lineage and find target rank
            lineage = tax.get_lineage(int(taxid))
            target_taxid = None

            for ancestor in reversed(lineage):  # Start from most specific
                if tax.get_rank(ancestor) == rank:
                    target_taxid = ancestor
                    break

            if target_taxid:
                rank_counts[target_taxid] = rank_counts.get(target_taxid, 0) + int(count)

        # Convert to Series
        rank_series = pd.Series(rank_counts, name=sample_name)
        all_data.append(rank_series)

    # Merge into OTU table
    otu_table = pd.DataFrame(all_data).fillna(0).astype(int)
    otu_table.index.name = "sample"

    log_info(f"✅ {rank}-level OTU table: {len(otu_table)} samples × {len(otu_table.columns)} {rank}s")

    # Normalize if requested
    if normalize:
        log_verbose("  Normalizing to relative abundances...")
        otu_table = otu_table.div(otu_table.sum(axis=1), axis=0)

    # Add taxonomy names
    log_verbose("  Adding taxonomy names...")
    taxid_to_name = {taxid: tax.get_name(int(taxid)) or f"taxid_{taxid}" for taxid in otu_table.columns}
    otu_table.columns = pd.Index([f"{taxid_to_name[tid]} ({tid})" for tid in otu_table.columns])

    return otu_table


def filter_otu_table(
    otu_table: pd.DataFrame,
    min_count: int = 0,
    min_samples: int = 1,
    min_prevalence: float = 0.0,
) -> pd.DataFrame:
    """
    Filter OTU table by count thresholds.

    Args:
        otu_table: OTU table to filter
        min_count: Minimum total count across all samples
        min_samples: Minimum number of samples where taxon must be present
        min_prevalence: Minimum prevalence (fraction of samples, 0-1)

    Returns:
        Filtered OTU table

    Example:
        >>> # Keep taxa present in at least 2 samples with ≥10 total reads
        >>> filtered = filter_otu_table(otu, min_count=10, min_samples=2)
    """
    log_verbose(f"🔍 Filtering OTU table (min_count={min_count}, min_samples={min_samples})...")

    original_taxa = len(otu_table.columns)

    # Filter by total count
    if min_count > 0:
        taxon_sums = otu_table.sum(axis=0)
        otu_table = otu_table.loc[:, taxon_sums >= min_count]

    # Filter by number of samples
    if min_samples > 1:
        presence = (otu_table > 0).sum(axis=0)
        otu_table = otu_table.loc[:, presence >= min_samples]

    # Filter by prevalence
    if min_prevalence > 0:
        prevalence = (otu_table > 0).sum(axis=0) / len(otu_table)
        otu_table = otu_table.loc[:, prevalence >= min_prevalence]

    remaining_taxa = len(otu_table.columns)
    pct = (remaining_taxa / original_taxa * 100) if original_taxa else 0.0
    log_info(f"✅ Filtered: {original_taxa} → {remaining_taxa} taxa ({pct:.1f}%)")

    return otu_table


def export_otu_table(
    otu_table: pd.DataFrame,
    output_file: Union[str, Path],
    format: str = "csv",
    transpose: bool = False,
) -> None:
    """
    Export OTU table to file.

    Args:
        otu_table: OTU table to export
        output_file: Output file path
        format: Output format ('csv', 'tsv', or 'biom')
        transpose: If True, export with taxa as rows (standard format)

    Example:
        >>> export_otu_table(otu, "otu_table.csv", transpose=True)
    """
    log_info(f"💾 Exporting OTU table to {output_file}...")

    output_path = Path(output_file)

    # Transpose if requested (taxa as rows is more standard)
    if transpose:
        otu_table = otu_table.T

    # Export based on format
    if format == "csv":
        otu_table.to_csv(output_path)
    elif format == "tsv":
        otu_table.to_csv(output_path, sep="\t")
    elif format == "biom":
        # BIOM format requires biom-format package
        try:
            import biom

            # Convert to BIOM format
            table = biom.Table(
                otu_table.values.T if not transpose else otu_table.values,
                observation_ids=otu_table.columns if not transpose else otu_table.index,
                sample_ids=otu_table.index if not transpose else otu_table.columns,
            )
            with biom.util.biom_open(str(output_path), "w") as f:
                table.to_hdf5(f, "OTU table generated by virome_classifier")
        except ImportError:
            log_info("⚠️  biom-format not installed, falling back to TSV")
            otu_table.to_csv(output_path, sep="\t")
    else:
        raise ValueError(f"Unknown format: {format}")

    log_info(f"✅ Exported OTU table: {len(otu_table)} × {len(otu_table.columns)}")


def create_otu_pipeline(
    lca_files: Dict[str, Union[str, Path]],
    tax: TaxonomyDB,
    output_dir: Union[str, Path],
    *,
    ranks: Optional[List[str]] = None,
    min_count: int = 10,
    min_samples: int = 1,
    normalize: bool = False,
) -> Dict[str, pd.DataFrame]:
    """
    Complete OTU table generation pipeline.

    Creates filtered OTU tables at multiple taxonomic ranks.

    Args:
        lca_files: Dictionary mapping sample names to LCA file paths
        tax: TaxonomyDB instance
        output_dir: Output directory
        ranks: List of ranks to generate tables for (default: genus, species)
        min_count: Minimum total count for filtering
        min_samples: Minimum number of samples for filtering
        normalize: Normalize to relative abundances

    Returns:
        Dictionary mapping rank names to OTU tables

    Example:
        >>> lca_files = {
        ...     "sample1": "sample1_lca.csv",
        ...     "sample2": "sample2_lca.csv",
        ... }
        >>> tables = create_otu_pipeline(
        ...     lca_files, tax, "./otu_results",
        ...     ranks=["genus", "species", "family"]
        ... )
    """
    if ranks is None:
        ranks = ["genus", "species"]

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    log_info("=" * 70)
    log_info("OTU TABLE GENERATION PIPELINE")
    log_info("=" * 70)

    # Drop samples with 0 identified reads BEFORE merging, so one empty sample
    # cannot break the whole table. Report which ones were dropped.
    n_total = len(lca_files)
    lca_files, empty_samples = partition_nonempty_samples(lca_files)
    if empty_samples:
        log_info(
            f"\n⚠️  Excluded {len(empty_samples)}/{n_total} sample(s) with 0 "
            f"identified reads from the OTU merge: {', '.join(sorted(empty_samples))}"
        )
        # also persist the list so it survives in the run outputs
        try:
            with open(output_path / "empty_samples.txt", "w") as fh:
                fh.write("# samples excluded from OTU merge (0 identified reads)\n")
                for s in sorted(empty_samples):
                    fh.write(f"{s}\n")
            log_info(f"   (logged to {output_path / 'empty_samples.txt'})")
        except Exception:
            pass
    if not lca_files:
        log_info("\n❌ All samples had 0 identified reads — no OTU table to build.")
        # write empty placeholder tables so downstream steps don't crash on missing files
        empty = pd.DataFrame()
        for name in ["raw", *ranks]:
            export_otu_table(empty, output_path / f"otu_table_{name}.csv", transpose=True)
        return {name: empty for name in ["raw", *ranks]}

    log_info(f"✅ Merging {len(lca_files)} sample(s) with classified reads.")

    results = {}

    # 1. Raw OTU table (LCA taxids)
    log_info("\n1. Generating raw OTU table (LCA taxids)...")
    raw_otu = build_otu_table(lca_files, tax=tax, normalize=False)
    filtered_raw = filter_otu_table(raw_otu, min_count=min_count, min_samples=min_samples)

    if normalize:
        filtered_raw = filtered_raw.div(filtered_raw.sum(axis=1), axis=0)

    export_otu_table(filtered_raw, output_path / "otu_table_raw.csv", transpose=True)
    results["raw"] = filtered_raw

    # 2. Rank-specific OTU tables
    for rank in ranks:
        log_info(f"\n2. Generating {rank}-level OTU table...")
        rank_otu = build_otu_table_at_rank(lca_files, tax, rank=rank, normalize=False)
        filtered_rank = filter_otu_table(rank_otu, min_count=min_count, min_samples=min_samples)

        if normalize:
            filtered_rank = filtered_rank.div(filtered_rank.sum(axis=1), axis=0)

        export_otu_table(filtered_rank, output_path / f"otu_table_{rank}.csv", transpose=True)
        results[rank] = filtered_rank

    log_info("\n" + "=" * 70)
    log_info("OTU TABLE GENERATION COMPLETE")
    log_info("=" * 70)
    log_info(f"✅ Generated {len(results)} OTU tables in {output_path}")

    return results
