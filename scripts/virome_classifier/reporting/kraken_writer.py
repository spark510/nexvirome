"""
Kraken-format report writer for taxonomy classification results.

This module provides functions to export LCA classification results
in Kraken-compatible formats:
- Per-read classification output (.kraken)
- Hierarchical taxonomic report (.kreport)
- Abundance table (.tsv)
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional
import numpy as np
import pandas as pd

from ..core import log_info, log_verbose
from ..taxonomy import TaxonomyDB
from .rank_utils import get_rank_code


def attach_taxonomy_columns(results_df: pd.DataFrame, tax: TaxonomyDB) -> pd.DataFrame:
    """
    Add taxonomy metadata columns to results DataFrame.

    Adds: taxon_name, taxon_rank, rank_code based on lca_taxid column.
    Uses caching for performance optimization.

    Args:
        results_df: DataFrame with 'lca_taxid' column
        tax: TaxonomyDB instance

    Returns:
        DataFrame with added taxonomy columns
    """
    if results_df is None or results_df.empty or "lca_taxid" not in results_df.columns:
        return results_df

    log_verbose("🔗 Attaching taxonomy metadata (name/rank/code) to results...")

    # Create caches
    name_cache: Dict[int, str] = {}
    rank_cache: Dict[int, str] = {}
    code_cache: Dict[int, str] = {}

    def _name(tid: int) -> str:
        if tid <= 0:
            return ""
        if tid in name_cache:
            return name_cache[tid]
        name_cache[tid] = tax.get_name(tid) or ""
        return name_cache[tid]

    def _rank(tid: int) -> str:
        if tid <= 0:
            return "no rank"
        if tid in rank_cache:
            return rank_cache[tid]
        rank_cache[tid] = tax.get_rank(tid) or "no rank"
        return rank_cache[tid]

    def _code(tid: int) -> str:
        if tid <= 0:
            return "-"
        if tid in code_cache:
            return code_cache[tid]
        code_cache[tid] = get_rank_code(_rank(tid))
        return code_cache[tid]

    out = results_df.copy()
    tids = out["lca_taxid"].fillna(0).astype(int)

    out["taxon_name"] = tids.map(_name)
    out["taxon_rank"] = tids.map(_rank)
    out["rank_code"] = tids.map(_code)

    return out


def build_abundance_from_results(
    results_df: pd.DataFrame,
    tax: TaxonomyDB,
    *,
    virus_root: bool = True,
) -> pd.DataFrame:
    """
    Build clade-level abundance table with hierarchical structure.

    Accumulates read counts up the taxonomic tree for each classified read.

    Args:
        results_df: DataFrame with 'lca_taxid' and optional 'read_count' columns
        tax: TaxonomyDB instance
        virus_root: If True, treat Viruses (taxid=10239) as effective root

    Returns:
        DataFrame with columns:
        - lca_taxid: Taxonomy ID
        - taxon_name: Taxon name
        - taxon_rank: Taxonomic rank
        - read_count: Total reads in clade (including descendants)
        - abundance: Fraction of total reads
        - depth: Hierarchical depth from root
        - parent_taxid: Parent taxon ID
    """
    if results_df is None or results_df.empty or "lca_taxid" not in results_df.columns:
        return pd.DataFrame(
            columns=["lca_taxid", "taxon_name", "taxon_rank", "read_count", "abundance", "depth", "parent_taxid"]
        )

    # Get read counts (default to 1 if not present)
    rc = results_df["read_count"] if "read_count" in results_df.columns else pd.Series(1, index=results_df.index)
    rc = rc.fillna(1).astype(int)
    tids = results_df["lca_taxid"].fillna(0).astype(int)
    total_reads = int(rc.sum())

    log_verbose(f"📈 Building clade abundance from {len(results_df):,} rows (total reads={total_reads:,})...")

    clade_counts: Dict[int, int] = {}
    all_taxids_in_lineages = set()

    # Accumulate counts up the lineage for each assigned taxid
    for tid, cnt in zip(tids, rc):
        if tid <= 0 or not tax.exists(int(tid)):
            continue

        # Get full lineage
        lineage = tax.get_lineage(int(tid))
        if not lineage:
            lineage = [int(tid)]

        # Remove general root (taxid=1) but keep virus root (taxid=10239)
        if virus_root:
            lineage = [t for t in lineage if t != 1]
        else:
            lineage = [t for t in lineage if t not in [1, 10239]]

        # Accumulate counts for all ancestors
        for anc in lineage:
            clade_counts[anc] = clade_counts.get(anc, 0) + int(cnt)
            all_taxids_in_lineages.add(anc)

    if not clade_counts:
        return pd.DataFrame(
            columns=["lca_taxid", "taxon_name", "taxon_rank", "read_count", "abundance", "depth", "parent_taxid"]
        )

    # Build rows with metadata, depth, and parent info
    rows: List[Dict] = []
    for tid in all_taxids_in_lineages:
        if tid not in clade_counts:
            continue

        count = clade_counts[tid]
        name = tax.get_name(int(tid)) or f"taxid_{tid}"
        rank = tax.get_rank(int(tid)) or "no rank"

        # Get lineage for depth calculation
        lineage = tax.get_lineage(int(tid))

        # Adjust lineage based on virus_root setting
        if virus_root:
            lineage = [t for t in lineage if t != 1]
            # If Viruses is in lineage, calculate depth from Viruses
            if 10239 in lineage:
                virus_idx = lineage.index(10239)
                depth = len(lineage) - virus_idx - 1
            else:
                depth = len(lineage) - 1
        else:
            lineage = [t for t in lineage if t not in [1, 10239]]
            depth = len(lineage) - 1

        # Get parent (from adjusted lineage)
        parent_tid = None
        if len(lineage) > 1:
            current_idx = lineage.index(tid) if tid in lineage else -1
            if current_idx > 0:
                parent_tid = lineage[current_idx - 1]

        rows.append(
            {
                "lca_taxid": int(tid),
                "taxon_name": name,
                "taxon_rank": rank,
                "read_count": int(count),
                "abundance": (float(count) / float(total_reads)) if total_reads > 0 else 0.0,
                "depth": depth,
                "parent_taxid": int(parent_tid) if parent_tid else None,
            }
        )

    df = pd.DataFrame(rows)
    df["parent_taxid"] = df["parent_taxid"].astype("Int64")

    # Sort by depth first, then by read count
    df = df.sort_values(["depth", "read_count"], ascending=[True, False]).reset_index(drop=True)
    return df


def write_kraken_output(
    results_df: pd.DataFrame,
    output_file: str,
    tax: TaxonomyDB,
) -> None:
    """
    Write per-read classification in Kraken format.

    Format:
        C/U <query_id> <lca_taxid> <query_length> <lineage>

    Where:
        - C/U: Classified or Unclassified
        - query_id: Read identifier
        - lca_taxid: Assigned taxonomy ID (0 if unclassified)
        - query_length: Read length
        - lineage: Pipe-separated taxids from root to assignment

    Args:
        results_df: DataFrame with columns: query, lca_taxid, qlen
        output_file: Output file path
        tax: TaxonomyDB instance
    """
    if results_df is None or results_df.empty:
        log_info("⚠️  No results to write (kraken output)")
        return

    try:
        # Vectorised. The old `for _, row in results_df.iterrows()` + per-read
        # tax.get_lineage() dominated deep-coverage samples (1.7M reads => ~130s
        # on Qiagen). Reads share very few distinct lca_taxids, so resolve the
        # lineage string ONCE per distinct taxid, map it onto the column, then
        # build all lines with a single vectorised join — byte-identical output.
        df = results_df
        taxids = pd.to_numeric(df.get("lca_taxid", 0), errors="coerce").fillna(0).astype(int)
        qlens = pd.to_numeric(df.get("qlen", 100), errors="coerce").fillna(100).astype(int)
        queries = df.get("query", pd.Series(["unknown"] * len(df))).astype(str)

        lineage_cache: dict = {}
        def _lineage_str(t):
            s = lineage_cache.get(t)
            if s is None:
                if t > 0 and tax.exists(t):
                    lin = tax.get_lineage(t)
                    s = "|".join(str(x) for x in lin) if lin else str(t)
                else:
                    s = "0"
                lineage_cache[t] = s
            return s

        lin_strs = taxids.map(_lineage_str).astype(str)
        cls = pd.Series(np.where(taxids.values > 0, "C", "U"), index=taxids.index)
        # build each tab-delimited line with pandas string concatenation
        # (numpy '+' on fixed-width <U arrays truncates / errors), then write once.
        lines = (cls.str.cat([queries.astype(str), taxids.astype(str),
                              qlens.astype(str), lin_strs], sep="\t") + "\n")
        with open(output_file, "w") as f:
            f.writelines(lines.tolist())

        log_info(f"📝 Kraken output written to: {output_file}")
    except Exception as e:
        log_info(f"❌ Error writing Kraken output: {e}")
        raise


def generate_kraken_report(abundance_df: pd.DataFrame, output_file: str) -> None:
    """
    Write Kraken-style hierarchical taxonomic report.

    Format:
        <percentage> <clade_reads> <taxon_reads> <rank_code> <taxid> <indented_name>

    Where:
        - percentage: Percentage of total reads in this clade
        - clade_reads: Reads assigned to this taxon + all descendants
        - taxon_reads: Reads assigned directly to this taxon only
        - rank_code: Single-letter rank code (D/K/P/C/O/F/G/S/-)
        - taxid: Taxonomy ID
        - indented_name: Taxon name with indentation showing hierarchy

    Args:
        abundance_df: DataFrame from build_abundance_from_results()
        output_file: Output file path
    """
    if abundance_df is None or abundance_df.empty:
        log_info("⚠️  No abundance to write (kraken report)")
        return

    try:
        # Find the effective root (highest level node with reads)
        abundance_df_sorted = abundance_df.sort_values("depth", ascending=True)

        # Get the root node (lowest depth, highest read count)
        root_candidates = abundance_df_sorted[abundance_df_sorted["depth"] == abundance_df_sorted["depth"].min()]
        if len(root_candidates) > 1:
            root_row = root_candidates.loc[root_candidates["read_count"].idxmax()]
        else:
            root_row = root_candidates.iloc[0]

        # Use root node's read count as 100% base
        total_reads = int(root_row["read_count"])
        root_taxid = int(root_row["lca_taxid"])

        log_verbose(f"📊 Using root taxid={root_taxid} ({root_row['taxon_name']}) with {total_reads:,} reads as 100% base")

        # Build parent-child relationships
        present = set(abundance_df["lca_taxid"].astype(int))
        parent_map: Dict[int, Optional[int]] = {}
        children: Dict[int, List[int]] = {tid: [] for tid in present}

        for _, row in abundance_df.iterrows():
            tid = int(row["lca_taxid"])
            parent_tid = row.get("parent_taxid")
            if pd.notna(parent_tid):
                parent_tid = int(parent_tid)
                if parent_tid in present:
                    parent_map[tid] = parent_tid
                    children[parent_tid].append(tid)
                else:
                    parent_map[tid] = None
            else:
                parent_map[tid] = None

        # Calculate reads_taxon (direct assignment vs clade total)
        read_clade: Dict[int, int] = {int(t): int(c) for t, c in zip(abundance_df["lca_taxid"], abundance_df["read_count"])}

        read_taxon: Dict[int, int] = {}
        # Process in reverse depth order (leaves first)
        order = abundance_df.sort_values("depth", ascending=False)["lca_taxid"].astype(int).tolist()
        for tid in order:
            child_sum = sum(read_clade.get(ch, 0) for ch in children.get(tid, []))
            read_taxon[tid] = max(read_clade.get(tid, 0) - child_sum, 0)

        # Create hierarchical ordering for output
        def get_hierarchical_order():
            """Get taxids in proper hierarchical display order"""
            ordered_tids = []
            visited = set()

            # Find roots (nodes with no parent or parent not in our set)
            roots = [tid for tid, parent in parent_map.items() if parent is None or parent not in present]

            def traverse_depth_first(tid, visited_set):
                if tid in visited_set:
                    return
                visited_set.add(tid)
                ordered_tids.append(tid)

                # Sort children by read count (descending)
                child_list = children.get(tid, [])
                child_list.sort(key=lambda x: read_clade.get(x, 0), reverse=True)

                for child in child_list:
                    traverse_depth_first(child, visited_set)

            # Sort roots by read count (descending)
            roots.sort(key=lambda x: read_clade.get(x, 0), reverse=True)

            for root in roots:
                traverse_depth_first(root, visited)

            return ordered_tids

        # Get hierarchically ordered taxids
        ordered_taxids = get_hierarchical_order()

        # Create lookup for row data
        row_data = {}
        for _, row in abundance_df.iterrows():
            tid = int(row["lca_taxid"])
            row_data[tid] = {"name": str(row["taxon_name"]), "rank": str(row["taxon_rank"]), "depth": int(row["depth"])}

        # Write report in hierarchical order
        with open(output_file, "w") as f:
            for tid in ordered_taxids:
                if tid not in row_data:
                    continue

                data = row_data[tid]
                name = data["name"]
                rank = data["rank"]
                depth = data["depth"]

                rank_code = get_rank_code(rank)
                clade = int(read_clade.get(tid, 0))
                taxonly = int(read_taxon.get(tid, 0))

                # Calculate percentage based on root node reads
                percentage = (float(clade) / float(total_reads) * 100.0) if total_reads > 0 else 0.0

                # Indentation based on depth
                indent = "  " * depth

                f.write(f"{percentage:.2f}\t{clade}\t{taxonly}\t{rank_code}\t{tid}\t{indent}{name}\n")

        log_info(f"📊 Kraken report written to: {output_file}")

    except Exception as e:
        log_info(f"❌ Error generating Kraken report: {e}")
        raise


def write_abundance_table(abundance_df: pd.DataFrame, output_file: str) -> None:
    """
    Write flat abundance table for downstream analysis.

    Simple TSV format with columns:
        lca_taxid, taxon_name, taxon_rank, read_count, abundance

    Args:
        abundance_df: DataFrame from build_abundance_from_results()
        output_file: Output file path
    """
    if abundance_df is None or abundance_df.empty:
        log_info("⚠️  No abundance to write (abundance table)")
        return

    try:
        cols = ["lca_taxid", "taxon_name", "taxon_rank", "read_count", "abundance"]
        existing = [c for c in cols if c in abundance_df.columns]
        abundance_df[existing].to_csv(output_file, sep="\t", index=False)
        log_info(f"📈 Abundance table written to: {output_file}")
    except Exception as e:
        log_info(f"❌ Error writing abundance table: {e}")
        raise


def write_all_outputs(
    results_df: pd.DataFrame,
    tax: TaxonomyDB,
    output_dir: str,
    sample_name: str,
    *,
    virus_root: bool = True,
) -> Dict[str, str]:
    """
    Write all Kraken-format outputs in one call.

    Generates three output files:
    1. {sample_name}.kraken - Per-read classification
    2. {sample_name}.kreport - Hierarchical report
    3. {sample_name}.abundance.tsv - Flat abundance table

    Args:
        results_df: DataFrame with columns: query, lca_taxid, qlen
        tax: TaxonomyDB instance
        output_dir: Output directory path
        sample_name: Sample name for output files
        virus_root: If True, use Viruses (taxid=10239) as root

    Returns:
        Dictionary mapping file types to paths:
        {'kraken_output': path, 'kraken_report': path, 'abundance_table': path}
    """
    if results_df is None or results_df.empty:
        log_info("⚠️  No results to write")
        return {}

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    files = {
        "kraken_output": str(Path(output_dir) / f"{sample_name}.kraken"),
        "kraken_report": str(Path(output_dir) / f"{sample_name}.kreport"),
        "abundance_table": str(Path(output_dir) / f"{sample_name}.abundance.tsv"),
    }

    log_info(f"📁 Writing Kraken outputs to: {output_dir}")

    # 1) Attach taxonomy metadata
    results_df = attach_taxonomy_columns(results_df, tax)

    # 2) Build clade-level abundance
    abundance_df = build_abundance_from_results(results_df, tax, virus_root=virus_root)

    # 3) Write all outputs
    write_kraken_output(results_df, files["kraken_output"], tax)

    if abundance_df is not None and not abundance_df.empty:
        generate_kraken_report(abundance_df, files["kraken_report"])
        write_abundance_table(abundance_df, files["abundance_table"])
    else:
        log_info("⚠️  Abundance is empty; skipping report/table")

    # 4) Taxon roll-ups (NCBI + ICTV × species/genus). Independent of the kreport
    #    so the manuscript and downstream tools can ingest a consistent unit;
    #    requires a TaxonomyDB with a sqlite db_path attribute.
    try:
        from .taxon_rollup_writer import write_taxon_rollups
        db_path = getattr(tax, "_db_path", None)
        if db_path is not None:
            rollups = write_taxon_rollups(results_df, tax, str(db_path),
                                          output_dir, sample_name)
            files.update(rollups)
    except Exception as e:
        log_info(f"⚠️  taxon roll-up writer skipped: {e}")

    log_info("✅ All Kraken output files written successfully")

    return files
