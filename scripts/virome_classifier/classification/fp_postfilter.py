"""
Mode-agnostic false-positive post-filter.

The unique-fraction and genus-competition gates were originally embedded inside
the coverage classifier, so they only worked in coverage/ml_filter modes. This
module re-implements them on the COMMON output schema (the per-read `lca_df`
emitted by every mode: lca / coverage / em), so the same FP gates apply uniformly
to any mode — `--mode lca --min-unique-fraction 0.2` now works too.

lca_df columns used: query, lca_taxid, n_unique_taxids (==1 => uniquely assigned
read). Per species we aggregate total reads and unique reads (n_unique_taxids==1),
then drop species failing the gates and remove their reads.
"""
from __future__ import annotations

from typing import Optional

import pandas as pd

from ..core import log_info


def apply_fp_postfilter(
    lca_df: pd.DataFrame,
    tax=None,
    min_unique_fraction: float = 0.0,
    min_rel_abundance: float = 0.0,
    min_read_count: int = 0,
) -> pd.DataFrame:
    """Filter lca_df by relative-abundance (default) and/or unique-fraction
    and/or an absolute per-taxon read floor (min_read_count).

    The DEFAULT false-positive control is the relative-abundance gate
    (min_rel_abundance): a taxon assigned < this fraction of all classified reads
    is dropped. On the KIT mock this removes same-genus cross-mapping shadows
    (e.g. Cytomegalovirus paninebeta2 at 0.001% vs human CMV; the true species sit
    at >=0.06%, a 40-600x margin) AND part of the contaminant tail, depth-robustly
    (a ratio, not an absolute read count). This is the standard metagenomics
    low-abundance cut and replaced the former genus-competition gate (which could
    demote a genuine same-genus co-detection and gave MORE FPs / fewer true taxa).
    unique-fraction remains as an optional extra gate, off by default.

    Returns lca_df with reads of rejected taxa removed. No-op when all gates off."""
    if lca_df is None or lca_df.empty:
        return lca_df
    if min_rel_abundance <= 0 and min_unique_fraction <= 0 and min_read_count <= 0:
        return lca_df

    df = lca_df.copy()
    df["lca_taxid"] = df["lca_taxid"].astype(int)
    is_unique = df.get("n_unique_taxids", 1).fillna(1).astype(int) == 1

    # per-species totals
    per = df.groupby("lca_taxid").size().rename("total")
    uniq = df[is_unique].groupby("lca_taxid").size().rename("unique")
    stats = pd.concat([per, uniq], axis=1).fillna(0)
    stats["ufrac"] = stats["unique"] / stats["total"].clip(lower=1)

    keep = set(stats.index)

    # 0) relative-abundance gate (DEFAULT): taxon reads / all classified reads
    if min_rel_abundance > 0:
        total_reads = int(stats["total"].sum())
        before = len(keep)
        keep = {t for t in keep
                if stats.loc[t, "total"] / max(total_reads, 1) >= min_rel_abundance}
        log_info(f"  [FP-postfilter] rel-abundance>={min_rel_abundance:.4%}: "
                 f"{before} -> {len(keep)} taxa")

    # 1) unique-fraction gate (optional)
    if min_unique_fraction > 0:
        before = len(keep)
        keep = {t for t in keep if stats.loc[t, "ufrac"] >= min_unique_fraction}
        log_info(f"  [FP-postfilter] unique-fraction>={min_unique_fraction}: "
                 f"{before} -> {len(keep)} species")

    # 2) absolute read-count floor (optional): drop a taxon supported by fewer
    # than min_read_count reads. Depth-DEPENDENT (unlike the rel-abundance gate),
    # so off by default; the LOCKED policy value is 10 (pass --min-read-count 10).
    if min_read_count > 0:
        before = len(keep)
        keep = {t for t in keep if stats.loc[t, "total"] >= min_read_count}
        log_info(f"  [FP-postfilter] read-count>={min_read_count}: "
                 f"{before} -> {len(keep)} taxa")

    out = df[df["lca_taxid"].isin(keep)].copy()
    log_info(f"  [FP-postfilter] reads {len(df):,} -> {len(out):,}, "
             f"species kept {len(keep)}")
    return out
