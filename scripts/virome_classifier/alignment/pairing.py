"""
Vectorized paired-end concordance annotation for combined R1+R2 hits.

The classifier merges R1 (strand '+') and R2 (strand '-') hits into one frame and
groups by normalized query, but it does NOT use the fact that a true species is
supported by BOTH mates while close-relative cross-mapping usually is not.

This module adds that signal cheaply (vectorized groupby, ~15s for 600k hits vs
~8min for a per-read Python loop). Validated on KIT: species-level "both-mate
support fraction" separates TP (median 0.84) from FP (median 0.16) strongly.

Adds columns to the hits frame (when both strands present):
  - mate_support : for this (query,target-species), how many mates (1 or 2) hit it
  - pair_both    : bool, mate_support == 2 (R1 and R2 agree on this species)
These can feed coverage judgment and the ML filter as features.
"""
from __future__ import annotations

import pandas as pd


def annotate_pairing(hits: pd.DataFrame, tax, species_col: str = "sp") -> pd.DataFrame:
    """Return hits with `mate_support` and `pair_both` columns (vectorized)."""
    if "strand" not in hits.columns or hits.empty:
        hits = hits.copy()
        hits["mate_support"] = 1
        hits["pair_both"] = False
        return hits

    h = hits.copy()
    if species_col not in h.columns:
        # normalize taxid -> species (cache-friendly via map over unique taxids)
        uniq = h["taxid"].astype(int).unique()
        sp_map = {t: (tax.get_taxid_at_rank(int(t), "species") or int(t)) for t in uniq}
        h[species_col] = h["taxid"].astype(int).map(sp_map)

    # how many distinct mates (strands) support each (query, species)
    sup = (h.groupby(["query", species_col])["strand"]
             .nunique().rename("mate_support").reset_index())
    h = h.merge(sup, on=["query", species_col], how="left")
    h["mate_support"] = h["mate_support"].fillna(1).astype(int)
    h["pair_both"] = h["mate_support"] >= 2
    return h


def species_pair_support(hits: pd.DataFrame, tax, species_col: str = "sp") -> pd.DataFrame:
    """Per-species pair-support summary: reads, both_mate_reads, both_frac.

    `both_frac` (fraction of a species' read-pairs supported by both mates) is the
    strong TP/FP discriminator (TP~0.84 vs FP~0.16 on KIT)."""
    h = annotate_pairing(hits, tax, species_col)
    # one row per (query, species): supported by 1 or 2 mates
    qs = (h.groupby(["query", species_col])["mate_support"].max().reset_index())
    agg = (qs.groupby(species_col)
             .agg(reads=("query", "count"),
                  both_mate_reads=("mate_support", lambda x: int((x >= 2).sum())))
             .reset_index())
    agg["both_frac"] = agg["both_mate_reads"] / agg["reads"].clip(lower=1)
    return agg.rename(columns={species_col: "species"})
