"""
Paired-end concordance classification (3 presets).

MMseqs gives multiple hits per read. For a read PAIR (R1 + R2, which share the
same normalized query id), we combine the two mates' candidate taxa to suppress
close-relative cross-mapping false positives: a true species is supported by
BOTH mates; a spillover relative usually is not.

Presets (chosen vs full grid to keep the design space tractable):
    pe_strict   : per-mate candidates within margin tau (=5) of that mate's best
                  bitscore; assign the INTERSECTION of R1 & R2 candidate taxa
                  (LCA if several). Empty intersection -> unclassified.
    pe_balanced : candidates within margin tau (=10); SUM bitscore per taxon over
                  both mates; assign argmax. If top1-top2 margin < tau -> LCA of
                  the tied taxa.
    pe_coord    : species-level intersection like pe_strict; if empty, fall back
                  to the LCA of the UNION (typically the genus) instead of dropping.

All operate on the combined hits DataFrame produced by classify.parse_alignments
(columns incl. query, target, taxid, bits, strand['+'=R1,'-'=R2]). Output is the
same lca_df schema the Kraken reporting expects.
"""
from __future__ import annotations

from typing import List

import pandas as pd

from ..core import log_info


PRESET_DEFAULTS = {
    "pe_strict":   {"tau": 5.0,  "combine": "intersection", "fallback": "drop"},
    "pe_balanced": {"tau": 10.0, "combine": "weighted_sum", "fallback": "lca"},
    "pe_coord":    {"tau": 10.0, "combine": "intersection", "fallback": "genus"},
}


def _candidates(mate_hits: pd.DataFrame, tau: float) -> set:
    """Taxa whose best bitscore is within `tau` of this mate's overall best."""
    if mate_hits.empty:
        return set()
    best = mate_hits["bits"].max()
    keep = mate_hits[mate_hits["bits"] >= best - tau]
    return set(keep["taxid"].astype(int).tolist())


def _assign_one(group: pd.DataFrame, tax, tau: float, combine: str, fallback: str):
    """Return assigned taxid for one read(pair), or None if unclassified."""
    r1 = group[group["strand"] == "+"]
    r2 = group[group["strand"] == "-"]
    c1, c2 = _candidates(r1, tau), _candidates(r2, tau)

    # Single-mate read: fall back to that mate's candidates (LCA)
    if not c1 or not c2:
        cand = c1 or c2
        return tax.compute_lca(list(cand)) if cand else None

    if combine == "weighted_sum":
        # sum bitscore per taxon over both mates; argmax with margin gating
        sums = group.groupby(group["taxid"].astype(int))["bits"].sum().sort_values(ascending=False)
        if len(sums) == 1:
            return int(sums.index[0])
        top1, top2 = sums.iloc[0], sums.iloc[1]
        if (top1 - top2) < tau:
            tied = [int(t) for t, v in sums.items() if (top1 - v) < tau]
            return tax.compute_lca(tied)
        return int(sums.index[0])

    # intersection-based (pe_strict / pe_coord)
    inter = c1 & c2
    if inter:
        return tax.compute_lca(list(inter)) if len(inter) > 1 else int(next(iter(inter)))
    # empty intersection
    if fallback == "drop":
        return None
    if fallback == "genus":
        lca = tax.compute_lca(list(c1 | c2))
        return lca  # union LCA is typically genus-or-higher
    if fallback == "lca":
        return tax.compute_lca(list(c1 | c2))
    return None


def classify_paired(combined_hits: pd.DataFrame, tax, preset: str = "pe_balanced",
                    tau: float = None, verbose: bool = False) -> pd.DataFrame:
    """Run a paired-end concordance preset; return lca_df (Kraken-ready)."""
    cfg = PRESET_DEFAULTS[preset]
    tau = cfg["tau"] if tau is None else tau
    combine, fallback = cfg["combine"], cfg["fallback"]
    log_info(f"\n🧬 Paired-end concordance [{preset}] (tau={tau}, combine={combine}, fallback={fallback})...")

    if "strand" not in combined_hits.columns:
        combined_hits = combined_hits.assign(strand="+")  # single-end safety

    rows: List[dict] = []
    for query, group in combined_hits.groupby("query"):
        taxid = _assign_one(group, tax, tau, combine, fallback)
        if taxid and taxid > 0:
            rows.append({
                "query": query,
                "lca_taxid": int(taxid),
                "lca_name": tax.get_name(taxid) or f"Unknown ({taxid})",
                "lca_rank": tax.get_rank(taxid) or "no rank",
                "qlen": int(group.iloc[0]["qlen"]) if "qlen" in group.columns else 100,
                "read_count": 1,
                "n_hits": len(group),
                "n_unique_taxids": group["taxid"].nunique(),
                "all_taxids": ",".join(map(str, group["taxid"].astype(int).unique())),
            })

    lca_df = pd.DataFrame(rows)
    log_info(f"✅ [{preset}] classified {len(lca_df):,} read pairs")
    return lca_df
