"""
Best-hit read classification.

Assigns each read the taxid of its single highest-scoring alignment (Strategy B):
read -> argmax(bits) hit's taxid. This is the user-core read-level assignment that
pins one taxon per read (lowest false-positive rate; cross-map reads do NOT drift to
a higher rank as they would under LCA). It is the alternative to LCAClassifier at
"stage 0" (read-level taxonomy); the downstream rank roll-up, phage-host roll-up,
FP post-filter and kraken/abundance writers all consume the SAME lca_df schema, so
this classifier emits exactly those columns.

Tie-break: highest bits, then highest fident, then smallest taxid (stable) — matching
the benchmark _best_hit so production and the three_strategy evaluation agree.
"""
from __future__ import annotations

import pandas as pd

from ..core import log_info
from ..taxonomy import TaxonomyDB


class BestHitClassifier:
    """Per-read best-hit (max-bitscore) taxonomic assignment."""

    def __init__(self, tax: TaxonomyDB, verbose: bool = False):
        self.tax = tax
        self.verbose = verbose

    def classify(self, hits: pd.DataFrame, fix_rank: str = "none") -> pd.DataFrame:
        """One row per read = its best-hit taxid. fix_rank lifts that taxid UP to
        species/genus (get_taxid_at_rank), 'none' keeps the hit's own rank. Output
        schema is identical to LCAClassifier.classify so all downstream steps work."""
        tax = self.tax
        log_info("\n🧬 Performing BEST-HIT classification...")
        if hits is None or hits.empty:
            return pd.DataFrame()

        has_qlen = "qlen" in hits.columns
        score_cols = [c for c in ("bits", "fident") if c in hits.columns]
        if not score_cols:
            # no score column -> fall back to first hit per read (stable)
            score_cols = []

        # best hit per read: sort by score asc then keep last (= max), stable so
        # ties resolve deterministically; then smallest taxid on remaining ties.
        df = hits
        if score_cols:
            df = df.sort_values(score_cols + ["taxid"],
                                ascending=[True] * len(score_cols) + [False],
                                kind="stable")
        best = df.drop_duplicates("query", keep="last")

        # per-read n_hits / unique-taxid count (for the FP post-filter's ufrac gate)
        grp = hits.groupby("query", sort=False)
        n_hits = grp["target"].size()
        n_uniq = grp["taxid"].nunique()

        taxids = best["taxid"].astype(int).values
        if fix_rank in ("species", "genus"):
            taxids = [self._lift(int(t), fix_rank) for t in taxids]

        keep = [(q, t) for q, t in zip(best["query"].values, taxids)
                if t is not None and t > 0]
        if not keep:
            return pd.DataFrame()
        queries = [q for q, _ in keep]
        taxids = [t for _, t in keep]

        name_cache, rank_cache = {}, {}
        def _nm(t):
            if t not in name_cache:
                name_cache[t] = tax.get_name(t)
            return name_cache[t]
        def _rk(t):
            if t not in rank_cache:
                rank_cache[t] = tax.get_rank(t)
            return rank_cache[t]

        best_q = best.set_index("query")
        out = pd.DataFrame({
            "query": queries,
            "lca_taxid": [int(t) for t in taxids],
            "lca_name": [_nm(int(t)) for t in taxids],
            "lca_rank": [_rk(int(t)) for t in taxids],
            "qlen": [int(best_q.at[q, "qlen"]) if has_qlen else 100 for q in queries],
            "read_count": 1,
            "n_hits": [int(n_hits.get(q, 1)) for q in queries],
            "n_unique_taxids": [int(n_uniq.get(q, 1)) for q in queries],
            "all_taxids": [str(best_q.at[q, "taxid"]) for q in queries],
        }).reset_index(drop=True)
        log_info(f"✅ Best-hit classified {len(out):,} reads")
        return out

    def _lift(self, taxid: int, rank: str):
        lifted = self.tax.get_taxid_at_rank(taxid, rank)
        return lifted if (lifted and lifted > 0) else taxid
