"""
LCA classification.

LCAClassifier owns the per-read lowest-common-ancestor assignment that was
previously implemented as free functions in cli/classify.py. Logic is moved
verbatim (byte-identical results) — only the packaging changed.

- classify()            : plain LCA (natural, or rank-fixed to species/genus)
- classify_conditional(): hierarchical LCA that retreats low-confidence species
                          calls up to genus
"""
from __future__ import annotations

import pandas as pd

from ..core import log_info
from ..taxonomy import TaxonomyDB


class LCAClassifier:
    """Lowest-common-ancestor read classifier over alignment hits."""

    def __init__(self, tax: TaxonomyDB, verbose: bool = False):
        self.tax = tax
        self.verbose = verbose

    def classify(
        self,
        hits: pd.DataFrame,
        fix_rank: str = "none",
    ) -> pd.DataFrame:
        """Perform LCA classification.

        fix_rank: 'none' (default) keeps the natural LCA — compute_lca stops at the
        real common ancestor, which varies per read (strain..root). 'species' or
        'genus' lifts each read's LCA UP to that rank (get_taxid_at_rank); a read
        whose natural LCA is already above the target rank (e.g. a cross-family read
        at Viruses root) is left unchanged, since that rank does not exist on its
        lineage. Rank-fixing makes the detection unit uniform: 'species' ≈ natural
        here (most reads already resolve near species), 'genus' merges same-genus
        relatives into one taxon (loses species resolution, dissolves the
        genus-competition gate — see fp_postfilter)."""
        tax = self.tax
        log_info("\n🧬 Performing LCA classification...")

        # Vectorised per-query aggregation. The previous `for query, group in
        # hits.groupby("query")` Python loop dominated deep-coverage samples (98万
        # groups => pandas get_iterator/_chop ~107s on Qiagen). Here we aggregate
        # each query's taxid list / qlen / hit count with a single groupby.agg, then
        # compute_lca ONCE per distinct taxid-tuple (most queries share a tuple), so
        # the result is byte-identical to the loop:
        #   taxids = first-seen-order unique taxids of the query  (== Series.unique())
        #   qlen   = first hit's qlen ; n_hits = len(group)
        has_qlen = "qlen" in hits.columns

        def _first_unique(s):
            # pandas Series.unique() preserves first-seen order, matching the old code
            return tuple(pd.unique(s).tolist())

        agg = {"taxid": _first_unique, "target": "size"}
        if has_qlen:
            agg["qlen"] = "first"
        g = hits.groupby("query", sort=False).agg(agg)
        g = g.rename(columns={"target": "n_hits"})

        # LCA once per distinct taxid-tuple
        lca_cache: dict = {}
        def _lca(tup):
            v = lca_cache.get(tup)
            if v is None:
                lt = tax.compute_lca(list(tup))
                if lt and lt > 0 and fix_rank in ("species", "genus"):
                    lt = tax.get_taxid_at_rank(int(lt), fix_rank)
                v = lca_cache[tup] = lt
            return v

        g["lca_taxid"] = g["taxid"].map(_lca)
        g = g[g["lca_taxid"].notna() & (g["lca_taxid"] > 0)]

        name_cache: dict = {}
        rank_cache: dict = {}
        def _nm(t):
            if t not in name_cache: name_cache[t] = tax.get_name(t)
            return name_cache[t]
        def _rk(t):
            if t not in rank_cache: rank_cache[t] = tax.get_rank(t)
            return rank_cache[t]

        lca_df = pd.DataFrame({
            "query": g.index,
            "lca_taxid": g["lca_taxid"].astype(int).values,
            "lca_name": [_nm(t) for t in g["lca_taxid"].astype(int)],
            "lca_rank": [_rk(t) for t in g["lca_taxid"].astype(int)],
            "qlen": (g["qlen"].astype(int).values if has_qlen else 100),
            "read_count": 1,
            "n_hits": g["n_hits"].astype(int).values,
            "n_unique_taxids": [len(t) for t in g["taxid"]],
            "all_taxids": [",".join(map(str, t)) for t in g["taxid"]],
        }).reset_index(drop=True)
        log_info(f"✅ Classified {len(lca_df):,} queries")

        if len(lca_df) > 0:
            log_info("\nLCA rank distribution:")
            for rank, count in lca_df["lca_rank"].value_counts().items():
                log_info(f"  {rank}: {count:,}")

        return lca_df

    def classify_conditional(
        self,
        hits: pd.DataFrame,
        min_species_confidence: float = 0.5,
    ) -> pd.DataFrame:
        """Hierarchical/conditional LCA: a read assigned at SPECIES rank is only kept
        at species if the species has enough confident support; otherwise that read is
        'held up' to the genus rank. This is the probabilistic generalisation of plain
        LCA (which never retreats) — it avoids fabricating low-confidence species calls,
        reducing false positives, while keeping the higher-rank signal.

        Confidence proxy (no training): per species, the fraction of its species-level
        reads that map UNIQUELY to it (n_unique_taxids == 1). Below the threshold the
        species' reads are reassigned to the species' genus (via the taxonomy tree).
        Keeps `lca` mode untouched; this is a separate path.
        """
        tax = self.tax
        log_info("\n🧬 Performing CONDITIONAL (hierarchical) LCA classification...")
        base = self.classify(hits)
        if base.empty:
            return base

        df = base.copy()

        # Confidence of a SPECIES-rank call = its species-rank read share among all
        # reads whose hits touch that species' GENUS. Rationale: a species-rank read
        # always maps to one species (that's why LCA reached species), so "unique
        # fraction" is trivially 1 and never retreats. Instead we ask: of all reads in
        # this genus's neighbourhood, how many resolve cleanly to THIS species vs. stay
        # ambiguous (genus/higher)? If the species captures only a small share of its
        # genus-level evidence, the species call is weakly supported → retreat to genus.
        df["_genus"] = df["lca_taxid"].map(
            lambda t: tax.get_taxid_at_rank(int(t), "genus") or 0)
        sp = df["lca_rank"] == "species"

        # reads per species (species-rank) and reads per genus (any rank under genus)
        sp_reads = df[sp].groupby("lca_taxid").size()
        genus_reads = df.groupby("_genus").size()
        # species confidence = species_reads / genus_reads(its genus)
        sp_genus = df[sp].groupby("lca_taxid")["_genus"].first()
        conf = {tid: (sp_reads[tid] / genus_reads.get(sp_genus[tid], sp_reads[tid]))
                for tid in sp_reads.index}
        low_conf = {tid for tid, c in conf.items() if c < min_species_confidence}

        retreated = 0
        for i in df.index[sp]:
            tid = df.at[i, "lca_taxid"]
            if tid in low_conf:
                genus = tax.get_taxid_at_rank(int(tid), "genus")
                if genus and genus > 0:
                    df.at[i, "lca_taxid"] = genus
                    df.at[i, "lca_rank"] = "genus"
                    df.at[i, "lca_name"] = tax.get_name(genus)
                    retreated += 1
        log_info(f"  conditional retreat: {retreated:,} species-reads -> genus "
                 f"({len(low_conf)} low-confidence species, conf<{min_species_confidence})")
        return df
