"""
Per-sample taxonomy roll-up reports.

The default kreport mixes whatever rank `compute_lca` happens to stop at
(strain / species / genus / ... / root) per read, which is useful raw output but
makes downstream comparison and ICTV-aligned reporting awkward. This module emits
FOUR additional flat TSVs alongside the kreport:

  {sample}.ncbi_species.tsv   reads rolled to NCBI species
  {sample}.ncbi_genus.tsv     reads rolled to NCBI genus
  {sample}.ictv_species.tsv   reads rolled to ICTV species (where mappable)
  {sample}.ictv_genus.tsv     reads rolled to ICTV genus

Each TSV: taxid, name, rank, reads, fraction (of total classified reads). Reads
that cannot be resolved to the requested rank (e.g. assigned at genus while we
ask for species, or an NCBI taxon with no ICTV link) are aggregated into a single
'_unmapped' row so the totals stay honest.

NCBI roll-up uses the in-memory TaxonomyDB (`get_taxid_at_rank`). ICTV roll-up
uses both available bridges and takes the first that resolves: (a) species-normalised
NCBI taxid -> `ictv_taxonomy.ncbi_taxid`, and (b) refseq_sequences.taxid ->
refseq_sequences.ictv_taxid -> ictv_taxonomy. ICTV walk-up then follows ICTV's own
parent_taxid/rank columns, so the result is in ICTV's namespace, not NCBI's.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Optional

import pandas as pd

from ..core import log_info


# ---------- ICTV bridge ----------

class IctvBridge:
    """NCBI taxid -> ICTV taxid lookup + ICTV walk-up to species/genus.

    Built lazily from the sqlite DB the classifier already loads. Cached so a
    single sample run touches the DB once per distinct taxid."""

    def __init__(self, db_path: str):
        self._con = sqlite3.connect(db_path)
        # NCBI -> ICTV index (direct match on ictv_taxonomy.ncbi_taxid)
        self._ncbi2ictv: dict = {}
        for nt, it in self._con.execute(
            "SELECT ncbi_taxid, taxid FROM ictv_taxonomy WHERE ncbi_taxid IS NOT NULL"
        ):
            if nt is not None and nt not in self._ncbi2ictv:
                self._ncbi2ictv[int(nt)] = it
        # NCBI -> ICTV via refseq_sequences (string ictv_taxid values are looked
        # up against ictv_taxonomy.ictv_id separately)
        self._ncbi2refseqictv: dict = {}
        for nt, it in self._con.execute(
            "SELECT taxid, ictv_taxid FROM refseq_sequences "
            "WHERE ictv_taxid IS NOT NULL AND taxid IS NOT NULL"
        ):
            if nt is not None and nt not in self._ncbi2refseqictv:
                self._ncbi2refseqictv[int(nt)] = it
        # ictv_id (string) -> taxid (int) translation, used by refseq bridge
        self._ictvid2taxid: dict = {
            r[0]: r[1] for r in self._con.execute(
                "SELECT ictv_id, taxid FROM ictv_taxonomy WHERE ictv_id IS NOT NULL"
            )
        }
        # ICTV parent / rank / name caches
        self._parent: dict = {}
        self._rank: dict = {}
        self._name: dict = {}
        for tid, p, rk, nm in self._con.execute(
            "SELECT taxid, parent_taxid, rank, name FROM ictv_taxonomy"
        ):
            self._parent[tid] = p
            self._rank[tid] = rk
            self._name[tid] = nm
        self._rank_walk: dict = {}  # (taxid, target_rank) -> ictv taxid

    def ncbi_to_ictv(self, ncbi_taxid: int) -> Optional[int]:
        """Return ICTV taxid for a given NCBI taxid, trying direct then refseq bridge."""
        if ncbi_taxid is None:
            return None
        t = self._ncbi2ictv.get(int(ncbi_taxid))
        if t is not None:
            return t
        ictv_id = self._ncbi2refseqictv.get(int(ncbi_taxid))
        if ictv_id is not None:
            return self._ictvid2taxid.get(ictv_id)
        return None

    def walk_to_rank(self, ictv_taxid: int, target_rank: str) -> Optional[int]:
        """Walk up ICTV lineage until rank matches target_rank; None if unreachable."""
        if ictv_taxid is None:
            return None
        key = (ictv_taxid, target_rank)
        if key in self._rank_walk:
            return self._rank_walk[key]
        cur = ictv_taxid
        for _ in range(30):
            if cur is None or cur not in self._rank:
                self._rank_walk[key] = None
                return None
            if self._rank[cur] == target_rank:
                self._rank_walk[key] = cur
                return cur
            p = self._parent.get(cur)
            if p is None or p == cur:
                self._rank_walk[key] = None
                return None
            cur = p
        self._rank_walk[key] = None
        return None

    def name(self, ictv_taxid: int) -> str:
        return self._name.get(ictv_taxid, "?")


# ---------- roll-up + write ----------

def _write_tsv(rows: dict, path: str, namespace: str, rank: str) -> None:
    total = sum(rows.values()) or 1
    df = pd.DataFrame(
        [(tid, name, reads, reads / total) for (tid, name), reads in rows.items()],
        columns=["taxid", "name", "reads", "fraction"],
    )
    df["namespace"] = namespace
    df["rank"] = rank
    df = df.sort_values("reads", ascending=False)
    df.to_csv(path, sep="\t", index=False)


def write_taxon_rollups(
    results_df: pd.DataFrame,
    tax,
    db_path: str,
    output_dir: str,
    sample_name: str,
) -> dict:
    """Emit NCBI species/genus + ICTV species/genus roll-up TSVs.

    Returns a dict mapping output kind to file path."""
    if results_df is None or results_df.empty:
        return {}
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    bridge = IctvBridge(db_path)
    # per-taxid total reads from the lca_df (already one row per query)
    counts = results_df.groupby("lca_taxid").size().to_dict()
    files = {
        "ncbi_species": str(out / f"{sample_name}.ncbi_species.tsv"),
        "ncbi_genus":   str(out / f"{sample_name}.ncbi_genus.tsv"),
        "ictv_species": str(out / f"{sample_name}.ictv_species.tsv"),
        "ictv_genus":   str(out / f"{sample_name}.ictv_genus.tsv"),
    }

    # NCBI roll-up: walk each lca_taxid up to species / genus
    ncbi_sp: dict = {}
    ncbi_gn: dict = {}
    for tid, n in counts.items():
        try:
            sp = tax.get_taxid_at_rank(int(tid), "species")
        except Exception:
            sp = None
        try:
            gn = tax.get_taxid_at_rank(int(tid), "genus")
        except Exception:
            gn = None
        key_sp = (int(sp), tax.get_name(sp)) if (sp and sp != int(tid) or
                  (sp and tax.get_rank(sp) == "species")) else (0, "_unmapped")
        key_gn = (int(gn), tax.get_name(gn)) if (gn and tax.get_rank(gn) == "genus") else (0, "_unmapped")
        ncbi_sp[key_sp] = ncbi_sp.get(key_sp, 0) + n
        ncbi_gn[key_gn] = ncbi_gn.get(key_gn, 0) + n

    _write_tsv(ncbi_sp, files["ncbi_species"], "NCBI", "species")
    _write_tsv(ncbi_gn, files["ncbi_genus"],   "NCBI", "genus")

    # ICTV roll-up: NCBI taxid -> bridge -> walk ICTV lineage
    ictv_sp: dict = {}
    ictv_gn: dict = {}
    for tid, n in counts.items():
        it = bridge.ncbi_to_ictv(int(tid))
        sp = bridge.walk_to_rank(it, "species") if it else None
        gn = bridge.walk_to_rank(it, "genus") if it else None
        key_sp = (int(sp), bridge.name(sp)) if sp else (0, "_unmapped")
        key_gn = (int(gn), bridge.name(gn)) if gn else (0, "_unmapped")
        ictv_sp[key_sp] = ictv_sp.get(key_sp, 0) + n
        ictv_gn[key_gn] = ictv_gn.get(key_gn, 0) + n

    _write_tsv(ictv_sp, files["ictv_species"], "ICTV", "species")
    _write_tsv(ictv_gn, files["ictv_genus"],   "ICTV", "genus")

    n_unmap_sp = ictv_sp.get((0, "_unmapped"), 0)
    n_unmap_gn = ictv_gn.get((0, "_unmapped"), 0)
    tot = sum(counts.values()) or 1
    log_info(
        f"📁 Taxon roll-ups written (NCBI + ICTV × species/genus). "
        f"ICTV unmapped: species {n_unmap_sp:,}/{tot:,} ({100*n_unmap_sp/tot:.1f}%), "
        f"genus {n_unmap_gn:,}/{tot:,} ({100*n_unmap_gn/tot:.1f}%)."
    )
    return files
