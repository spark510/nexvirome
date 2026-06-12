"""
Phage -> host-genus roll-up (optional, mode-agnostic post-step).

Phage cross-map disperses one read across many phage references of the SAME host
(genus-of-same-phage => same host in 93% of cases). Rolling phage detections up to
their bacterial/archaeal HOST GENUS collapses that dispersion and yields a directly
interpretable read-out ("phages_of_Streptococcus" ~ Streptococcus present). It is
applied AFTER the FP post-filter, only when --phage-host-rollup is set.

Rules (see docs/phage_host_processing.md):
  - phage with a KNOWN host       -> relabel lca_taxid to a synthetic host node and
                                     lca_name to "phages_of_<HostGenus>" (reads of all
                                     same-host phage merge into one taxon).
  - phage with an UNKNOWN host     -> kept PER-SPECIES (NOT merged): they are distinct
                                     phages we merely lack a host for; collapsing them
                                     into one "(host_unknown)" bin would fabricate a
                                     false merge. (left unchanged.)
  - non-phage (herpes/human/...)   -> left at species (never merged).

Host source is VMR-INDEPENDENT: refseq_metadata.host (NCBI-Virus "Host" field) joined
on version-stripped accession, plus the title-parsed phage_host_from_title table.
Environmental/metagenome host strings collapse to unknown (kept per-species).

The synthetic host node uses a NEGATIVE id (-(hash) space) so it never collides with a
real NCBI taxid; kraken/abundance writers key on lca_taxid/lca_name and treat it as a
leaf. is_phage uses the same lineage rule as the benchmark (Caudoviricetes / ssRNA /
ssDNA phage clades + the word 'phage'); Duplodnaviria is NOT treated as phage (it also
contains herpesviruses).
"""
from __future__ import annotations

import sqlite3
from typing import Optional

import pandas as pd

from ..core import log_info

# phage clades (lineage names, lowercased). Mirrors fp/benchmark; excludes
# Duplodnaviria (realm shared with herpesviruses).
_PHAGE_CLADES = {
    "caudoviricetes", "caudovirales", "microviridae", "inoviridae",
    "leviviricetes", "microviricetes", "tectiviridae", "corticoviridae",
    "tubulavirales", "faserviricetes",
}
_ENV_HOST = ("metagenome", "sludge", "seawater", "sediment", "environment",
             "soil", "wastewater", "uncultured")
# lineage clades with an established host the title doesn't spell out
_CLADE_HOST = {"crassvirales": "Bacteroides", "suoliviridae": "Bacteroides"}


def _host_genus(h: Optional[str]) -> Optional[str]:
    """First token of the host name; None for empty / environmental sources, and
    None for a human host — a bacteriophage cannot have Homo sapiens as its true
    host, so a 'Homo sapiens' entry is a RefSeq metadata error (e.g. taxid 38018
    'Bacteriophage sp.') and is dropped rather than rolled up to a 'Homo phage'."""
    if not h or not h.strip():
        return None
    low = h.lower()
    if any(k in low for k in _ENV_HOST):
        return None
    if low.startswith("homo sapiens") or low.startswith("homo "):
        return None
    return h.split()[0]


def build_phage_host_map(db_path: str):
    """Build {taxid: host_genus} (known host only) and a phage-taxid set from the
    taxonomy DB. Host = refseq_metadata.host (primary, version-stripped join) +
    phage_host_from_title (fill) + Crassvirales->Bacteroides lineage rule."""
    con = sqlite3.connect(db_path)
    acc2host = {a: h.strip() for a, h in con.execute(
        "SELECT accession, host FROM refseq_metadata "
        "WHERE host IS NOT NULL AND trim(host)!=''")}
    try:
        for a, h in con.execute("SELECT base_accession, host FROM phage_host_from_title"):
            acc2host.setdefault(a, h)
    except sqlite3.OperationalError:
        pass

    parent, name = {}, {}
    for tx, pt, nm in con.execute(
            "SELECT taxid, parent_taxid, scientific_name FROM ncbi_taxonomy"):
        parent[tx] = pt
        name[tx] = nm

    def lineage_lower(tx):
        out, seen = [], 0
        while tx in parent and seen < 90:
            out.append((name.get(tx, "") or "").lower())
            if tx == parent[tx]:
                break
            tx = parent[tx]
            seen += 1
        return out

    tax2host, phage = {}, set()
    seen_tax = set()
    for acc, tx, title in con.execute(
            "SELECT accession, taxid, title FROM refseq_sequences"):
        ln = set(lineage_lower(tx))
        if ("phage" in (title or "").lower()) or (ln & _PHAGE_CLADES):
            phage.add(tx)
        if tx in seen_tax:
            continue
        seen_tax.add(tx)
        # host: lineage clade rule first (crAssphage), then metadata/title
        hg = None
        for clade, host in _CLADE_HOST.items():
            if clade in ln:
                hg = host
                break
        if hg is None:
            hg = _host_genus(acc2host.get(acc.split(".")[0]))
        if hg:
            tax2host[tx] = hg
    con.close()
    return tax2host, phage


def _synthetic_host_taxid(host_genus: str) -> int:
    """Stable NEGATIVE pseudo-taxid for a host-genus node (never collides with a
    real positive NCBI taxid). Deterministic per host name."""
    return -(abs(hash(("phages_of", host_genus))) % 2_000_000_000) - 1


def apply_phage_host_rollup(lca_df: pd.DataFrame, tax2host: dict, phage: set) -> pd.DataFrame:
    """Relabel phage rows (known host) to a synthetic host-genus taxon; leave
    unknown-host phage per-species and all non-phage untouched. Returns a new df
    with lca_taxid/lca_name/lca_rank rewritten for rolled-up rows.

    No-op-safe: rows whose taxid is not a known-host phage pass through unchanged,
    so downstream kraken/abundance writers see one merged leaf per host genus."""
    if lca_df is None or lca_df.empty:
        return lca_df
    df = lca_df.copy()
    df["lca_taxid"] = df["lca_taxid"].astype(int)

    def remap(tid):
        if tid in phage and tid in tax2host:
            hg = tax2host[tid]
            return _synthetic_host_taxid(hg), f"phages_of_{hg}", "host_genus"
        return None

    n_rolled = 0
    new_tid, new_nm, new_rk = [], [], []
    for tid, nm, rk in zip(df["lca_taxid"], df.get("lca_name", ""), df.get("lca_rank", "")):
        r = remap(int(tid))
        if r:
            new_tid.append(r[0]); new_nm.append(r[1]); new_rk.append(r[2]); n_rolled += 1
        else:
            new_tid.append(int(tid)); new_nm.append(nm); new_rk.append(rk)
    df["lca_taxid"] = new_tid
    df["lca_name"] = new_nm
    df["lca_rank"] = new_rk

    n_hosts = len({t for t in new_tid if t < 0})
    log_info(f"  [phage-host-rollup] {n_rolled:,} phage reads -> {n_hosts} host-genus "
             f"taxa (unknown-host phage & non-phage left per-species)")
    return df
