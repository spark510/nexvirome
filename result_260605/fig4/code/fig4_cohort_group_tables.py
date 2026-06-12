#!/usr/bin/env python3
"""
Fig4-asthma source tables — phage host-genus (B') and non-phage genus (C')
composition split by cohort x clinical group:
    DNA·Control, DNA·Asthma, RNA·Control, RNA·Asthma   (4 bars each)

This is the asthma/control-stratified expansion of the main Fig4 panels B and C.
Self-contained: reads the same single detection source as Fig2/Fig3/Fig4
(NexVirome method B, mask_v3_full) and computes TPM relative abundance the same
way (genome-length-normalised, per-sample then averaged).

Inputs
------
- Detections (NexVirome method B): result_260603/fig3/tables/fig3_species_long.csv
    (tool==NexVirome; cohort, sample, taxid, name, reads, rel_abund_TPM)
- Sample group manifest: paper/figures/Fig3/source_data/{dna,rna}_asthma_groups.csv
- Genome length (TPM): paper/figures/Fig5_extra/tables/species_genome_length.csv
- Taxonomy DB (genus rollup, phage set, host-genus map):
    resources/db/custom/tax_seq_v20260526_MSL41.db

Outputs (result_260603/fig4/tables/)
-------------------------------------------
- hostgenus_cohort_group.csv       cohort, group, host_genus, frac   (B')
- nonphage_genus_cohort_group.csv  cohort, group, label, frac        (C')

Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260603/fig4/code/build_tables.py
"""
from __future__ import annotations
import os, sqlite3, sys
import numpy as np
import pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/result_260605")
from golden_rule import keep_samples
DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
# result_260605: source the complete, asthma-inclusive, vir17-EXCLUDED detection
# table that THIS folder owns (built by fig4_build_species_long.py), not Fig3's
# control-only fig3_species_long.csv.
MIX = f"{NX}/result_260605/fig4/tables/fig4_species_long.csv"
GLEN = f"{NX}/paper/figures/Fig5_extra/tables/species_genome_length.csv"
DNA_GROUPS = f"{NX}/paper/figures/Fig3/source_data/dna_asthma_groups.csv"
RNA_GROUPS = f"{NX}/paper/figures/Fig3/source_data/rna_asthma_groups.csv"
OUT = f"{NX}/result_260605/fig4/tables"

# 7 named categories + 'Other' so the per-panel legend fits the panel height even
# at the larger (1.5x) legend font; the rest are pooled into Other.
TOP_N = 7
BY = ["cohort", "group"]


# ----------------------------- inputs + TPM -----------------------------
def load_detections():
    mix = pd.read_csv(MIX)
    nv = mix[mix["tool"].isin(["NexVirome", "NexVirome lca"])].copy() \
        if "tool" in mix.columns else mix.copy()
    if "read_count" not in nv.columns and "reads" in nv.columns:
        nv = nv.rename(columns={"reads": "read_count"})
    nv["taxid"] = nv["taxid"].astype(int)
    nv["read_count"] = nv["read_count"].astype(int)

    # per-detection TPM weight (reads / (genome_len/1000), per-sample-normalised)
    glen = dict(pd.read_csv(GLEN).itertuples(index=False))
    L = nv["taxid"].map(glen).fillna(0).astype(float)
    rpk = np.where(L > 0, nv["read_count"] / (L / 1000.0), 0.0)
    nv["tpm"] = rpk
    st = nv.groupby(["cohort", "sample"])["tpm"].transform("sum")
    nv["tpm"] = np.where(st > 0, nv["tpm"] / st, 0.0)

    dna = pd.read_csv(DNA_GROUPS)[["sample", "group"]].assign(cohort="dna")
    rna = pd.read_csv(RNA_GROUPS)[["sample", "group"]].assign(cohort="rna")
    man = pd.concat([dna, rna], ignore_index=True)
    # result_260605: drop EXCLUDE_SAMPLES (vir17) from both the manifest and the
    # detections so it can never enter a cohort x group aggregate.
    man = man[man["sample"].isin(keep_samples(man["sample"]))].reset_index(drop=True)
    nv = nv[nv["sample"].isin(keep_samples(nv["sample"]))].copy()
    # fig4_species_long already carries a `group` column from the manifest; drop it
    # so the authoritative manifest merge below does not create group_x/group_y.
    nv = nv.drop(columns=[c for c in ("group",) if c in nv.columns])
    nv = nv.merge(man[["cohort", "sample", "group"]], on=["cohort", "sample"], how="left")
    return nv


def _mean_tpm(df, label_col, by=BY):
    """Mean TPM composition over label_col, grouped by `by`. tpm sums to 1 per
    sample; sum per label per sample, average over n_samples of each `by` group,
    renormalise so each `by` group sums to 1."""
    g = df.groupby(by + ["sample", label_col])["tpm"].sum().reset_index()
    nsamp = df.groupby(by)["sample"].nunique()
    cm = g.groupby(by + [label_col])["tpm"].sum().reset_index()
    cm = cm.merge(nsamp.rename("n").reset_index(), on=by)
    cm["frac"] = cm["tpm"] / cm["n"]
    tot = cm.groupby(by)["frac"].transform("sum")
    cm["frac"] = cm["frac"] / tot
    return cm[by + [label_col, "frac"]]


def _topn(cm, label_col, other_label, top_n=TOP_N):
    rank = cm.groupby(label_col)["frac"].max().sort_values(ascending=False)
    keep = list(rank.head(top_n).index)
    cm[label_col] = cm[label_col].where(cm[label_col].isin(keep), other_label)
    return cm.groupby(BY + [label_col])["frac"].sum().reset_index()


# ----------------------------- taxonomy lookups -----------------------------
def phage_set(db):
    sys.path.insert(0, f"{NX}/scripts")
    from virome_classifier.classification.phage_host_rollup import build_phage_host_map
    _t, phage = build_phage_host_map(db)
    return phage


def hostgenus_map(db):
    env = ("metagenome", "sludge", "seawater", "sediment", "environment", "soil",
           "wastewater", "uncultured")
    def hg(h):
        if not h or not h.strip():
            return None
        low = h.lower()
        if any(k in low for k in env):
            return "(metagenome)"
        # phage can't truly have a human host: 'Homo sapiens' is a RefSeq metadata
        # error (taxid 38018 'Bacteriophage sp.') -> drop, don't make 'Homo phage'
        if low.startswith("homo sapiens") or low.startswith("homo "):
            return None
        return h.split()[0]
    con = sqlite3.connect(db)
    acc2host = {a: h.strip() for a, h in con.execute(
        "SELECT accession, host FROM refseq_metadata WHERE host IS NOT NULL AND trim(host)!=''")}
    try:
        for a, h in con.execute("SELECT base_accession, host FROM phage_host_from_title"):
            acc2host.setdefault(a, h)
    except sqlite3.OperationalError:
        pass
    t2hg = {}
    for acc, tx in con.execute("SELECT accession, taxid FROM refseq_sequences"):
        if tx in t2hg:
            continue
        g = hg(acc2host.get(acc.split(".")[0]))
        if g:
            t2hg[tx] = g
    con.close()
    return t2hg


def genus_map(taxids, db):
    sys.path.insert(0, f"{NX}/scripts")
    from virome_classifier.taxonomy import TaxonomyDB
    tax = TaxonomyDB.from_sqlite(db)
    out = {}
    for t in set(int(x) for x in taxids):
        try:
            g = tax.get_taxid_at_rank(t, "genus")
        except Exception:
            g = None
        out[t] = (tax.get_name(int(g)) if g and g > 0 else None) or tax.get_name(t) or f"taxid_{t}"
    return out


# ----------------------------- table builders -----------------------------
def t_hostgenus(nv, phage, db):
    """B' — phage host-genus composition (TPM) by cohort x group."""
    t2hg = hostgenus_map(db)
    df = nv[nv["taxid"].isin(phage)].copy()
    df["host_genus"] = df["taxid"].map(t2hg)
    df = df[df["host_genus"].notna()]
    st = df.groupby(["cohort", "sample"])["tpm"].transform("sum")
    df["tpm"] = np.where(st > 0, df["tpm"] / st, 0.0)        # renorm within phage
    cm = _mean_tpm(df, "host_genus")
    return _topn(cm, "host_genus", "Other host")


def t_nonphage_genus(nv, phage, db):
    """C' — non-phage (eukaryotic virus) genus composition (TPM) by cohort x group."""
    df = nv[~nv["taxid"].isin(phage)].copy()
    df["label"] = df["taxid"].map(genus_map(df["taxid"].unique(), db))
    st = df.groupby(["cohort", "sample"])["tpm"].transform("sum")
    df["tpm"] = np.where(st > 0, df["tpm"] / st, 0.0)        # renorm within non-phage
    cm = _mean_tpm(df, "label")
    return _topn(cm, "label", "Other")


def main():
    os.makedirs(OUT, exist_ok=True)
    print("Loading detections (NexVirome method B, from fig3_species_long)...")
    nv = load_detections()
    print(f"  rows={len(nv)}, cohort x group sample counts:")
    print(nv.groupby(["cohort", "group"])["sample"].nunique().to_string())

    phage = phage_set(DB)
    print(f"  phage taxids: {len(phage):,}")

    hg = t_hostgenus(nv, phage, DB)
    ng = t_nonphage_genus(nv, phage, DB)
    hg.to_csv(f"{OUT}/hostgenus_cohort_group.csv", index=False)
    ng.to_csv(f"{OUT}/nonphage_genus_cohort_group.csv", index=False)
    print(f"-> {OUT}/hostgenus_cohort_group.csv ({len(hg)} rows)")
    print(f"-> {OUT}/nonphage_genus_cohort_group.csv ({len(ng)} rows)")

    for name, t, lc in (("B' host-genus", hg, "host_genus"),
                        ("C' non-phage genus", ng, "label")):
        t = t.copy(); t["pct"] = (t["frac"] * 100).round(1)
        t["cg"] = t["cohort"].str.upper() + "-" + t["group"]
        pv = t.pivot_table(index=lc, columns="cg", values="pct", fill_value=0.0)
        pv["mx"] = pv.max(axis=1)
        print(f"\n=== {name} (TPM %) ===")
        print(pv.sort_values("mx", ascending=False).drop(columns="mx").to_string())


if __name__ == "__main__":
    main()
