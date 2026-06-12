#!/usr/bin/env python3
"""
Fig4 dedicated species table — ALL 81 samples (asthma + control), NexVirome only.

WHY THIS EXISTS
  Fig4's panels were sourcing detections from result_260603/fig3/tables/
  fig3_species_long.csv, but that table is the HEALTHY (control-only) Fig3 cohort:
  it contains DNA vir21-25 (control) and RNA A1024O0001-0031 (control) and NO
  asthma samples. So every Fig4 asthma bar (DNA·Asthma n=20, RNA·Asthma n=25) was
  silently empty / stale. This builder classifies EVERY HQ-cache sample with the
  canonical Method-B call so Fig4 has its own complete, asthma-inclusive source.

METHOD  (GOLDEN_RULE.md — identical to fig3_realdata_aggregate.nexvirome_counts):
  best-hit per read -> unmasked breadth >= 0.01 (mask_v3_full) -> per-taxon n >= 3.
  rel_abund = TPM-style genome-length-normalised, sum=1 per sample.

OUTPUT  result_260603/fig4/tables/
  fig4_species_long.csv   cohort, group, sample, taxid, name, reads, rel_abund
  fig4_sample_richness.csv  cohort, group, sample, n_species, total_reads (incl. 0)

Run: conda run -n shotgun_virome python result_260603/fig4/code/fig4_build_species_long.py
"""
from __future__ import annotations
import os, sys, glob
import numpy as np, pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/scripts"); sys.path.insert(0, f"{NX}/scripts/benchmark")
sys.path.insert(0, f"{NX}/result_260605")
import three_strategy_breadth_lca as TS
from virome_classifier.alignment.filters.filter import MaskingFilter
from virome_classifier.taxonomy import TaxonomyDB
from golden_rule import (DB, MASK, GENOME_LENGTH_CSV, HQ_CACHE,
                         BREADTH_CUT, MIN_TAXON_READS, DNA_GROUPS, RNA_GROUPS,
                         apply_method_b, keep_samples, EXCLUDE_SAMPLES)

OUT = f"{NX}/result_260605/fig4/tables"
os.makedirs(OUT, exist_ok=True)

# ---- wire up the classifier exactly like fig3_realdata_aggregate ----
TS.DB = DB; TS._TAX_PATH = DB; TS._TLEN = TS._load_tlen()
TAX = TaxonomyDB.from_sqlite(DB)
GLEN_MAP = dict(pd.read_csv(GENOME_LENGTH_CSV).itertuples(index=False))
_MF = MaskingFilter.from_dataframe(
    pd.read_csv(MASK, sep="\t", header=None, usecols=[0, 1, 2],
                names=["target", "start", "end"]))

_name = {}
def nm(t):
    if t not in _name:
        _name[t] = TAX.get_name(int(t)) or f"taxid_{t}"
    return _name[t]


def nexvirome_counts(sample):
    """{taxid: reads} via the canonical Method-B call on this sample's HQ parquet."""
    df = TS._add_tlen(pd.read_parquet(f"{HQ_CACHE}/{sample}.parquet"))
    return apply_method_b(df, _MF, TS, breadth_cut=BREADTH_CUT,
                          min_taxon_reads=MIN_TAXON_READS)


def rpk_abund(counts):
    """TPM-style genome-length-normalised relative abundance, sum=1 per sample."""
    base = {}
    for t, r in counts.items():
        glen = GLEN_MAP.get(int(t))
        base[t] = (r / (glen / 1000.0)) if glen and glen > 0 else float(r)
    s = sum(base.values())
    return {t: v / s for t, v in base.items()} if s else {}


def main():
    # sample -> (cohort, group) from the manifests; cohort also = dna/rna by id
    dna_g = pd.read_csv(DNA_GROUPS)[["sample", "group"]].assign(cohort="dna")
    rna_g = pd.read_csv(RNA_GROUPS)[["sample", "group"]].assign(cohort="rna")
    manifest = pd.concat([dna_g, rna_g], ignore_index=True)
    # result_260605: drop EXCLUDE_SAMPLES (vir17) from the manifest so it is never
    # iterated/classified and cannot appear in any output row or richness count.
    manifest = manifest[manifest["sample"].isin(keep_samples(manifest["sample"]))] \
        .reset_index(drop=True)
    grp = {(r.cohort, r["sample"]): r["group"] for _, r in manifest.iterrows()}

    # every sample present in BOTH the manifest and the HQ cache
    have = {os.path.basename(p)[:-8] for p in glob.glob(f"{HQ_CACHE}/*.parquet")}
    rows, summ = [], []
    missing_cache = []
    for _, mr in manifest.iterrows():
        cohort, s, g = mr["cohort"], mr["sample"], mr["group"]
        if s not in have:
            missing_cache.append((cohort, g, s))
            summ.append(dict(cohort=cohort, group=g, sample=s,
                             n_species=0, total_reads=0, in_cache=False))
            continue
        c = {int(t): int(r) for t, r in nexvirome_counts(s).items() if r > 0}
        ra = rpk_abund(c)
        for t, r in c.items():
            rows.append(dict(cohort=cohort, group=g, sample=s, taxid=t,
                             name=nm(t), reads=r, rel_abund=round(ra.get(t, 0.0), 6)))
        summ.append(dict(cohort=cohort, group=g, sample=s,
                         n_species=len(c), total_reads=int(sum(c.values())),
                         in_cache=True))

    sp = pd.DataFrame(rows)
    sp.to_csv(f"{OUT}/fig4_species_long.csv", index=False)
    sm = pd.DataFrame(summ)
    sm.to_csv(f"{OUT}/fig4_sample_richness.csv", index=False)

    print(f"-> {OUT}/fig4_species_long.csv  ({len(sp)} rows)")
    print(f"-> {OUT}/fig4_sample_richness.csv  ({len(sm)} samples)")
    if missing_cache:
        print(f"\n[!] {len(missing_cache)} manifest samples have NO HQ-cache parquet "
              f"(counted as 0-detection):")
        for cohort, g, s in missing_cache:
            print(f"      {cohort:3s} {g:8s} {s}")
    print("\n=== samples WITH detection, per cohort x group ===")
    det = sm[sm["n_species"] > 0]
    print(sm.groupby(["cohort", "group"]).agg(
        n_samples=("sample", "nunique"),
        n_in_cache=("in_cache", "sum"),
        n_with_detection=("n_species", lambda x: (x > 0).sum()),
        mean_species=("n_species", "mean"),
        total_reads=("total_reads", "sum")).round(1).to_string())


if __name__ == "__main__":
    main()
