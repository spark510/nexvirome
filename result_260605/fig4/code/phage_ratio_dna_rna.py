#!/usr/bin/env python3
"""
DNA library vs RNA library — is the phage : non-phage ratio different?

Per sample we compute the phage fraction (phage / (phage + non-phage)) two ways:
  - READ-based  : raw NexVirome method-B read counts
  - TPM-based   : genome-length-normalised (the basis used by the Fig4 panels)
then compare DNA cohort vs RNA cohort with a Wilcoxon rank-sum test.

KEY FINDING (after the GOLDEN n>=3 per-taxon floor): the DNA library has a
significantly higher NON-phage fraction than RNA, by BOTH metrics:
  read-based : DNA non-phage ~2x RNA;  Wilcoxon p ~ 0.0015 (r=-0.42)
  TPM-based  : DNA still higher non-phage;  p ~ 0.023 (r=-0.30)
(Before the n>=3 floor the TPM test was n.s. (p~0.53); removing 1-2 read
cross-map noise sharpened the signal so it now holds under TPM too.) DNA's
non-phage are herpes/anello (mid-size dsDNA), abundant enough to survive
length normalisation once noise is gone.

Inputs
------
- result_260603/fig3/tables/fig3_species_long.csv  (tool==NexVirome)
- paper/figures/Fig5_extra/tables/species_genome_length.csv
- resources/db/custom/tax_seq_v20260526_MSL41.db  (phage set)

Outputs (result_260603/fig4/tables/)
-------------------------------------------
- phage_ratio_per_sample_DNA_RNA.csv   cohort, sample, metric, phage_frac, nonphage_frac
- phage_ratio_test_DNA_RNA.csv         metric, dna_median, rna_median, dna_mean,
                                       rna_mean, U, p, rank_biserial, n_dna, n_rna

Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260603/fig4/code/phage_ratio_dna_rna.py
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/result_260605")
from golden_rule import keep_samples
DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
MIX = f"{NX}/result_260605/fig4/tables/fig4_species_long.csv"
GLEN = f"{NX}/paper/figures/Fig5_extra/tables/species_genome_length.csv"
OUT = f"{NX}/result_260605/fig4/tables"


def load():
    sys.path.insert(0, f"{NX}/scripts")
    from virome_classifier.classification.phage_host_rollup import build_phage_host_map
    mix = pd.read_csv(MIX)
    nv = mix[mix["tool"].isin(["NexVirome", "NexVirome lca"])].copy() \
        if "tool" in mix.columns else mix.copy()
    if "read_count" not in nv.columns and "reads" in nv.columns:
        nv = nv.rename(columns={"reads": "read_count"})
    # result_260605: drop EXCLUDE_SAMPLES (vir17).
    nv = nv[nv["sample"].isin(keep_samples(nv["sample"]))].copy()
    nv["taxid"] = nv["taxid"].astype(int)
    nv["read_count"] = nv["read_count"].astype(int)
    _t, phage = build_phage_host_map(DB)
    nv["is_phage"] = nv["taxid"].isin(phage)
    glen = dict(pd.read_csv(GLEN).itertuples(index=False))
    L = nv["taxid"].map(glen).fillna(0).astype(float)
    nv["rpk"] = np.where(L > 0, nv["read_count"] / (L / 1000.0), 0.0)
    return nv


def per_sample(nv, value):
    g = nv.groupby(["cohort", "sample", "is_phage"])[value].sum().unstack(fill_value=0.0)
    g = g.rename(columns={False: "nonphage", True: "phage"})
    for c in ("phage", "nonphage"):
        if c not in g.columns:
            g[c] = 0.0
    g["tot"] = g["phage"] + g["nonphage"]
    g = g[g["tot"] > 0].copy()
    g["phage_frac"] = g["phage"] / g["tot"]
    g["nonphage_frac"] = g["nonphage"] / g["tot"]
    return g.reset_index()[["cohort", "sample", "phage_frac", "nonphage_frac"]]


def main():
    os.makedirs(OUT, exist_ok=True)
    nv = load()

    long_rows, test_rows = [], []
    for value, metric in (("read_count", "read"), ("rpk", "TPM")):
        d = per_sample(nv, value)
        d["metric"] = metric
        long_rows.append(d)
        dna = d[d["cohort"] == "dna"]["phage_frac"].values
        rna = d[d["cohort"] == "rna"]["phage_frac"].values
        U, p = mannwhitneyu(dna, rna, alternative="two-sided")
        r = 2 * U / (len(dna) * len(rna)) - 1
        test_rows.append(dict(
            metric=metric, n_dna=len(dna), n_rna=len(rna),
            dna_median=round(float(np.median(dna)), 4),
            rna_median=round(float(np.median(rna)), 4),
            dna_mean=round(float(dna.mean()), 4), rna_mean=round(float(rna.mean()), 4),
            U=float(U), p=round(float(p), 4), rank_biserial=round(float(r), 3)))

    long_df = pd.concat(long_rows, ignore_index=True)
    test_df = pd.DataFrame(test_rows)
    long_df.to_csv(f"{OUT}/phage_ratio_per_sample_DNA_RNA.csv", index=False)
    test_df.to_csv(f"{OUT}/phage_ratio_test_DNA_RNA.csv", index=False)

    print("=== phage fraction: DNA vs RNA (Wilcoxon) ===")
    print(test_df.to_string(index=False))
    print(f"\n-> {OUT}/phage_ratio_per_sample_DNA_RNA.csv")
    print(f"-> {OUT}/phage_ratio_test_DNA_RNA.csv")


if __name__ == "__main__":
    main()
