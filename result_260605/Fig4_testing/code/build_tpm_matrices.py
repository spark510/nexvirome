#!/usr/bin/env python3
"""
Fig4_testing — TPM relative-abundance matrices for MaAsLin2 differential
abundance (Control vs Asthma, within each library protocol).

MaAsLin2 takes a RELATIVE-ABUNDANCE feature table (sample x feature) + a metadata
table (sample x variable). It applies its own normalisation/transform (default
TSS + LOG). We feed the genome-length-normalised TPM (rel_abund) from
fig4_species_long.csv, renormalised within the analysed subset so each sample
sums to 1.

Four feature levels (each as sample x feature TPM, one file per DNA/RNA):
  1. kind        phage vs non-phage
  2. host_genus  phage rolled to host genus (phage subset, renorm within phage)
  3. genus       non-phage virus genus      (non-phage subset, renorm)
  4. species     non-phage virus species    (non-phage subset, renorm)

Outputs (result_260605/Fig4_testing/data/):
  tpm_kind_{DNA,RNA}.csv               (sample x feature)
  tpm_hostgenus_{DNA,RNA}.csv
  tpm_nonphage_genus_{DNA,RNA}.csv
  tpm_nonphage_species_{DNA,RNA}.csv
  meta_{DNA,RNA}.csv                   (sample x group)   [shared with counts]

Run: conda run -n shotgun_virome python \
       result_260605/Fig4_testing/code/build_tpm_matrices.py
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/result_260605/fig4/code")
import _fig4_persample_common as C

OUT = f"{NX}/result_260605/Fig4_testing/data"
SPECIES_LONG = f"{NX}/result_260605/fig4/tables/fig4_species_long.csv"


def _tpm_matrix(df, feature_col, samples):
    """sample x feature TPM matrix, renormalised so each sample sums to 1 over the
    features present in `df` (the subset). Zero-filled to the full sample set."""
    # raw genome-length-normalised weight = rel_abund (already per-sample TPM over
    # ALL detections); within a subset we renormalise to sum=1 per sample.
    d = df.copy()
    st = d.groupby("sample")["rel_abund"].transform("sum")
    d["w"] = np.where(st > 0, d["rel_abund"] / st, 0.0)
    m = (d.groupby([feature_col, "sample"])["w"].sum().unstack(fill_value=0.0))
    m = m.reindex(columns=samples, fill_value=0.0)
    m = m.loc[m.sum(axis=1) > 0]              # drop all-zero features
    return m.T                                 # -> sample x feature (MaAsLin2)


def main():
    os.makedirs(OUT, exist_ok=True)
    raw = pd.read_csv(SPECIES_LONG)
    raw["taxid"] = raw["taxid"].astype(int)
    phage, t2hg = C.build_phage_and_host()
    raw["is_phage"] = raw["taxid"].isin(phage)
    raw["kind"] = np.where(raw["is_phage"], "phage", "non-phage")
    raw["host_genus"] = raw["taxid"].map(t2hg)
    raw["nonphage_genus"] = raw["taxid"].map(lambda t: C.genus_name_of(t))

    man = C.manifest()
    for proto, cohort in [("DNA", "dna"), ("RNA", "rna")]:
        sm = man[man["cohort"] == cohort][["sample", "group"]].copy()
        sm = sm.sort_values(["group", "sample"]).reset_index(drop=True)
        sm.set_index("sample").to_csv(f"{OUT}/meta_{proto}.csv")
        samples = sm["sample"].tolist()
        sub = raw[raw["cohort"] == cohort]

        _tpm_matrix(sub, "kind", samples).to_csv(f"{OUT}/tpm_kind_{proto}.csv")
        ph = sub[sub["is_phage"] & sub["host_genus"].notna()]
        _tpm_matrix(ph, "host_genus", samples).to_csv(f"{OUT}/tpm_hostgenus_{proto}.csv")
        npg = sub[~sub["is_phage"]]
        _tpm_matrix(npg, "nonphage_genus", samples).to_csv(f"{OUT}/tpm_nonphage_genus_{proto}.csv")
        _tpm_matrix(npg, "name", samples).to_csv(f"{OUT}/tpm_nonphage_species_{proto}.csv")

        nc = (sm["group"] == "Control").sum(); na = (sm["group"] == "Asthma").sum()
        print(f"{proto}: {len(samples)} samples (Control {nc} / Asthma {na})")
        for lvl in ["kind", "hostgenus", "nonphage_genus", "nonphage_species"]:
            n = pd.read_csv(f"{OUT}/tpm_{lvl}_{proto}.csv", index_col=0).shape[1]
            print(f"    {lvl:18s}: {n} features")
    print(f"\n-> {OUT}/  (tpm_*_{{DNA,RNA}}.csv + meta_{{DNA,RNA}}.csv)")


if __name__ == "__main__":
    main()
