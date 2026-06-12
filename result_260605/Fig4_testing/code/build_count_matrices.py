#!/usr/bin/env python3
"""
Fig4_testing — build the count matrices for the four differential-abundance
explorations (Control vs Asthma, within each library protocol).

Source: result_260605/fig4/tables/fig4_species_long.csv  (vir17-excluded, 80
samples, NexVirome method-B; `reads` = raw counts, `rel_abund` = per-sample TPM).

ALDEx2 takes RAW INTEGER COUNTS (it does CLR + Dirichlet Monte-Carlo internally),
so every matrix below is integer reads (feature x sample). Four feature levels:

  1. kind        phage vs non-phage         (2 rows)
  2. host_genus  phage rolled to host genus (phage subset)
  3. genus       non-phage virus genus      (non-phage subset)
  4. species     non-phage virus species    (non-phage subset)

For each level we write ONE matrix per protocol (DNA, RNA), columns = that
protocol's samples (Control + Asthma), plus a sample->group table.

Outputs (result_260605/Fig4_testing/data/):
  counts_kind_{DNA,RNA}.csv
  counts_hostgenus_{DNA,RNA}.csv
  counts_nonphage_genus_{DNA,RNA}.csv
  counts_nonphage_species_{DNA,RNA}.csv
  sample_group_{DNA,RNA}.csv     sample, group   (Control / Asthma)

Run: conda run -n shotgun_virome python \
       result_260605/Fig4_testing/code/build_count_matrices.py
"""
from __future__ import annotations
import os, sys
import numpy as np, pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/result_260605/fig4/code")
import _fig4_persample_common as C

OUT = f"{NX}/result_260605/Fig4_testing/data"
SPECIES_LONG = f"{NX}/result_260605/fig4/tables/fig4_species_long.csv"


def _matrix(df, feature_col, samples):
    """Integer feature x sample read-count matrix over `samples` (zero-filled)."""
    m = (df.groupby([feature_col, "sample"])["reads"].sum()
           .unstack(fill_value=0))
    m = m.reindex(columns=samples, fill_value=0)        # full sample set incl zeros
    m = m.loc[m.sum(axis=1) > 0]                          # drop all-zero features
    return m.astype(int)


def main():
    os.makedirs(OUT, exist_ok=True)
    raw = pd.read_csv(SPECIES_LONG)
    raw["taxid"] = raw["taxid"].astype(int)
    phage, t2hg = C.build_phage_and_host()
    raw["is_phage"] = raw["taxid"].isin(phage)

    # labels
    raw["kind"] = np.where(raw["is_phage"], "phage", "non-phage")
    raw["host_genus"] = raw["taxid"].map(t2hg)
    raw["nonphage_genus"] = raw["taxid"].map(lambda t: C.genus_name_of(t))
    # species name is `name`

    # full sample manifest per protocol (incl. zero-detection samples)
    man = C.manifest()   # cohort, sample, group  (vir17 already excluded)

    for proto, cohort in [("DNA", "dna"), ("RNA", "rna")]:
        sm = man[man["cohort"] == cohort][["sample", "group"]].copy()
        sm = sm.sort_values(["group", "sample"]).reset_index(drop=True)
        sm.to_csv(f"{OUT}/sample_group_{proto}.csv", index=False)
        samples = sm["sample"].tolist()
        sub = raw[raw["cohort"] == cohort]

        # 1. phage vs non-phage (all detections)
        _matrix(sub, "kind", samples).to_csv(f"{OUT}/counts_kind_{proto}.csv")
        # 2. phage host genus (phage only)
        ph = sub[sub["is_phage"] & sub["host_genus"].notna()]
        _matrix(ph, "host_genus", samples).to_csv(f"{OUT}/counts_hostgenus_{proto}.csv")
        # 3. non-phage genus
        npg = sub[~sub["is_phage"]]
        _matrix(npg, "nonphage_genus", samples).to_csv(f"{OUT}/counts_nonphage_genus_{proto}.csv")
        # 4. non-phage species
        _matrix(npg, "name", samples).to_csv(f"{OUT}/counts_nonphage_species_{proto}.csv")

        nc = (sm["group"] == "Control").sum(); na = (sm["group"] == "Asthma").sum()
        print(f"{proto}: {len(samples)} samples (Control {nc} / Asthma {na})")
        for lvl in ["kind", "hostgenus", "nonphage_genus", "nonphage_species"]:
            fn = f"{OUT}/counts_{lvl}_{proto}.csv"
            n = pd.read_csv(fn, index_col=0).shape[0]
            print(f"    {lvl:18s}: {n} features")

    print(f"\n-> {OUT}/  (counts_*_{{DNA,RNA}}.csv + sample_group_{{DNA,RNA}}.csv)")


if __name__ == "__main__":
    main()
