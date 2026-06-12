#!/usr/bin/env python3
"""
Fig4_testing — COMBINED (phage + non-phage in ONE feature table) TPM matrices for
MaAsLin2 Control-vs-Asthma. phage -> '<HostGenus> phage'; non-phage -> its own
genus (or species). Renormalised per sample over the WHOLE integrated set (sum=1).

Outputs (data/): tpm_combined_genus_{DNA,RNA}.csv, tpm_combined_species_{DNA,RNA}.csv
"""
import os, sys
import numpy as np, pandas as pd
NX="/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/result_260605/fig4/code")
import _fig4_persample_common as C
OUT=f"{NX}/result_260605/Fig4_testing/data"
SL=f"{NX}/result_260605/fig4/tables/fig4_species_long.csv"

def tpm_combined(df, level, samples):
    d=df.copy()
    st=d.groupby("sample")["rel_abund"].transform("sum")
    d["w"]=np.where(st>0, d["rel_abund"]/st, 0.0)
    m=d.groupby(["label_"+level,"sample"])["w"].sum().unstack(fill_value=0.0)
    m=m.reindex(columns=samples, fill_value=0.0)
    m=m.loc[m.sum(axis=1)>0]
    return m.T  # sample x feature

raw=pd.read_csv(SL); raw["taxid"]=raw["taxid"].astype(int)
phage,t2hg=C.build_phage_and_host()
raw["is_phage"]=raw["taxid"].isin(phage)
raw["host_genus"]=raw["taxid"].map(t2hg)
def lab(level,row):
    if row["is_phage"]:
        hg=row["host_genus"]; return f"{hg} phage" if pd.notna(hg) else "Other-host phage"
    return C.genus_name_of(row["taxid"]) if level=="genus" else str(row["name"])
raw["label_genus"]=raw.apply(lambda r: lab("genus",r),axis=1)
raw["label_species"]=raw.apply(lambda r: lab("species",r),axis=1)

man=C.manifest()
for proto,cohort in [("DNA","dna"),("RNA","rna")]:
    sm=man[man["cohort"]==cohort][["sample","group"]].sort_values(["group","sample"])
    samples=sm["sample"].tolist()
    sub=raw[raw["cohort"]==cohort]
    for level in ["genus","species"]:
        tpm_combined(sub,level,samples).to_csv(f"{OUT}/tpm_combined_{level}_{proto}.csv")
        n=pd.read_csv(f"{OUT}/tpm_combined_{level}_{proto}.csv",index_col=0).shape[1]
        print(f"{proto} combined {level}: {n} features")
print("done")
