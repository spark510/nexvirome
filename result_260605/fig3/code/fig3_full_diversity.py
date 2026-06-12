#!/usr/bin/env python3
"""
Fig3-FULL diversity — ALL samples (Control + Asthma), DNA-vs-RNA comparison.

Identical analysis to healthy_diversity.py but WITHOUT the Control-only filter,
so it uses the entire cohort (DNA n=24, RNA n=56; vir17 already excluded upstream
by golden_rule.keep_samples in fig3_realdata_aggregate.py). The comparison axis is
still the library protocol (DNA vs RNA), not the disease group.

alpha per (cohort, tool, sample); cross-tool Kruskal-Wallis per cohort/index;
per-tool Bray-Curtis PCoA pooling DNA+RNA; PERMANOVA(DNA vs RNA, 999 perm).

Reads result_260605/fig3/tables/fig3_species_long.csv
Writes fig3 alpha/pcoa/kw CSVs with a _full suffix (the full-cohort sibling of the
Healthy fig3_*.csv).
Run: conda run -n shotgun_virome python result_260605/fig3/code/fig3_full_diversity.py
"""
from __future__ import annotations
import os, sys
NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/scripts"); sys.path.insert(0, f"{NX}/scripts/benchmark")
from diversity_metrics import bray_curtis, alpha_indices, pcoa, permanova
import numpy as np, pandas as pd
from scipy.stats import kruskal, mannwhitneyu

T = f"{NX}/result_260605/fig3/tables"
TOOLS = ["NexVirome", "Kraken2", "Phanta", "Ganon", "Metabuli"]
SUF = "_full"   # full-cohort outputs, sibling of the Healthy fig3_*.csv


def main():
    sp = pd.read_csv(f"{T}/fig3_species_long.csv")
    # NO Control filter: use the entire cohort (Control + Asthma).
    alpha_rows, pcoa_rows, kw_rows = [], [], []
    for cohort in ("rna", "dna"):
        for tool in TOOLS:
            d = sp[(sp.cohort == cohort) & (sp.tool == tool)]
            if d.empty:
                continue
            mat = d.pivot_table(index="sample", columns="taxid", values="reads",
                                aggfunc="sum", fill_value=0)
            for s in mat.index:
                a = alpha_indices(mat.loc[s].to_numpy(dtype=float))
                alpha_rows.append(dict(cohort=cohort, tool=tool, sample=s,
                                       richness=int((mat.loc[s] > 0).sum()),
                                       shannon=round(a["shannon"], 4),
                                       invsimpson=round(a["invsimpson"], 4),
                                       pielou=round(a["pielou"], 4)))
            # per-tool PCoA (all samples)
            keep = mat.sum(axis=1) > 0
            M = mat[keep].to_numpy(dtype=float); ss = list(mat[keep].index)
            if M.shape[0] >= 3 and M.shape[1] >= 2:
                coords, var = pcoa(bray_curtis(M), 2)
                for s, xy in zip(ss, coords):
                    pcoa_rows.append(dict(cohort=cohort, tool=tool, sample=s,
                                          PC1=round(float(xy[0]), 4), PC2=round(float(xy[1]), 4),
                                          var1=round(float(var[0]), 4), var2=round(float(var[1]), 4)))
        # cross-tool Kruskal-Wallis per index (do tools differ in alpha?)
        adf = pd.DataFrame([r for r in alpha_rows if r["cohort"] == cohort])
        for idx in ("richness", "shannon", "invsimpson", "pielou"):
            groups = [adf[adf.tool == t][idx].dropna().values for t in TOOLS]
            groups = [g for g in groups if len(g) >= 2]
            if len(groups) >= 2:
                try: p = kruskal(*groups).pvalue
                except Exception: p = np.nan
            else: p = np.nan
            kw_rows.append(dict(cohort=cohort, index=idx, n_tools=len(groups),
                                p=round(float(p), 5) if p == p else np.nan))

    # DNA vs RNA Mann-Whitney per (tool, index) — non-paired
    a = pd.DataFrame(alpha_rows)
    mw_rows = []
    for tool in TOOLS:
        for idx in ("richness", "shannon", "invsimpson", "pielou"):
            dna = a[(a.tool == tool) & (a.cohort == "dna")][idx].dropna().values
            rna = a[(a.tool == tool) & (a.cohort == "rna")][idx].dropna().values
            if len(dna) >= 2 and len(rna) >= 2:
                try: p = mannwhitneyu(dna, rna, alternative="two-sided").pvalue
                except Exception: p = np.nan
            else: p = np.nan
            mw_rows.append(dict(tool=tool, index=idx, n_dna=len(dna), n_rna=len(rna),
                                p=round(float(p), 5) if p == p else np.nan))
    pd.DataFrame(mw_rows).to_csv(f"{T}/fig3_alpha_dna_vs_rna{SUF}.csv", index=False)

    # per-TOOL PCoA pooling DNA+RNA (all) samples + PERMANOVA(DNA vs RNA, 999 perm)
    pooled_rows, perm_rows = [], []
    for tool in TOOLS:
        d = sp[sp.tool == tool]
        if d.empty:
            continue
        d = d.copy(); d["skey"] = d["cohort"] + ":" + d["sample"]
        mat = d.pivot_table(index="skey", columns="taxid", values="reads",
                            aggfunc="sum", fill_value=0)
        keep = mat.sum(axis=1) > 0
        M = mat[keep].to_numpy(dtype=float); ss = list(mat[keep].index)
        labels = [sk.split(":", 1)[0] for sk in ss]   # cohort label per sample
        if M.shape[0] >= 3 and M.shape[1] >= 2:
            D = bray_curtis(M)
            coords, var = pcoa(D, 2)
            for sk, xy in zip(ss, coords):
                coh, sname = sk.split(":", 1)
                pooled_rows.append(dict(tool=tool, cohort=coh, sample=sname,
                                        PC1=round(float(xy[0]), 4), PC2=round(float(xy[1]), 4),
                                        var1=round(float(var[0]), 4), var2=round(float(var[1]), 4)))
            if labels.count("dna") >= 2 and labels.count("rna") >= 2:
                try:
                    F, p = permanova(D, labels, n_perm=999, seed=0)
                    perm_rows.append(dict(tool=tool, F=round(float(F), 4), p=round(float(p), 4),
                                          n_dna=labels.count("dna"), n_rna=labels.count("rna")))
                except Exception:
                    pass
    pd.DataFrame(pooled_rows).to_csv(f"{T}/fig3_pcoa_pooled{SUF}.csv", index=False)
    pd.DataFrame(perm_rows).to_csv(f"{T}/fig3_permanova_dna_vs_rna{SUF}.csv", index=False)
    if perm_rows:
        print("\n=== PERMANOVA DNA vs RNA (999 perm, per tool) — FULL cohort ===")
        print(pd.DataFrame(perm_rows).to_string(index=False))

    pd.DataFrame(alpha_rows).to_csv(f"{T}/fig3_alpha{SUF}.csv", index=False)
    pd.DataFrame(pcoa_rows).to_csv(f"{T}/fig3_pcoa{SUF}.csv", index=False)
    pd.DataFrame(kw_rows).to_csv(f"{T}/fig3_alpha_kruskal{SUF}.csv", index=False)
    pd.set_option("display.width", 140)
    print("=== FULL alpha mean (cohort x tool) ===")
    a = pd.DataFrame(alpha_rows)
    print(a.groupby(["cohort", "tool"])[["richness", "shannon", "invsimpson"]].mean().round(2).to_string())
    print("\n=== n samples per (cohort, tool) ===")
    print(a.groupby(["cohort", "tool"])["sample"].nunique().to_string())
    print(f"\n-> {T}/  (*{SUF}.csv)")


if __name__ == "__main__":
    main()
