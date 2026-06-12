#!/usr/bin/env python3
"""
Fig3 (FULL cohort) PCoA — coloured by the FOUR protocol x disease groups:
  DNA-Control, DNA-Asthma, RNA-Control, RNA-Asthma.

Same Bray-Curtis PCoA ordination as render_fig3_full.py panel B (per tool, pooling
all DNA+RNA samples; coordinates read straight from fig3_pcoa_pooled_full.csv so the
geometry is identical to the DNA-vs-RNA figure), but each point is coloured by its
(protocol, disease-group) combination instead of by protocol alone, with a 95%
ellipse per group. PERMANOVA is recomputed on the 4-level grouping (999 perm).

Group membership: paper/figures/Fig3/source_data/{dna,rna}_asthma_groups.csv.
Writes fig3/Fig3_pcoa_4group.{eps,png}.
Run: conda run -n shotgun_virome python result_260605/fig3/code/render_fig3_pcoa_4group.py
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.2,
})
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.patches import Patch, Ellipse

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/scripts"); sys.path.insert(0, f"{NX}/scripts/benchmark")
from diversity_metrics import bray_curtis, pcoa, permanova

T = f"{NX}/result_260605/fig3/tables"
OUT = f"{NX}/result_260605/fig3"
GROUPS = f"{NX}/paper/figures/Fig3/source_data"
TOOLS = ["NexVirome", "Ganon", "Phanta", "Kraken2", "Metabuli"]

# four (protocol, disease) groups; orange/green family = DNA/RNA, dark/light = Control/Asthma
GKEY = ["DNA-Control", "DNA-Asthma", "RNA-Control", "RNA-Asthma"]
GCOL = {"DNA-Control": "#b35900", "DNA-Asthma": "#f0a500",   # dark / light orange
        "RNA-Control": "#006653", "RNA-Asthma": "#42c79a"}   # dark / light green


def group_map():
    """sample -> 'DNA-Control' etc. (cohort prefix from the manifest seq_type)."""
    m = {}
    for f, coh in (("dna_asthma_groups.csv", "dna"), ("rna_asthma_groups.csv", "rna")):
        g = pd.read_csv(f"{GROUPS}/{f}")
        for _, r in g.iterrows():
            m[(coh, str(r["sample"]))] = f"{coh.upper()}-{r['group']}"
    return m


def _ellipse(ax, x, y, color, n_std=2.0):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)) or cov[0, 0] <= 0 or cov[1, 1] <= 0:
        return
    pear = cov[0, 1] / np.sqrt(cov[0, 0]*cov[1, 1])
    rx, ry = np.sqrt(1+pear), np.sqrt(1-pear)
    e = Ellipse((0, 0), width=rx*2, height=ry*2, facecolor="none",
                edgecolor=color, lw=1.8, zorder=1)
    tr = (mtransforms.Affine2D().rotate_deg(45)
          .scale(np.sqrt(cov[0, 0])*n_std, np.sqrt(cov[1, 1])*n_std)
          .translate(np.mean(x), np.mean(y)))
    e.set_transform(tr + ax.transData); ax.add_patch(e)


def cell(ax, d, ptxt, tool):
    ax.set_box_aspect(1)
    if d.empty:
        ax.text(0.5, 0.5, "(no viral)", ha="center", va="center", transform=ax.transAxes,
                fontsize=20, color="grey")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(tool, fontsize=21, labelpad=8); return
    for gk in GKEY:
        g = d[d.grp == gk]
        if g.empty:
            continue
        _ellipse(ax, g.PC1.values, g.PC2.values, GCOL[gk])
        ax.scatter(g.PC1, g.PC2, s=34, color=GCOL[gk], alpha=0.85, edgecolor="white",
                   linewidth=0.4, zorder=3)
    v1, v2 = d.var1.iloc[0]*100, d.var2.iloc[0]*100
    ax.set_xlabel(f"PC1 ({v1:.0f}%)", fontsize=18); ax.set_ylabel(f"PC2 ({v2:.0f}%)", fontsize=18)
    ax.tick_params(labelsize=15)
    ax.text(0.04, 0.96, ptxt, transform=ax.transAxes, ha="left", va="top", fontsize=15)
    ax.text(0.5, -0.27, tool, transform=ax.transAxes, ha="center", va="top", fontsize=21)


def main():
    pcoa_df = pd.read_csv(f"{T}/fig3_pcoa_pooled_full.csv")
    gm = group_map()
    pcoa_df["grp"] = pcoa_df.apply(lambda r: gm.get((r.cohort, str(r["sample"]))), axis=1)
    miss = pcoa_df[pcoa_df.grp.isna()]
    if len(miss):
        print(f"WARN: {len(miss)} PCoA points have no group label:",
              sorted(set(zip(miss.cohort, miss["sample"]))))
        pcoa_df = pcoa_df.dropna(subset=["grp"])

    # PERMANOVA on the 4-level grouping (rebuild the Bray-Curtis matrix per tool from
    # species_long so the test matches the coords; coords themselves come from the CSV).
    sp = pd.read_csv(f"{T}/fig3_species_long.csv")
    perm = {}
    for tool in TOOLS:
        d = sp[sp.tool == tool].copy()
        if d.empty:
            continue
        d["skey"] = d["cohort"] + ":" + d["sample"]
        mat = d.pivot_table(index="skey", columns="taxid", values="reads",
                            aggfunc="sum", fill_value=0)
        keep = mat.sum(axis=1) > 0
        M = mat[keep].to_numpy(dtype=float); ss = list(mat[keep].index)
        labels = []
        for sk in ss:
            coh, sn = sk.split(":", 1)
            labels.append(gm.get((coh, sn)))
        if M.shape[0] >= 4 and M.shape[1] >= 2 and len(set(labels)) >= 2:
            try:
                F, p = permanova(bray_curtis(M), labels, n_perm=999, seed=0)
                perm[tool] = (F, p)
            except Exception:
                pass

    perm_rows = [dict(tool=t, F=round(F, 4), p=round(p, 4)) for t, (F, p) in perm.items()]
    pd.DataFrame(perm_rows).to_csv(f"{T}/fig3_permanova_4group_full.csv", index=False)
    print("=== PERMANOVA (4 groups, 999 perm) — FULL cohort ===")
    print(pd.DataFrame(perm_rows).to_string(index=False))

    fig, axes = plt.subplots(1, len(TOOLS), figsize=(18.0, 4.7))
    for ax, tool in zip(axes, TOOLS):
        d = pcoa_df[pcoa_df.tool == tool]
        F, p = perm.get(tool, (np.nan, np.nan))
        ptxt = f"PERM p={p:.3f}" if p == p else ""
        cell(ax, d, ptxt, tool)

    handles = [Patch(facecolor=GCOL[gk]) for gk in GKEY]
    fig.legend(handles, GKEY, loc="lower center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, -0.01), fontsize=17, columnspacing=2.0,
               handletextpad=0.6)
    fig.subplots_adjust(left=0.045, right=0.99, top=0.93, bottom=0.34, wspace=0.32)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/Fig3_pcoa_4group.{ext}", dpi=300)
    print(f"-> {OUT}/Fig3_pcoa_4group.png , .eps")
    print("n per group:")
    print(pcoa_df.groupby(["tool", "grp"])["sample"].nunique().to_string())


if __name__ == "__main__":
    main()
