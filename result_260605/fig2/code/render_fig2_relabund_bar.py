#!/usr/bin/env python3
"""
Fig2 relative-abundance stacked bar: 6 KIT ground-truth species in color + ALL
false positives lumped into one GREY 'FP' block, summing to 100% per tool.
TPM-normalized (genome-length corrected), SPECIES level, mean over 4 KIT samples.

Two panels in ONE figure:
  (a) same-DB : external tools given OUR viral DB (+ NexVirome)
  (b) native  : external tools on their own integrated DB (+ NexVirome, whose
                viral-only DB is its native, so its bar is the same in both)
The two panels share the legend -> one shared legend below the figure.

Abundance sources: Kraken2/Metabuli = Bracken re-estimated; Ganon = own reads;
Phanta = built-in bracken; NexVirome = best-hit reads.

Input : result_260605/fig2/fig2_kit_composition_with_FP.csv  (category incl. 'FP')
Output: result_260605/fig2/fig2_relabund_bar_withFP.{png,pdf}
Run: conda run -n shotgun_virome python result_260605/fig2/code/render_fig2_relabund_bar.py
"""
from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Unified style — identical to Fig2_tpfpf1_per_sample (the figure this sits below).
matplotlib.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.size": 22,
    "axes.titlesize": 24,
    "axes.labelsize": 22,
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "legend.fontsize": 20,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.2,
})

NX = "/home/share/programs/nexvirome"
CSV = f"{NX}/result_260605/fig2/fig2_kit_composition_with_FP.csv"
OUT = f"{NX}/result_260605/fig2"

SPECIES = ["Adenovirus 40", "HHV-5 (CMV)", "Human RSV", "Influenza B", "Reovirus 3", "Zika virus"]
SP_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#CCB974"]
FP_COLOR = "#9e9e9e"
EXPECTED = 100.0 / 6
TOOLS = ["NexVirome", "Ganon", "Kraken2", "Metabuli", "Phanta"]
SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
# R = Roche (MagNA Pure kit); internal sample id MagNA_1 unchanged in the CSVs.
SAMPLE_LABEL = {"MagNA_1": "R1", "MagNA_2": "R2", "Qiagen_1": "Q1", "Qiagen_2": "Q2"}

# Per panel: each tool's db_mode. NexVirome only has its own DB (same-DB in both).
PANELS = [
    ("same-DB", {"NexVirome": "same-DB", "Ganon": "same-DB", "Kraken2": "same-DB",
                 "Metabuli": "same-DB", "Phanta": "same-DB"}),
    ("native",  {"NexVirome": "same-DB", "Ganon": "native", "Kraken2": "native",
                 "Metabuli": "native", "Phanta": "native"}),
]


def panel(ax, piv, mode_by_tool, show_ylabel):
    """Grouped stacked bars: 5 tool groups, each with its 4 samples side by side."""
    n_s = len(SAMPLES)
    group_w, gap = 1.0, 0.1375  # halved again (tool groups closer together)
    bar_w = group_w / n_s * 0.9
    xticks, group_centers = [], []
    for gi, tool in enumerate(TOOLS):
        base = gi * (group_w + gap)
        mode = mode_by_tool[tool]
        for si, s in enumerate(SAMPLES):
            xx = base + si * (group_w / n_s)
            key = (tool, mode, s)
            bottom = 0.0
            for sp, col in zip(SPECIES, SP_COLORS):
                v = piv.loc[key, sp] if (key in piv.index and sp in piv.columns) else 0.0
                ax.bar(xx, v, bar_w, bottom=bottom, color=col, edgecolor="white",
                       linewidth=0.3, label=sp if (gi == 0 and si == 0) else None)
                bottom += v
            fp = piv.loc[key, "FP"] if (key in piv.index and "FP" in piv.columns) else 0.0
            ax.bar(xx, fp, bar_w, bottom=bottom, color=FP_COLOR, edgecolor="white",
                   linewidth=0.3, hatch="//",
                   label="FP (false positives)" if (gi == 0 and si == 0) else None)
            xticks.append(xx)
        group_centers.append((tool, base + (group_w / n_s) * (n_s - 1) / 2))
    for i in range(1, 6):
        ax.axhline(EXPECTED * i, color="grey", ls=":", lw=0.5, alpha=0.4, zorder=0)
    ax.set_ylim(0, 100)
    ax.set_xticks(xticks)
    # both the sample tick labels and the tool group labels = 15 (per request)
    ax.set_xticklabels([SAMPLE_LABEL[s] for s in SAMPLES] * len(TOOLS), fontsize=15)
    for tool, cx in group_centers:
        ax.text(cx, -0.20, tool, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=15, weight="bold")
    ax.tick_params(axis="x", length=0)
    if show_ylabel:
        ax.set_ylabel("Relative abundance (%)")  # rcParams axes.labelsize (22)


def main():
    d = pd.read_csv(CSV)
    m = d[d["level"] == "species"]
    piv = m.pivot_table(index=["tool", "db_mode", "sample"], columns="category",
                        values="TPM_pct").fillna(0)
    piv = piv.div(piv.sum(axis=1), axis=0) * 100  # renormalize each bar to 100

    # The plot box (axes) height must EQUAL one row of Fig2_tpfpf1_per_sample so the
    # two figures align when stacked. That row's measured axes height = 3.135 in.
    # We set the figure height and fix top/bottom in INCHES so the axes are exactly
    # AXES_H tall regardless of label/legend padding.
    AXES_H = 3.135          # measured axes height of one tpfpf1 row
    TOP_PAD = 0.45          # inches above axes (for the ylabel top / ticks)
    BOT_PAD = 1.35          # inches below axes (group labels + close legend)
    fig_h = AXES_H + TOP_PAD + BOT_PAD
    fig, axes = plt.subplots(1, 2, figsize=(18, fig_h), sharey=True)
    for ax, (mode, mode_by_tool) in zip(axes, PANELS):
        panel(ax, piv, mode_by_tool, show_ylabel=(mode == "same-DB"))

    # fix axes box: bottom/top as fractions giving exactly AXES_H inches
    fig.subplots_adjust(top=1 - TOP_PAD / fig_h, bottom=BOT_PAD / fig_h,
                        left=0.05, right=0.985, wspace=0.06)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=7, frameon=False,
               fontsize=14, bbox_to_anchor=(0.5, 0.01),
               columnspacing=1.0, handletextpad=0.4)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/fig2_relabund_bar_withFP.{ext}", dpi=300)
    print(f"-> {OUT}/fig2_relabund_bar_withFP.png , .pdf")


if __name__ == "__main__":
    main()
