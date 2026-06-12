#!/usr/bin/env python3
"""
Fig2 VAE plot — abundance accuracy vs the ideal equal mock (VAE = 0).

VAE = mean |observed_TPM_frac - 1/6| over the 6 GT species (TPM-normalized,
GT-only). VAE=0 means a perfectly even 6-species composition; larger = more
skewed. A dotted line at y=0 marks the ideal.

x = 5 tools; each tool shows two points: same-DB (open) and native (filled).
NexVirome's viral-only DB is its native, so its native point = its same-DB value.
The 4 per-sample VAE values are drawn as small jittered dots; the mean as a big
marker. Lower is better.

Input : result_260605/fig2/fig2_kit_scores.csv  (VAE column, per sample)
Output: result_260605/fig2/fig2_vae.{png,pdf}
Run: conda run -n shotgun_virome python result_260605/fig2/code/render_fig2_vae.py
"""
from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Unified style — identical to Fig2_tpfpf1_per_sample / relabund (same Fig2 stack).
matplotlib.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42, "ps.fonttype": 42,
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
CSV = f"{NX}/result_260605/fig2/fig2_kit_scores.csv"
OUT = f"{NX}/result_260605/fig2"

TOOLS = ["NexVirome", "Ganon", "Kraken2", "Metabuli", "Phanta"]
TOOL_COLOR = {
    "NexVirome": "#1D4E89", "Ganon": "#00B2CA", "Kraken2": "#FF595E",
    "Metabuli": "#FBD1A2", "Phanta": "#7DCFB6",
}
# same-DB = open circle, native = filled circle
MODE_STYLE = {"same-DB": dict(marker="o", facecolor="white", label="same-DB"),
              "native":  dict(marker="o", facecolor=None,    label="native")}


def main():
    d = pd.read_csv(CSV)
    d = d[d["sample"] != "MEAN"].copy()
    # NexVirome native = same-DB (its viral-only DB is its native)
    nv = d[(d["tool"] == "NexVirome") & (d["db_mode"] == "same-DB")].copy()
    nv["db_mode"] = "native"
    d = pd.concat([d, nv], ignore_index=True)

    # plot box (axes) height = one row of Fig2_tpfpf1_per_sample (3.135 in) so the
    # figures align when stacked; x/y fonts already follow the shared rcParams.
    AXES_H = 3.135
    TOP_PAD = 0.55          # title + headroom
    BOT_PAD = 0.75          # tool labels (no legend in this figure)
    fig_h = AXES_H + TOP_PAD + BOT_PAD
    ymax = max(0.26, d["VAE"].max() * 1.18)
    fig, axes = plt.subplots(1, 2, figsize=(18, fig_h), sharey=True)

    PANELS = [("same-DB", "(a) same-DB"), ("native", "(b) native")]
    for ax, (mode, title) in zip(axes, PANELS):
        x = np.arange(len(TOOLS))
        for ti, tool in enumerate(TOOLS):
            vals = d[(d["tool"] == tool) & (d["db_mode"] == mode)]["VAE"].dropna().values
            if len(vals) == 0:
                continue
            col = TOOL_COLOR[tool]
            # unfilled box, thicker lines
            ax.boxplot([vals], positions=[ti], widths=0.55, patch_artist=True,
                       showfliers=False, medianprops=dict(color=col, lw=3.0),
                       whiskerprops=dict(color=col, lw=2.6),
                       capprops=dict(color=col, lw=2.6),
                       boxprops=dict(facecolor="none", edgecolor=col, lw=2.8))
            ax.scatter(np.full(len(vals), ti), vals, s=42, color=col,
                       edgecolor="white", linewidth=0.4, zorder=3)
        # ideal reference line — grey dotted across the panel
        ax.axhline(0, color="grey", ls=":", lw=1.6, zorder=0)
        ax.text(len(TOOLS) - 0.5, 0.006, "ideal (VAE = 0)", ha="right", va="bottom",
                color="grey", fontsize=15)
        ax.set_xticks(x); ax.set_xticklabels(TOOLS)  # horizontal, rcParams xtick size (15)
        ax.set_xlim(-0.6, len(TOOLS) - 0.4)
        ax.set_ylim(-0.015, ymax)
        ax.set_title(title)
        ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    axes[0].set_ylabel("VAE")
    # fix axes box to exactly AXES_H inches (match the stack)
    fig.subplots_adjust(top=1 - TOP_PAD / fig_h, bottom=BOT_PAD / fig_h,
                        left=0.05, right=0.985, wspace=0.06)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/fig2_vae.{ext}", dpi=300)
    print(f"-> {OUT}/fig2_vae.png , .pdf")


if __name__ == "__main__":
    main()
