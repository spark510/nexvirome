#!/usr/bin/env python3
"""
Fig2 metric plot (FINALIZED 2026-06-03) — 2 rows x 2 cols.

Row 1 (F1):  per-sample grouped BARS, 4 KIT samples x 5 tools (unchanged).
Row 2 (TP+FP combined): x = 5 tools, TWIN y-axis —
    left  = TP (species, 0-6)  ONE green line across the tools
    right = FP (species, log)  ONE red line across the tools
  Each point = mean over the 4 KIT samples; the 4 per-sample values are shown as
  faint jittered dots. So exactly two lines per panel: TP (green) and FP (red).

Columns: (a) same-DB / (b) native. NexVirome's viral-only DB is its native, so
its same-DB values are copied into the native column.

Input : result_260605/fig2/fig2_kit_scores.csv  (sp_TP, sp_FP, sp_F1 per sample)
Output: result_260605/fig2/Fig2_tpfpf1_per_sample.{png,pdf}
Run: conda run -n shotgun_virome python result_260605/fig2/code/render_fig2_tpfpf1_per_sample.py
"""
from __future__ import annotations
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl

NX = "/home/share/programs/nexvirome"
SRC = f"{NX}/result_260605/fig2/fig2_kit_scores.csv"
OUT = f"{NX}/result_260605/fig2"

mpl.rcParams.update({
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

SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
# display label only (the MagNA extraction kit is Roche MagNA Pure); the internal
# sample id MagNA_1 stays unchanged in the CSVs/pipeline.
SAMPLE_DISPLAY = {"MagNA_1": "Roche 1", "MagNA_2": "Roche 2",
                  "Qiagen_1": "Qiagen 1", "Qiagen_2": "Qiagen 2"}
TOOL_ORDER = ["NexVirome", "Ganon", "Phanta", "Kraken2", "Metabuli"]
TOOL_COLOR = {
    "NexVirome": "#1D4E89", "Ganon": "#00B2CA", "Phanta": "#7DCFB6",
    "Kraken2": "#FF595E", "Metabuli": "#FBD1A2",
}
TP_LINE = "black"      # TP connector line (dots stay tool-coloured)
FP_LINE = "#de2d26"    # FP connector line (solid)
F1_LABEL_FS = 12       # F1 bar value labels (1.2x of the old 10)
LABEL_FS = 14          # TP/FP per-dot value labels (1.4x of the old 10)


def load_data():
    df = pd.read_csv(SRC)
    df = df[df["sample"] != "MEAN"].copy()
    df = df[df["tool"].isin(TOOL_ORDER)].copy()
    df = df.rename(columns={"db_mode": "bench"})
    nv = df[(df["tool"] == "NexVirome") & (df["bench"] == "same-DB")].copy()
    nv["bench"] = "native"
    return pd.concat([df, nv], ignore_index=True)


def f1_panel(ax, df, bench, show_ylabel):
    """Row 1: per-sample F1 bars (x = 4 samples, bars per tool)."""
    sub = df[df["bench"] == bench]
    n = len(TOOL_ORDER)
    width = 0.95 / n * 0.85
    x = np.arange(len(SAMPLES))
    for i, tool in enumerate(TOOL_ORDER):
        vals = [float(sub[(sub.tool == tool) & (sub["sample"] == s)]["sp_F1"].iloc[0])
                if not sub[(sub.tool == tool) & (sub["sample"] == s)].empty else 0.0
                for s in SAMPLES]
        pos = x + (i - n / 2 + 0.5) * width
        ax.bar(pos, vals, width, color=TOOL_COLOR[tool], edgecolor="black",
               linewidth=0.6, label=tool if bench == "same-DB" else None)
        for px, v in zip(pos, vals):
            ax.text(px, v + 0.02, f"{v:.2f}", ha="center", va="bottom", fontsize=F1_LABEL_FS)
    ax.set_xticks(x); ax.set_xticklabels([SAMPLE_DISPLAY[s] for s in SAMPLES])
    ax.set_xlim(-0.5, len(SAMPLES) - 0.5)   # align x with the TP/FP line panel below
    ax.set_ylim(0, 1.08)
    if show_ylabel:
        ax.set_ylabel("F1 score")


def tpfp_panel(ax_tp, df, bench, show_tp_label, show_fp_label):
    """Row 2: the 5 tools are laid out left-to-right WITHIN each sample block; a
    line connects those 5 tool dots inside the block. One block per sample
    (Roche 1 / Roche 2 / Qiagen 1 / Qiagen 2). Each block has a green TP line (left
    axis) and a red FP line (right axis). Twin y. Every dot is labelled."""
    sub = df[df["bench"] == bench]
    ax_fp = ax_tp.twinx()
    ax_fp.spines["top"].set_visible(False)

    # Use the SAME x geometry as the F1 bar panel above so the 5 dots sit exactly
    # over the 5 bars, and sample blocks line up 1:1. F1 bars: x=arange(4 samples),
    # bar i at x + (i - n/2 + 0.5)*width, width = 0.95/n*0.85.
    n_tool = len(TOOL_ORDER)
    width = 0.95 / n_tool * 0.85
    tool_dx = np.array([(i - n_tool / 2 + 0.5) * width for i in range(n_tool)])
    centers, labels = [], []
    for bi, s in enumerate(SAMPLES):
        base = float(bi)                 # same as F1 panel x = arange(samples)
        xs = base + tool_dx
        tp = [float(sub[(sub.tool == t) & (sub["sample"] == s)]["sp_TP"].iloc[0])
              if not sub[(sub.tool == t) & (sub["sample"] == s)].empty else np.nan
              for t in TOOL_ORDER]
        fp = [float(sub[(sub.tool == t) & (sub["sample"] == s)]["sp_FP"].iloc[0])
              if not sub[(sub.tool == t) & (sub["sample"] == s)].empty else np.nan
              for t in TOOL_ORDER]
        # connect the 5 tool dots within this sample block: TP line black,
        # FP line solid red; dots coloured per tool to match the F1 bars above.
        ax_tp.plot(xs, tp, "-", color=TP_LINE, lw=2.4, zorder=3)
        ax_fp.plot(xs, [v + 1 for v in fp], "-", color=FP_LINE, lw=2.0, zorder=2)
        for xi, v, t in zip(xs, tp, TOOL_ORDER):
            ax_tp.scatter(xi, v, s=70, color=TOOL_COLOR[t], edgecolor="black",
                          linewidth=0.6, zorder=4)
        for xi, v, t in zip(xs, fp, TOOL_ORDER):
            ax_fp.scatter(xi, v + 1, s=60, marker="s", color=TOOL_COLOR[t],
                          edgecolor="black", linewidth=0.6, zorder=3)
        for xi, v in zip(xs, tp):
            if not np.isnan(v):
                ax_tp.annotate(f"{int(round(v))}", (xi, v), textcoords="offset points",
                               xytext=(0, 7), ha="center", fontsize=LABEL_FS, color="black")
        for xi, v in zip(xs, fp):
            if not np.isnan(v):
                # FP label sits directly above its square. The FP log-axis top is
                # raised (set_ylim below) so even the largest FP points stay well
                # under the TP=6 dots, giving the FP label room above the square
                # without colliding with the TP "6" label.
                ax_fp.annotate(f"{int(round(v))}", (xi, v + 1), textcoords="offset points",
                               xytext=(0, 8), ha="center", va="bottom",
                               fontsize=LABEL_FS, color=FP_LINE)
        centers.append(base); labels.append(SAMPLE_DISPLAY[s])

    ax_tp.set_ylim(0, 6.5); ax_tp.set_yticks(range(0, 7))
    # raise the FP log-axis top so the largest FP points (~800) sit well below the
    # TP=6 dots, leaving headroom for the FP value label above each square
    ax_fp.set_yscale("log"); ax_fp.set_ylim(1, 8000)
    ax_tp.set_xticks(centers); ax_tp.set_xticklabels(labels)
    ax_tp.set_xlim(-0.5, len(SAMPLES) - 0.5)   # match the F1 bar panel above
    if show_tp_label:
        ax_tp.set_ylabel("TP (species)", color="black")
    ax_tp.tick_params(axis="y", colors="black")
    ax_tp.spines["left"].set_color("black")
    if show_fp_label:
        ax_fp.set_ylabel("FP (species) : log-scale", color=FP_LINE)
    ax_fp.tick_params(axis="y", colors=FP_LINE)
    ax_fp.spines["right"].set_color(FP_LINE)


def main():
    df = load_data()
    os.makedirs(OUT, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(18, 9.2))
    for col, bench in enumerate(["same-DB", "native"]):
        f1_panel(axes[0, col], df, bench, show_ylabel=(col == 0))
        tpfp_panel(axes[1, col], df, bench, show_tp_label=(col == 0), show_fp_label=(col == 1))
        axes[0, col].text(-0.02, 1.15, "(a)" if col == 0 else "(b)",
                          transform=axes[0, col].transAxes, ha="left", va="bottom",
                          fontsize=26, weight="bold")
        axes[0, col].set_title("same-DB" if col == 0 else "native")

    # legend: F1 tool bars (row1) + TP/FP line key (row2)
    h1, l1 = axes[0, 0].get_legend_handles_labels()
    from matplotlib.lines import Line2D
    h2 = [Line2D([0], [0], color=TP_LINE, lw=3, marker="o", label="TP (left axis)"),
          Line2D([0], [0], color=FP_LINE, lw=3, marker="s", label="FP (right axis)")]
    fig.legend(h1 + h2, l1 + ["TP (left axis)", "FP (right axis)"],
               loc="lower center", ncol=7, frameon=False, fontsize=16,
               bbox_to_anchor=(0.5, -0.01))
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    for ext in ("png", "eps"):
        out = f"{OUT}/Fig2_tpfpf1_per_sample.{ext}"
        fig.savefig(out, dpi=300)
        print(f"saved {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
