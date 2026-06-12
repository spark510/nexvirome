#!/usr/bin/env python3
"""
Fig2 combined — *per-panel legend BELOW each panel* variant.

Same figure as render_fig2_combined.py (rows A=F1, B=TP/FP, C=composition, D=VAE;
columns same-DB | native), but each ROW carries its OWN legend in a thin strip
directly BELOW that row's panels, instead of one shared block at the very bottom.
The bottom-legend file (render_fig2_combined.py) is left untouched for comparison.

Legend placement — each legend sits just under the LEFT-column panel of its row,
spanning horizontally so it reads under both columns:
  A  — classifier colour key (5 tools)
  B  — TP(o)/FP(square) line key
  C  — 6 GT species (abbreviated) + FP block
  D  — classifier colour key (same as A)

Panel C uses a LOCAL wide-spacing bar layout (panel_wide) so the M1/M2/Q1/Q2
sample labels don't collide; the shared P2.panel is left unchanged.

Reuses f1/tpfp panel functions from render_fig2_tpfpf1_per_sample (P1) and the
species/colour constants from render_fig2_relabund_bar (P2).

Output: result_260605/fig2/Fig2_combined_inlegend.{png,eps}
Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260605/fig2/code/render_fig2_combined_inlegend.py
"""
from __future__ import annotations
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import render_fig2_tpfpf1_per_sample as P1
import render_fig2_relabund_bar as P2

mpl.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    # axis labels + ticks scaled 1.2x (16->19.2, 13->15.6); legend 1.5x (11->16.5)
    "font.size": 16, "axes.titlesize": 18, "axes.labelsize": 19.2,
    "xtick.labelsize": 15.6, "ytick.labelsize": 15.6, "legend.fontsize": 16.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.1,
})

NX = "/home/share/programs/nexvirome"
OUT = f"{NX}/result_260605/fig2"
SAMPLES = P1.SAMPLES
TOOL_ORDER = P1.TOOL_ORDER
TOOL_COLOR = P1.TOOL_COLOR
BENCHES = [("same-DB", "Reference-controlled DB"), ("native", "native DB")]

# abbreviated species names for panel C (full names are very long)
SP_ABBR = {
    "Adenovirus 40": "AdV-40", "HHV-5 (CMV)": "CMV", "Human RSV": "RSV",
    "Influenza B": "FluB", "Reovirus 3": "Reo-3", "Zika virus": "ZIKV",
}


def _vae_panel(ax, d, mode, show_ylabel):
    """VAE boxplot per tool — filled tool-colour box, black outline."""
    for ti, tool in enumerate(TOOL_ORDER):
        vals = d[(d["tool"] == tool) & (d["db_mode"] == mode)]["VAE"].dropna().values
        if len(vals) == 0:
            continue
        col = TOOL_COLOR[tool]
        ax.boxplot([vals], positions=[ti], widths=0.55, patch_artist=True,
                   showfliers=False, medianprops=dict(color="black", lw=2.4),
                   whiskerprops=dict(color="black", lw=1.8), capprops=dict(color="black", lw=1.8),
                   boxprops=dict(facecolor=col, edgecolor="black", lw=1.6, alpha=0.9))
        ax.scatter(np.full(len(vals), ti), vals, s=34, color=col,
                   edgecolor="black", linewidth=0.5, zorder=3)
    ax.set_xticks(range(len(TOOL_ORDER)))
    ax.set_xticklabels(TOOL_ORDER, rotation=0, ha="center")
    ax.set_xlim(-0.6, len(TOOL_ORDER) - 0.4); ax.set_ylim(0.10, 0.26)
    if show_ylabel:
        ax.set_ylabel("VAE")


def _tool_handles():
    return [Patch(facecolor=TOOL_COLOR[t], edgecolor="black", label=t)
            for t in TOOL_ORDER]


def panel_wide(ax, piv, mode_by_tool, show_ylabel):
    """Local copy of P2.panel with WIDER within-group sample spacing and bigger
    inter-group gaps so the M1/M2/Q1/Q2 labels don't collide. Shared P2.panel is
    left untouched."""
    SP, SPC, FPC = P2.SPECIES, P2.SP_COLORS, P2.FP_COLOR
    TOOLS, SAMP, SLAB = P2.TOOLS, P2.SAMPLES, P2.SAMPLE_LABEL
    n_s = len(SAMP)
    step = 0.66          # within-group spacing between samples (tighter: less
                         # whitespace between same-tool sample bars)
    bar_w = step * 0.86  # bar nearly fills the step -> minimal within-tool gaps
    group_w = step * (n_s - 1)
    gap = 0.9            # space between tool groups (was 0.1375)
    xticks, group_centers = [], []
    for gi, tool in enumerate(TOOLS):
        base = gi * (group_w + gap)
        mode = mode_by_tool[tool]
        for si, s in enumerate(SAMP):
            xx = base + si * step
            key = (tool, mode, s)
            bottom = 0.0
            for sp, col in zip(SP, SPC):
                v = piv.loc[key, sp] if (key in piv.index and sp in piv.columns) else 0.0
                ax.bar(xx, v, bar_w, bottom=bottom, color=col, edgecolor="white",
                       linewidth=0.3)
                bottom += v
            fp = piv.loc[key, "FP"] if (key in piv.index and "FP" in piv.columns) else 0.0
            ax.bar(xx, fp, bar_w, bottom=bottom, color=FPC, edgecolor="white",
                   linewidth=0.3, hatch="//")
            xticks.append(xx)
        group_centers.append((tool, base + group_w / 2))
    for i in range(1, 6):
        ax.axhline(P2.EXPECTED * i, color="grey", ls=":", lw=0.5, alpha=0.4, zorder=0)
    ax.set_ylim(0, 100)
    ax.set_xticks(xticks)
    ax.set_xticklabels([SLAB[s] for s in SAMP] * len(TOOLS), fontsize=15.6)
    for tool, cx in group_centers:
        ax.text(cx, -0.12, tool, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=18, weight="bold")
    ax.tick_params(axis="x", length=0)
    x_lo = -bar_w / 2 - 0.18
    x_hi = (len(TOOLS) - 1) * (group_w + gap) + group_w + bar_w / 2 + 0.18
    ax.set_xlim(x_lo, x_hi)
    if show_ylabel:
        ax.set_ylabel("Relative abundance (%)")


def _legend_row_center(fig, row_axes, handles, labels, ncol, gap=0.012, fs=16.5):
    """Legend centred horizontally on the FIGURE, in a thin strip directly below
    the given row's panels. `gap` is the figure-fraction offset below the lowest
    panel edge so the key sits close to the plots (no title)."""
    fig.canvas.draw()  # positions are valid after a draw
    y_bottom = min(ax.get_position().y0 for ax in row_axes)
    leg = fig.legend(handles, labels, loc="upper center",
                     bbox_to_anchor=(0.5, y_bottom - gap),
                     ncol=ncol, frameon=False, fontsize=fs, handlelength=1.3,
                     handletextpad=0.5, columnspacing=1.1, labelspacing=0.3,
                     borderaxespad=0.0)
    leg._legend_box.align = "center"
    fig.add_artist(leg)
    return leg


def _add_below_legends(fig, row0, row1, row2):
    """One titleless legend per row (A/B/C), centred on the figure just below that
    row. Row D carries no legend — its x-axis already labels each box by tool."""
    th = _tool_handles()
    tl = [h.get_label() for h in th]
    # A — classifier (extra gap so the key sits clear of the F1 bars)
    _legend_row_center(fig, row0, th, tl, ncol=5, gap=0.030)
    # B — TP / FP line key, then each tool as a colour-filled dot
    b_handles = [
        Line2D([0], [0], color=P1.TP_LINE, lw=2.6, marker="o", mfc="white",
               mec="black", ms=9, label="TP (left axis, linear)"),
        Line2D([0], [0], color=P1.FP_LINE, lw=2.6, marker="s", mfc="white",
               mec="black", ms=9, label="FP (right axis, log)"),
    ]
    b_handles += [
        Line2D([0], [0], marker="o", color="none", lw=0, mfc=TOOL_COLOR[t],
               mec="black", mew=0.8, ms=11, label=t)
        for t in TOOL_ORDER
    ]
    _legend_row_center(fig, row1, b_handles, [h.get_label() for h in b_handles],
                       ncol=7, gap=0.030)
    # C — species (abbreviated) + FP; sits BELOW the bold tool names (extra gap)
    sp_handles = [Patch(facecolor=c, edgecolor="white", label=SP_ABBR.get(s, s))
                  for s, c in zip(P2.SPECIES, P2.SP_COLORS)]
    sp_handles.append(Patch(facecolor=P2.FP_COLOR, edgecolor="white",
                            hatch="//", label="FP"))
    _legend_row_center(fig, row2, sp_handles, [h.get_label() for h in sp_handles],
                       ncol=7, gap=0.045)
    # D — no legend: the x-axis already labels each box by classifier.


def main():
    df1 = P1.load_data()
    d = pd.read_csv(P2.CSV)
    m = d[d["level"] == "species"]
    piv = m.pivot_table(index=["tool", "db_mode", "sample"], columns="category",
                        values="TPM_pct").fillna(0)
    piv = piv.div(piv.sum(axis=1), axis=0) * 100
    vae = pd.read_csv(P1.SRC)
    vae = vae[vae["sample"] != "MEAN"].copy()
    nv = vae[(vae.tool == "NexVirome") & (vae.db_mode == "same-DB")].copy()
    nv["db_mode"] = "native"
    vae = pd.concat([vae, nv], ignore_index=True)

    fig = plt.figure(figsize=(16, 16.5))
    # tighter hspace: panels closer together; the per-row legend still fits in the
    # strip below each row
    gs = GridSpec(4, 2, figure=fig, hspace=0.50, wspace=0.16,
                  height_ratios=[1.0, 1.0, 1.1, 1.0])

    rows = {0: [], 1: [], 2: [], 3: []}   # row -> [left ax, right ax]
    first_col = {}
    for ci, (bench, title) in enumerate(BENCHES):
        ax_f1 = fig.add_subplot(gs[0, ci])
        P1.f1_panel(ax_f1, df1, bench, show_ylabel=(ci == 0))
        ax_f1.set_title(title, fontsize=18, weight="bold")

        ax_tp = fig.add_subplot(gs[1, ci])
        P1.tpfp_panel(ax_tp, df1, bench,
                      show_tp_label=(ci == 0), show_fp_label=(ci == 1))

        ax_rb = fig.add_subplot(gs[2, ci])
        mode_by_tool = dict(P2.PANELS)[bench]
        panel_wide(ax_rb, piv, mode_by_tool, show_ylabel=(ci == 0))

        ax_vae = fig.add_subplot(gs[3, ci])
        _vae_panel(ax_vae, vae, bench, show_ylabel=(ci == 0))

        for r, a in zip((0, 1, 2, 3), (ax_f1, ax_tp, ax_rb, ax_vae)):
            rows[r].append(a)
        if ci == 0:
            first_col = {0: ax_f1, 1: ax_tp, 2: ax_rb, 3: ax_vae}

    for row, lab in [(0, "A"), (1, "B"), (2, "C"), (3, "D")]:
        first_col[row].text(-0.16, 1.10, lab, transform=first_col[row].transAxes,
                            fontsize=26, weight="bold", va="bottom")

    fig.subplots_adjust(left=0.08, right=0.97, top=0.97, bottom=0.05)
    # per-row legends, titleless and centred on the figure under each row (A/B/C)
    _add_below_legends(fig, rows[0], rows[1], rows[2])

    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/Fig2_combined_inlegend.{ext}", dpi=300, bbox_inches="tight")
    print(f"-> {OUT}/Fig2_combined_inlegend.png , .eps")


if __name__ == "__main__":
    main()
