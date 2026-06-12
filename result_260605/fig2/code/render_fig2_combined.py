#!/usr/bin/env python3
"""
Fig2 combined — single stacked figure (A / B / C) reusing the three standalone
panel functions so the published Fig2 is one file, with unified Arial style,
shared tool palette and aligned axes.

  (A) Detection performance  — per-sample F1 bars (row 1) + TP/FP twin-axis line
                               (row 2), columns = same-DB | native
  (B) Community composition  — TPM stacked bars, 6 GT species + FP, same-DB|native
  (C) Abundance accuracy     — VAE boxplots vs ideal=0, same-DB | native

Panel functions are imported from:
  render_fig2_tpfpf1_per_sample  (load_data, f1_panel, tpfp_panel, TOOL_*)
  render_fig2_relabund_bar       (panel as relabund_panel, load + SP_COLORS)
  render_fig2_vae                (data load inline; reuse TOOL_COLOR)
Output: result_260605/fig2/Fig2_combined.{png,eps}
Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260605/fig2/code/render_fig2_combined.py
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

# unified style (matches the standalone scripts)
mpl.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 16, "axes.titlesize": 18, "axes.labelsize": 16,
    "xtick.labelsize": 13, "ytick.labelsize": 13, "legend.fontsize": 13,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.1,
})

NX = "/home/share/programs/nexvirome"
OUT = f"{NX}/result_260605/fig2"
SAMPLES = P1.SAMPLES
TOOL_ORDER = P1.TOOL_ORDER
TOOL_COLOR = P1.TOOL_COLOR
BENCHES = [("same-DB", "same-DB"), ("native", "native")]   # column subtitles (no a/b)


def _vae_panel(ax, d, mode, show_ylabel):
    """VAE boxplot per tool for one bench column (reuses TOOL_COLOR)."""
    for ti, tool in enumerate(TOOL_ORDER):
        vals = d[(d["tool"] == tool) & (d["db_mode"] == mode)]["VAE"].dropna().values
        if len(vals) == 0:
            continue
        col = TOOL_COLOR[tool]
        # filled box in the tool colour, black outline / median / whiskers
        ax.boxplot([vals], positions=[ti], widths=0.55, patch_artist=True,
                   showfliers=False, medianprops=dict(color="black", lw=2.4),
                   whiskerprops=dict(color="black", lw=1.8), capprops=dict(color="black", lw=1.8),
                   boxprops=dict(facecolor=col, edgecolor="black", lw=1.6, alpha=0.9))
        ax.scatter(np.full(len(vals), ti), vals, s=34, color=col,
                   edgecolor="black", linewidth=0.5, zorder=3)
    ax.axhline(0, color="grey", ls=":", lw=1.4, zorder=0)
    ax.set_xticks(range(len(TOOL_ORDER)))
    ax.set_xticklabels(TOOL_ORDER, rotation=20, ha="right")
    ax.set_xlim(-0.6, len(TOOL_ORDER) - 0.4); ax.set_ylim(-0.015, 0.26)
    if show_ylabel:
        ax.set_ylabel("VAE")


def main():
    df1 = P1.load_data()
    # relabund inputs (replicate P2.main()'s pivot; P2 has no load() helper)
    d = pd.read_csv(P2.CSV)
    m = d[d["level"] == "species"]
    piv = m.pivot_table(index=["tool", "db_mode", "sample"], columns="category",
                        values="TPM_pct").fillna(0)
    piv = piv.div(piv.sum(axis=1), axis=0) * 100
    # vae data
    vae = pd.read_csv(P1.SRC)
    vae = vae[vae["sample"] != "MEAN"].copy()
    nv = vae[(vae.tool == "NexVirome") & (vae.db_mode == "same-DB")].copy()
    nv["db_mode"] = "native"
    vae = pd.concat([vae, nv], ignore_index=True)

    # layout: 4 content rows (A-f1, A-tpfp, B-relabund, C-vae) x 2 bench cols
    fig = plt.figure(figsize=(16, 16))
    gs = GridSpec(4, 2, figure=fig, hspace=0.55, wspace=0.16,
                  height_ratios=[1.0, 1.0, 1.1, 1.0])

    first_col = {}   # row -> left-column ax, for panel-letter placement
    for ci, (bench, title) in enumerate(BENCHES):
        ax_f1 = fig.add_subplot(gs[0, ci])
        P1.f1_panel(ax_f1, df1, bench, show_ylabel=(ci == 0))
        ax_f1.set_title(title, fontsize=18, weight="bold")

        ax_tp = fig.add_subplot(gs[1, ci])
        P1.tpfp_panel(ax_tp, df1, bench,
                      show_tp_label=(ci == 0), show_fp_label=(ci == 1))

        ax_rb = fig.add_subplot(gs[2, ci])
        mode_by_tool = dict(P2.PANELS)[bench]   # per-bench tool->db_mode map
        P2.panel(ax_rb, piv, mode_by_tool, show_ylabel=(ci == 0))
        # match panel B's plotting width: B uses 12.5 % whitespace each side of
        # its bar block; replicate the same fractional margin here so the C bars
        # don't run edge-to-edge while B sits inset.
        n_s, n_t = len(P2.SAMPLES), len(P2.TOOLS)
        group_w, gap, bar_w = 1.0, 0.1375, 1.0 / len(P2.SAMPLES) * 0.9
        x_lo = 0.0 - bar_w / 2
        x_hi = (n_t - 1) * (group_w + gap) + (n_s - 1) * (group_w / n_s) + bar_w / 2
        pad = 0.125 * (x_hi - x_lo)
        ax_rb.set_xlim(x_lo - pad, x_hi + pad)

        ax_vae = fig.add_subplot(gs[3, ci])
        _vae_panel(ax_vae, vae, bench, show_ylabel=(ci == 0))

        if ci == 0:
            first_col = {0: ax_f1, 1: ax_tp, 2: ax_rb, 3: ax_vae}

    # one letter per row: A=F1, B=TP/FP, C=composition, D=VAE
    for row, lab in [(0, "A"), (1, "B"), (2, "C"), (3, "D")]:
        first_col[row].text(-0.16, 1.10, lab, transform=first_col[row].transAxes,
                            fontsize=26, weight="bold", va="bottom")

    # ---- per-panel legends, all stacked at the very bottom ----
    # abbreviated species names for panel C (full names are very long)
    SP_ABBR = {
        "Adenovirus 40": "AdV-40", "HHV-5 (CMV)": "CMV", "Human RSV": "RSV",
        "Influenza B": "FluB", "Reovirus 3": "Reo-3", "Zika virus": "ZIKV",
    }

    # (A, D) classifier colour key — A bars and D boxes are coloured by tool
    tool_handles = [Patch(facecolor=TOOL_COLOR[t], edgecolor="black", label=t)
                    for t in TOOL_ORDER]
    leg_ad = fig.legend(tool_handles, [h.get_label() for h in tool_handles],
                        loc="lower center", ncol=len(TOOL_ORDER), frameon=False,
                        bbox_to_anchor=(0.5, 0.040), fontsize=13,
                        title="A, D — classifier", title_fontsize=12)
    leg_ad._legend_box.align = "center"
    fig.add_artist(leg_ad)

    # (B) TP / FP line key — line colour = metric, marker shape = TP(○)/FP(□),
    #     dot fill = tool (same colours as A/D)
    b_handles = [
        Line2D([0], [0], color=P1.TP_LINE, lw=2.6, marker="o", mfc="white",
               mec="black", ms=9, label="TP (left axis, linear)"),
        Line2D([0], [0], color=P1.FP_LINE, lw=2.6, marker="s", mfc="white",
               mec="black", ms=9, label="FP (right axis, log)"),
        Line2D([0], [0], marker="o", color="0.4", lw=0, mfc="0.4", mec="black",
               ms=9, label="dot fill = classifier (see A, D)"),
    ]
    leg_b = fig.legend(b_handles, [h.get_label() for h in b_handles],
                       loc="lower center", ncol=3, frameon=False,
                       bbox_to_anchor=(0.5, -0.005), fontsize=13,
                       title="B — TP / FP", title_fontsize=12)
    leg_b._legend_box.align = "center"
    fig.add_artist(leg_b)

    # (C) mock community species (abbreviated) + combined false-positive block
    sp_handles = [Patch(facecolor=c, edgecolor="white",
                        label=SP_ABBR.get(s, s) + f"  ({s})")
                  for s, c in zip(P2.SPECIES, P2.SP_COLORS)]
    sp_handles.append(Patch(facecolor=P2.FP_COLOR, edgecolor="white",
                            hatch="//", label="FP  (false positives, combined)"))
    leg_c = fig.legend(sp_handles, [h.get_label() for h in sp_handles],
                       loc="lower center", ncol=4, frameon=False,
                       bbox_to_anchor=(0.5, -0.075), fontsize=12,
                       title="C — mock-community ground-truth species", title_fontsize=12)
    leg_c._legend_box.align = "center"

    fig.subplots_adjust(left=0.08, right=0.97, top=0.96, bottom=0.135)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/Fig2_combined.{ext}", dpi=300, bbox_inches="tight")
    print(f"-> {OUT}/Fig2_combined.png , .eps")


if __name__ == "__main__":
    main()
