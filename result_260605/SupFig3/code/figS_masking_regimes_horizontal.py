#!/usr/bin/env python3
"""
Supplementary figure (horizontal layout): masking regime effect on detection
metrics, LCA only.

  rows = level (species top, genus bottom)
  cols = 6 metrics (TP, FP, FN, Precision, Recall, F1)
  per panel: 3 regime bars (white -> grey -> black gradient) + 4 sample dots

All strict greyscale: bar fills white / grey / black,
outlines black, dot edges black, no transparency.

Run:
  conda run -n shotgun_virome python \
      scripts/benchmark/figS_masking_regimes_horizontal.py
Outputs:
  paper/figures/FigSx_masking_regimes_horizontal.{png,tiff,eps}
  paper/tables/masking_regimes_per_sample.csv      (re-used)
  paper/tables/masking_regimes_summary.csv         (re-used)
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, "/home/share/programs/nexvirome/notebooks")
sys.path.insert(0, "/home/share/programs/nexvirome/scripts/benchmark")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import paper_style  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from benchmark_utils import (GROUND_TRUTH, evaluate_sample, evaluate_at_rank,
                             parse_kreport, kreport_to_species_counts,
                             save_figure)

GROUND_TRUTH_GENUS = [
    {'name': g['name'], 'taxids': {g['genus']}, 'expected_frac': g['expected_frac']}
    for g in GROUND_TRUTH
]
NX = "/home/share/programs/nexvirome"
# GOLDEN method-B run output (run_masking_regimes.py): /tmp/masking_regimes_golden
ROOT = "/tmp/masking_regimes_golden"
SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
REGIMES = ["No mask", "Pre-mask (DB-level)", "Post-mask (full)"]
REGIME_DIRS = {"No mask": "No_mask",
               "Pre-mask (DB-level)": "Pre-mask_DB-level",
               "Post-mask (full)": "Post-mask_full"}
REGIME_LABELS = ["No mask", "Pre-mask", "Post-mask"]
MODE = "method B"          # GOLDEN single mode (was LCA)
MIN_READS = 3              # GOLDEN per-taxon read floor (was 1)

METRICS = [          # (col suffix, display name, is_fraction)
    ("TP",        "TP",        False),
    ("FP",        "FP",        False),
    # FN dropped: TP is 6/6 in every regime, so FN is uniformly 0 (empty panel).
    ("Precision", "Precision", True),
    ("Recall",    "Recall",    True),
    ("F1",        "F1",        True),
]
LEVELS = [("species", "Species level"), ("genus", "Genus level")]

# Sample = colour (constant greyscale fill per sample, no hatch).
# x-axis primary category = masking regime; the 3 regime groups are separated
# by a visual gap so the no/pre/post grouping is obvious at a glance.
SAMPLE_FILL = {"MagNA_1":  "#000000",   # solid black
               "MagNA_2":  "#666666",   # dark grey
               "Qiagen_1": "#bdbdbd",   # light grey
               "Qiagen_2": "#ffffff"}   # white


def collect():
    rows = []
    for rg in REGIMES:
        for s in SAMPLES:
            kp = f"{ROOT}/{REGIME_DIRS[rg]}/{MODE}/{s}/{s}.kreport"
            if not os.path.exists(kp): continue
            kdf = parse_kreport(kp)
            c = kreport_to_species_counts(kdf)
            ev_sp = evaluate_sample(c, GROUND_TRUTH, min_reads=MIN_READS)
            ev_g  = evaluate_at_rank(kdf, GROUND_TRUTH_GENUS,
                                     rank_code='G', min_reads=MIN_READS)
            rows.append(dict(regime=rg, sample=s,
                             TP_species=ev_sp["TP"], FP_species=ev_sp["FP"],
                             FN_species=ev_sp["FN"],
                             Precision_species=ev_sp["precision"],
                             Recall_species=ev_sp["recall"],
                             F1_species=ev_sp["f1"],
                             TP_genus=ev_g["TP"], FP_genus=ev_g["FP"],
                             FN_genus=ev_g["FN"],
                             Precision_genus=ev_g["precision"],
                             Recall_genus=ev_g["recall"],
                             F1_genus=ev_g["f1"]))
    return pd.DataFrame(rows)


def render(df):
    # 2 rows (species, genus) x 6 cols (TP,FP,FN,Precision,Recall,F1)
    ncol = len(METRICS)
    fig, axes = plt.subplots(2, ncol, figsize=(3.25 * ncol, 8.5), squeeze=False,
                             gridspec_kw=dict(hspace=0.40, wspace=0.32,
                                              left=0.06, right=0.985,
                                              top=0.90, bottom=0.16))
    # FP/FN can be large; TP is bounded by the 6 GT species, so give TP its own
    # 0–6 axis instead of the shared FP/FN scale (else TP bars look squashed).
    n_gt = len(GROUND_TRUTH)                       # 6 ground-truth species
    fpfn_max = max(df[["FP_species", "FP_genus"]].max().max(), 1)  # FP scale (FN dropped)

    # x layout: 4 sample bars per regime group, separated by a visible gap.
    n_s    = len(SAMPLES)
    bar_w  = 1.0
    gap    = 1.5                                 # blank space between regimes
    group_w = n_s * bar_w + gap
    # x-coord of the j-th sample bar inside regime k (0..n_s-1, 0..2)
    bar_x  = lambda k, j: k * group_w + j * bar_w
    # centre of each regime group (for x ticks / labels)
    group_centres = [k * group_w + (n_s - 1) * bar_w / 2 for k in range(len(REGIMES))]
    xmin = -bar_w / 2 - 0.4
    xmax = group_centres[-1] + (n_s - 1) * bar_w / 2 + bar_w / 2 + 0.4

    for r, (lvl, lvl_label) in enumerate(LEVELS):
        for c, (suf, mname, is_frac) in enumerate(METRICS):
            ax = axes[r, c]
            metric = f"{suf}_{lvl}"
            for k, rg in enumerate(REGIMES):
                for j, s in enumerate(SAMPLES):
                    v = df[(df.regime == rg) & (df["sample"] == s)][metric].mean()
                    ax.bar(bar_x(k, j), v, width=bar_w * 0.92,
                           facecolor=SAMPLE_FILL[s],
                           edgecolor="black", linewidth=0.9, zorder=2)

            ax.set_xlim(xmin, xmax)
            ax.set_xticks(group_centres)
            ax.set_xticklabels(REGIME_LABELS, fontsize=13,
                               rotation=20, ha="right",
                               rotation_mode="anchor")
            if is_frac:
                ax.set_ylim(0, 1.08)
                ax.set_yticks(np.arange(0, 1.05, 0.2))
            elif suf == "TP":
                ax.set_ylim(0, n_gt * 1.10)        # TP bounded by 6 GT species
                ax.set_yticks(range(0, n_gt + 1, 2))
            else:
                ax.set_ylim(0, fpfn_max * 1.18)    # FP / FN share their own scale
            ax.tick_params(axis="y", labelsize=11)
            ax.grid(axis="y", linestyle=":", linewidth=0.6,
                    color="black", alpha=0.30, zorder=0)
            for nm, sp in ax.spines.items():
                if nm in ("top", "right"):
                    sp.set_visible(False)
                else:
                    sp.set_visible(True); sp.set_color("black"); sp.set_linewidth(1.0)
            ax.tick_params(axis="both", colors="black", width=1.0)

            # column header: metric name (only top row)
            if r == 0:
                ax.set_title(mname, fontsize=16, fontweight="bold", pad=8)

    # Big row labels on the left margin (rotated)
    fig.text(0.012, 0.685, "Species level", rotation=90,
             ha="center", va="center", fontsize=17, fontweight="bold")
    fig.text(0.012, 0.270, "Genus level",   rotation=90,
             ha="center", va="center", fontsize=17, fontweight="bold")

    # Light horizontal divider between species and genus rows
    sep_y = (axes[0, 0].get_position().y0 + axes[1, 0].get_position().y1) / 2
    fig.add_artist(plt.Line2D([0.04, 0.99], [sep_y, sep_y],
                              color="black", linewidth=0.7,
                              transform=fig.transFigure))

    # Legend: 4 sample fills (sample = colour, the only legend axis needed)
    handles = [Patch(facecolor=SAMPLE_FILL[s], edgecolor="black",
                     linewidth=1.0, label=s)
               for s in SAMPLES]
    fig.legend(handles=handles, loc="lower center", ncol=4,
               fontsize=14, frameon=False, bbox_to_anchor=(0.5, 0.005))

    return fig


def main():
    print("collecting…")
    df = collect()
    print(df.to_string(index=False))
    print("\nrendering…")
    fig = render(df)
    save_figure(fig, "FigSx_masking_regimes_horizontal", close=True)


if __name__ == "__main__":
    main()
