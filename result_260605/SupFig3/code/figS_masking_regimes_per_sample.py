#!/usr/bin/env python3
"""
Supplementary figure: per-sample breakdown of the masking regime comparison.

Main Fig3 reports mean F1/FP across the 4 KIT samples. This Supp figure shows
the underlying per-sample dispersion: each sample contributes one dot, bars
show the mean. Two side-by-side panels (F1, FP) × three regimes × three modes.

Companion table (Table Sx) lists exact F1/FP mean ± SD across samples.

Inputs (existing artefacts from `run_masking_regimes.py`):
  /tmp/masking_regimes/<regime>/<mode>/<sample>/<sample>.kreport
  ground truth from benchmark_utils.GROUND_TRUTH

Run:
  conda run -n shotgun_virome python scripts/benchmark/figS_masking_regimes_per_sample.py
Outputs:
  paper/figures/FigSx_masking_regimes_per_sample.{png,tiff,eps}
  paper/tables/masking_regimes_per_sample.csv
  paper/tables/masking_regimes_summary.csv  (mean +/- SD)
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, "/home/share/programs/nexvirome/notebooks")
sys.path.insert(0, "/home/share/programs/nexvirome/scripts/benchmark")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import paper_style  # noqa: F401  (Arial style)
import matplotlib.pyplot as plt
from benchmark_utils import (GROUND_TRUTH, evaluate_sample, evaluate_at_rank,
                             parse_kreport, kreport_to_species_counts,
                             save_figure)

# Genus-level ground truth: same 6 viruses, taxids replaced by genus taxid.
# evaluate_at_rank() at rank_code='G' counts kreport rows at G level against
# the 'taxids' set of each GT group — so we put the genus taxid here.
GROUND_TRUTH_GENUS = [
    {'name': g['name'], 'taxids': {g['genus']}, 'expected_frac': g['expected_frac']}
    for g in GROUND_TRUTH
]

NX = "/home/share/programs/nexvirome"
# GOLDEN method-B run output (run_masking_regimes.py)
ROOT = "/tmp/masking_regimes_golden"
SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
REGIMES = ["No mask", "Pre-mask (DB-level)", "Post-mask (full)"]
REGIME_DIRS = {
    "No mask": "No_mask",
    "Pre-mask (DB-level)": "Pre-mask_DB-level",
    "Post-mask (full)": "Post-mask_full",
}
MODES = ["method B"]        # GOLDEN single mode (was LCA/Coverage/EM)
MIN_READS = 3              # GOLDEN per-taxon read floor (was 1)


def collect():
    rows = []
    for rg in REGIMES:
        for md in MODES:
            for s in SAMPLES:
                kp = f"{ROOT}/{REGIME_DIRS[rg]}/{md}/{s}/{s}.kreport"
                if not os.path.exists(kp):
                    continue
                kdf = parse_kreport(kp)
                # species level
                c = kreport_to_species_counts(kdf)
                ev_sp = evaluate_sample(c, GROUND_TRUTH, min_reads=MIN_READS)
                # genus level
                ev_g  = evaluate_at_rank(kdf, GROUND_TRUTH_GENUS,
                                         rank_code='G', min_reads=MIN_READS)
                rows.append(dict(
                    regime=rg, mode=md, sample=s,
                    # species
                    TP_species=ev_sp["TP"],       FP_species=ev_sp["FP"],
                    FN_species=ev_sp["FN"],
                    Precision_species=ev_sp["precision"],
                    Recall_species=ev_sp["recall"],
                    F1_species=ev_sp["f1"],
                    # genus
                    TP_genus=ev_g["TP"],          FP_genus=ev_g["FP"],
                    FN_genus=ev_g["FN"],
                    Precision_genus=ev_g["precision"],
                    Recall_genus=ev_g["recall"],
                    F1_genus=ev_g["f1"],
                ))
    return pd.DataFrame(rows)


def render(df):
    # 4 rows: F1_species, FP_species, F1_genus, FP_genus
    # 1 col: method B (GOLDEN single mode)
    fig, axes = plt.subplots(4, len(MODES), figsize=(5.0 * len(MODES), 14.5),
                             sharex=True, squeeze=False,
                             gridspec_kw=dict(hspace=0.28, wspace=0.20,
                                              left=0.16, right=0.95,
                                              top=0.92, bottom=0.07))

    # Grayscale-safe per-sample styling.
    # MagNA (kit A) = darker greys      Qiagen (kit B) = lighter greys
    # Within each kit, the two replicates differ by hatch pattern,
    # so the figure remains strictly black-and-white but every bar is
    # visually distinguishable.
    # Strictly black-and-white styles.
    # All bar outlines, hatch lines and legend edges are BLACK.
    # MagNA = grey fill;  Qiagen = white fill.
    # Within each kit, replicate _1 = solid, replicate _2 = diagonal hatch.
    SAMPLE_STYLE = {
        "MagNA_1":  dict(facecolor="#9a9a9a", hatch=""),     # solid grey
        "MagNA_2":  dict(facecolor="#9a9a9a", hatch="////"), # grey  + black diag hatch
        "Qiagen_1": dict(facecolor="#ffffff", hatch=""),     # solid white
        "Qiagen_2": dict(facecolor="#ffffff", hatch="////"), # white + black diag hatch
    }

    plt.rcParams["hatch.linewidth"] = 0.9

    n_samples = len(SAMPLES)
    group_w   = 0.86          # total width of the 4-bar group at each regime
    bar_w     = group_w / n_samples
    x_pos = np.arange(len(REGIMES))
    # offsets so the 4 bars are centred on x_pos[i]
    offs = (np.arange(n_samples) - (n_samples - 1) / 2) * bar_w

    # rows: (metric column name, y-label, is_f1)
    ROW_SPEC = [
        ("F1_species", "F1", True),
        ("FP_species", "FP", False),
        ("F1_genus",   "F1", True),
        ("FP_genus",   "FP", False),
    ]
    fp_max = max(df[["FP_species", "FP_genus"]].max().max(), 1)

    for col, md in enumerate(MODES):
        for row, (metric, ylab, is_f1) in enumerate(ROW_SPEC):
            ax = axes[row, col]
            sub = df[df["mode"] == md]
            for j, s in enumerate(SAMPLES):
                style = SAMPLE_STYLE[s]
                vals = [
                    sub[(sub["regime"] == rg) & (sub["sample"] == s)][metric].mean()
                    for rg in REGIMES
                ]
                ax.bar(x_pos + offs[j], vals, width=bar_w * 0.92,
                       facecolor=style["facecolor"],
                       hatch=style["hatch"],
                       edgecolor="black", linewidth=0.9, zorder=2,
                       label=s if (row == 0 and col == 0) else None)
            ax.set_xticks(x_pos)
            ax.set_xticklabels(["No mask", "Pre-mask", "Post-mask"], fontsize=14)
            if row == 0:
                ax.set_title(md, fontsize=18, pad=10, fontweight="bold")
            if is_f1:
                ax.set_ylim(0.5, 1.05)
                ax.set_yticks(np.arange(0.5, 1.05, 0.1))
            else:
                ax.set_ylim(0, fp_max * 1.18)
            if col == 0:
                ax.set_ylabel(ylab, fontsize=16)
            ax.tick_params(axis="y", labelsize=13)
            ax.grid(axis="y", linestyle=":", linewidth=0.6,
                    color="black", alpha=0.35, zorder=0)
            # all four spines BLACK and visible — explicit
            for sp_name, sp in ax.spines.items():
                if sp_name in ("top", "right"):
                    sp.set_visible(False)
                else:
                    sp.set_visible(True)
                    sp.set_color("black")
                    sp.set_linewidth(1.0)
            # tick lines black
            ax.tick_params(axis="both", colors="black", width=1.0)

    # Bottom-row axis label
    for ax in axes[-1]:
        ax.set_xlabel("Masking regime", fontsize=15)

    # Big group labels on the left side: "species" for rows 0-1, "genus" for rows 2-3
    # Place rotated text aligned to the y-axis area, far enough left to not collide
    fig.text(0.018, 0.74, "Species level", rotation=90,
             ha="center", va="center", fontsize=20, fontweight="bold")
    fig.text(0.018, 0.28, "Genus level",   rotation=90,
             ha="center", va="center", fontsize=20, fontweight="bold")

    # Thin horizontal separator between species block (rows 0-1) and genus block (rows 2-3)
    sep_y = (axes[1, 0].get_position().y0 + axes[2, 0].get_position().y1) / 2
    fig.add_artist(plt.Line2D([0.05, 0.98], [sep_y, sep_y],
                              color="black", linewidth=0.8,
                              transform=fig.transFigure))

    # Shared legend below — 4 sample bars (all black outline + black hatch)
    from matplotlib.patches import Patch
    handles = [
        Patch(facecolor=SAMPLE_STYLE[s]["facecolor"],
              hatch=SAMPLE_STYLE[s]["hatch"],
              edgecolor="black", linewidth=0.7, label=s)
        for s in SAMPLES
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4,
               fontsize=14, frameon=False, bbox_to_anchor=(0.5, 0.005))

    return fig


def main():
    print("collecting per-sample metrics…")
    df = collect()
    print(df.to_string(index=False))

    df.to_csv(f"{NX}/result_260605/SupFig3/tables/masking_regimes_per_sample.csv", index=False)
    summary = (df.groupby(["regime", "mode"])
                 .agg(F1_species_mean=("F1_species", "mean"),
                      F1_species_sd  =("F1_species", "std"),
                      FP_species_mean=("FP_species", "mean"),
                      FP_species_sd  =("FP_species", "std"),
                      F1_genus_mean  =("F1_genus",   "mean"),
                      F1_genus_sd    =("F1_genus",   "std"),
                      FP_genus_mean  =("FP_genus",   "mean"),
                      FP_genus_sd    =("FP_genus",   "std"))
                 .round(3).reset_index())
    # order
    summary["regime"] = pd.Categorical(summary["regime"], categories=REGIMES, ordered=True)
    summary["mode"]   = pd.Categorical(summary["mode"],   categories=MODES,   ordered=True)
    summary = summary.sort_values(["regime", "mode"]).reset_index(drop=True)
    summary.to_csv(f"{NX}/result_260605/SupFig3/tables/masking_regimes_summary.csv", index=False)
    print("\n=== summary (mean +/- SD across 4 samples) ===")
    print(summary.to_string(index=False))

    print("\nrendering figure…")
    fig = render(df)
    save_figure(fig, "FigSx_masking_regimes_per_sample", close=True)


if __name__ == "__main__":
    main()
