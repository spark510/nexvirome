#!/usr/bin/env python3
"""
Supplementary Figure for Table S2 — read-level parameter sensitivity (no
downstream gate). Three columns (identity / query coverage / aligned length).

Row 1  KIT mock (ground truth): F1 (line, left axis) and FP (line, right axis).
Row 2  Real cohorts: detected-species richness vs the parameter — DNA and RNA
       (means over 5 representative samples each).

Message is SENSITIVITY, not optimality: detection depends strongly on each
read-level filter, so the threshold must be chosen deliberately. The production
default is marked with a dashed vertical line in every panel (shown for reference,
NOT claimed to be the F1 optimum — without any downstream gate the KIT F1 keeps
rising past the default).

Input: result_260605/SuppTab2/tables/param_sensitivity_summary.csv
Output: result_260605/SuppTab2/FigS_param_sensitivity.{png,eps}
Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260605/SuppTab2/code/render_param_sensitivity.py
"""
from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

mpl.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 13, "axes.titlesize": 15, "axes.labelsize": 13,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 11,
    "axes.spines.top": False, "axes.linewidth": 1.0,
})

NX = "/home/share/programs/nexvirome"
TBL = f"{NX}/result_260605/SuppTab2/tables"
OUT = f"{NX}/result_260605/SuppTab2"

PARAMS = [("identity", "Min identity", 0.85),
          ("qcov", "Min query coverage", 0.50),
          ("length", "Min aligned length (bp)", 60)]
PREC_COL = "#1D4E89"      # KIT precision (left)
REC_COL = "#2E8B57"       # KIT recall (right)
DNA_COL = "#1D4E89"
RNA_COL = "#E07A5F"
DEF_LS = dict(color="grey", ls="--", lw=1.3, zorder=0)


def _kit_precision_recall():
    """Per-(param,value) KIT precision & recall, computed from the per-sample
    TP/FP (the summary CSV has no precision column). precision = TP/(TP+FP)
    averaged over the 4 KIT samples; recall is averaged likewise."""
    p = pd.read_csv(f"{TBL}/param_sensitivity_per_sample.csv")
    k = p[p["cohort"] == "KIT"].copy()
    denom = k["TP"] + k["FP"]
    k["precision"] = (k["TP"] / denom).where(denom > 0, 0.0)
    g = (k.groupby(["param", "value"])
           .agg(precision_mean=("precision", "mean"),
                recall_mean=("recall", "mean"))
           .reset_index())
    return g


def main():
    s = pd.read_csv(f"{TBL}/param_sensitivity_summary.csv")
    pr = _kit_precision_recall()
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.0),
                             gridspec_kw=dict(hspace=0.42, wspace=0.5,
                                              left=0.07, right=0.93,
                                              top=0.92, bottom=0.10))

    for j, (param, xlabel, default) in enumerate(PARAMS):
        # ---- Row 1: KIT precision (left) + recall (right) ----
        ax = axes[0, j]
        k = pr[pr["param"] == param].sort_values("value")
        x = k["value"].values
        ax.plot(x, k["precision_mean"], "-o", color=PREC_COL, lw=2.0, ms=5,
                label="precision")
        ax.set_ylim(0, 1.05); ax.set_ylabel("KIT precision", color=PREC_COL)
        ax.tick_params(axis="y", labelcolor=PREC_COL)
        ax2 = ax.twinx(); ax2.spines["top"].set_visible(False)
        ax2.plot(x, k["recall_mean"], "-s", color=REC_COL, lw=2.0, ms=5,
                 label="recall")
        ax2.set_ylim(0, 1.05)
        ax2.set_ylabel("KIT recall", color=REC_COL)
        ax2.tick_params(axis="y", labelcolor=REC_COL)
        ax.axvline(default, **DEF_LS)
        ax.set_title(xlabel)
        ax.set_xlabel(xlabel)

        # ---- Row 2: DNA / RNA detected-species richness ----
        ax = axes[1, j]
        for coh, col in (("DNA", DNA_COL), ("RNA", RNA_COL)):
            d = s[(s.cohort == coh) & (s.param == param)].sort_values("value")
            ax.plot(d["value"], d["n_species_mean"], "-o", color=col, lw=2.0,
                    ms=5, label=coh)
        ax.axvline(default, **DEF_LS)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Detected species\n(real data, mean)")
        if j == 0:
            ax.legend(frameon=False, loc="upper right")

    # row labels
    fig.text(0.012, 0.71, "KIT mock\n(ground truth)", rotation=90,
             ha="center", va="center", fontsize=13, fontweight="bold")
    fig.text(0.012, 0.29, "Real cohorts\n(richness)", rotation=90,
             ha="center", va="center", fontsize=13, fontweight="bold")
    # KIT legend (precision/recall) once, top-left panel
    h1, l1 = axes[0, 0].get_legend_handles_labels()
    axes[0, 0].legend(h1 + [plt.Line2D([], [], color=REC_COL, marker="s", lw=2)],
                      l1 + ["recall"], frameon=False, loc="lower right")
    fig.suptitle("Detection sensitivity to read-level filters "
                 "(no breadth / read-floor / abundance gate; dashed = production default)",
                 fontsize=13, y=0.985)

    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/FigS_param_sensitivity.{ext}", dpi=300, bbox_inches="tight")
    print(f"-> {OUT}/FigS_param_sensitivity.png , .eps")


if __name__ == "__main__":
    main()
