#!/usr/bin/env python3
"""
SuppTab2 — combined parameter-stability heatmaps.

Three 2D heatmaps (one per parameter pair: identity×qcov, identity×length,
qcov×length), colour = KIT-mock precision at that grid cell (third parameter at
its production default). A broad flat high-precision plateau = a stable operating
region; the production default (identity 0.85, qcov 0.5, length 60) is marked.

Reads result_260605/SuppTab2/tables/param_grid_KIT_summary.csv (from
build_param_grid.py).

Output: result_260605/SuppTab2/test/FigS_param_grid_heatmap.{png,eps}
Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260605/SuppTab2/code/render_param_grid_heatmap.py
"""
from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

mpl.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
    "xtick.labelsize": 10, "ytick.labelsize": 10,
})

NX = "/home/share/programs/nexvirome"
TBL = f"{NX}/result_260605/SuppTab2/tables"
OUT = f"{NX}/result_260605/SuppTab2/test"

DEF = {"identity": 0.85, "qcov": 0.50, "length": 60}
LABEL = {"identity": "Min identity", "qcov": "Min query coverage",
         "length": "Min aligned length (bp)"}
PAIRS = [("identity", "qcov"), ("identity", "length"), ("qcov", "length")]
METRIC = "precision_mean"     # colour = KIT precision


def main():
    s = pd.read_csv(f"{TBL}/param_grid_KIT_summary.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6),
                             gridspec_kw=dict(wspace=0.42, left=0.06, right=0.97,
                                              top=0.86, bottom=0.16))
    vmin, vmax = 0.0, float(s[METRIC].max())
    im = None
    for ax, (xp, yp) in zip(axes, PAIRS):
        sub = s[s["pair"] == f"{xp}x{yp}"]
        piv = sub.pivot(index="y_value", columns="x_value", values=METRIC)
        piv = piv.sort_index(ascending=True).sort_index(axis=1)
        im = ax.imshow(piv.values, origin="lower", aspect="auto", cmap="viridis",
                       vmin=vmin, vmax=vmax)
        ax.set_xticks(range(len(piv.columns)))
        ax.set_xticklabels([f"{c:g}" for c in piv.columns], rotation=0)
        ax.set_yticks(range(len(piv.index)))
        ax.set_yticklabels([f"{r:g}" for r in piv.index])
        ax.set_xlabel(LABEL[xp]); ax.set_ylabel(LABEL[yp])
        ax.set_title(f"{LABEL[xp]} × {LABEL[yp]}\n(third param at default)",
                     fontsize=11)
        # annotate each cell with the value
        for iy in range(len(piv.index)):
            for ix in range(len(piv.columns)):
                v = piv.values[iy, ix]
                if not np.isnan(v):
                    ax.text(ix, iy, f"{v:.2f}", ha="center", va="center",
                            fontsize=14,
                            color="white" if v < 0.55 * vmax else "black")
        # mark the production default cell
        if DEF[xp] in list(piv.columns) and DEF[yp] in list(piv.index):
            cx = list(piv.columns).index(DEF[xp])
            cy = list(piv.index).index(DEF[yp])
            ax.add_patch(plt.Rectangle((cx - 0.5, cy - 0.5), 1, 1, fill=False,
                                       edgecolor="red", lw=2.2))

    cbar = fig.colorbar(im, ax=axes, fraction=0.025, pad=0.015)
    cbar.set_label("KIT precision  (TP / (TP+FP))")
    fig.suptitle("Combined read-level filter stability — KIT-mock precision over "
                 "parameter-pair grids (red box = production default)",
                 fontsize=12.5, y=0.99)

    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/FigS_param_grid_heatmap.{ext}", dpi=300,
                    bbox_inches="tight")
    print(f"-> {OUT}/FigS_param_grid_heatmap.png , .eps")


if __name__ == "__main__":
    main()
