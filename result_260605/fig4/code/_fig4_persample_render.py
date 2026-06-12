#!/usr/bin/env python3
"""
Single render routine for the Fig4 per-sample 4-view family. A driver script in
each view folder calls `render_view(level, metric, outdir, tag)`.

Each figure has two stacked-bar rows, per sample, grouped (DNA·Control /
DNA·Asthma / RNA·Control / RNA·Asthma) with gaps + labels:
  (A) genome type (Baltimore class) — always read-fraction per sample
  (B) integrated composition        — phage->host genus + non-phage->taxon,
                                       at `level`, weighted by `metric`
"""
from __future__ import annotations
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import _fig4_persample_common as C
import fig4_build_tables as B   # for label_detections (genome-type / molecule)
import render_fig4_dna_rna as R  # CLASS_ORDER / CLASS_COLOR

mpl.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13,
    "xtick.labelsize": 6, "ytick.labelsize": 11, "legend.fontsize": 13.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.0,
})

OTHER_GT = "#B0B0B0"


def _genometype_matrix(df):
    """Per-sample read-fraction over Baltimore class (phage + non-phage together)."""
    d = B.label_detections(df.rename(columns={}), B.DB)   # adds `molecule`
    d["w"] = d["reads"].astype(float)
    st = d.groupby("sample")["w"].transform("sum")
    d["w"] = np.where(st > 0, d["w"] / st, 0.0)
    piv = d.groupby(["molecule", "sample"])["w"].sum().unstack(fill_value=0.0)
    order = [c for c in R.CLASS_ORDER if c in piv.index]
    return piv, order


def _draw(ax, piv, order, colors, xs, spans, ylabel, ymax=None,
          show_xticklabels=False, logy=False, letter=None):
    bar_w = 0.9
    FLOOR = 0.8   # log-axis baseline (reads are integers >= 1)
    for s, xp in xs:
        col = piv[s] if s in piv.columns else None
        bottom = FLOOR if logy else 0.0
        for lab in order:
            v = float(col[lab]) if (col is not None and lab in col.index) else 0.0
            if v <= 0:
                continue
            if logy:
                # stack on a log axis: segment spans [bottom, bottom+v]; draw it as
                # a bar from `bottom` of height v (matplotlib handles the log mapping)
                ax.bar(xp, v, bar_w, bottom=bottom, color=colors[lab],
                       edgecolor="white", linewidth=0.15)
                bottom += v
            else:
                ax.bar(xp, v, bar_w, bottom=bottom, color=colors[lab],
                       edgecolor="white", linewidth=0.2)
                bottom += v
    sample_tot = [float(piv[s].sum()) for s, _ in xs if s in piv.columns] + [1.0]
    if logy:
        ax.set_yscale("log")
        top = max(sample_tot)
        ax.set_ylim(FLOOR, top * 1.6)
        lab_y = top * 1.7
    else:
        top = ymax if ymax is not None else max(sample_tot)
        ax.set_ylim(0, top)
        lab_y = top * 1.02
    for glab, x0, x1, n in spans:
        ax.axvspan(x0 - 0.6, x1 + 0.6, color="0.96", zorder=0)
        ax.text((x0 + x1) / 2, lab_y, glab, ha="center",
                va="bottom", fontsize=10.5, weight="bold")
    ax.set_xticks([xp for _, xp in xs])
    ax.set_xticklabels([s for s, _ in xs] if show_xticklabels else [],
                       rotation=90, fontsize=5)
    ax.set_xlim(-1.0, xs[-1][1] + 1.0)
    ax.set_ylabel(ylabel)
    if letter:
        ax.text(-0.065, 1.06, letter, transform=ax.transAxes, fontsize=20,
                weight="bold", va="bottom", ha="right")


def _legend(ax, order, colors, title=None):
    # title intentionally dropped from all legends (per request)
    h = [Patch(facecolor=colors[l], edgecolor="white", label=C.short_name(l))
         for l in order]
    ax.legend(handles=h, loc="center left", bbox_to_anchor=(1.004, 0.5),
              frameon=False, fontsize=12.5, handlelength=1.1, labelspacing=0.28)


def render_view(level, metric, outdir, tag):
    """level: 'species'|'genus' ; metric: 'relabund'|'reads'."""
    df = C.load_labeled(level)
    man = C.manifest()
    xs, spans = C.sample_layout(df, man)

    # row 1 — genome type (always read fraction)
    pivA, ordA = _genometype_matrix(df)
    colA = {c: R.CLASS_COLOR.get(c, OTHER_GT) for c in ordA}

    # row 2 — phage vs non-phage split (read fraction per sample)
    pivP, ordP = C.phage_split_matrix(df, metric="reads")
    colP = C.PHAGE_KIND_COLOR

    # row 3 — PHAGE only (rolled to host genus), renormalised within phage subset.
    # top_n=9 (+Other) keeps the larger legend font within the panel height.
    pivB, ordB = C.build_matrix(df, metric, top_n=9, subset="phage")
    colB = C.colour_map(ordB)

    # row 4 — NON-PHAGE only, renormalised within the non-phage subset
    pivN, ordN = C.build_matrix(df, metric, top_n=9, subset="non-phage")
    colN = C.colour_map(ordN, palette=C.D_PALETTE)   # distinct palette from C

    metric_lab = ("relative abundance" if metric == "relabund" else "read count")
    ymaxB = 1.0 if metric == "relabund" else None
    logyB = (metric == "reads")   # read counts span orders of magnitude -> log y
    suffix = " (log scale)" if logyB else ""
    ylabB = (f"phage host-genus\n{metric_lab}{suffix}")
    ylabN = (f"non-phage {level}\n{metric_lab}{suffix}")

    width = max(15, 0.14 * len(xs) + 5)
    # A,B made shorter again (~4/9 of C,D); tighter vertical spacing
    # overall height scaled to 3/4 (13.0 -> 9.75)
    fig = plt.figure(figsize=(width, 9.75))
    gs = GridSpec(4, 1, figure=fig, hspace=0.20, height_ratios=[1.35, 1.35, 3.0, 3.0])

    axA = fig.add_subplot(gs[0])
    _draw(axA, pivA, ordA, colA, xs, spans,
          "genome type\n(read fraction)", ymax=1.0, show_xticklabels=False, letter="A")
    _legend(axA, ordA, colA, "genome type")

    axB = fig.add_subplot(gs[1])
    _draw(axB, pivP, ordP, colP, xs, spans,
          "phage / non-phage\n(read fraction)", ymax=1.0, show_xticklabels=False, letter="B")
    _legend(axB, ordP, colP, "kind")

    axC = fig.add_subplot(gs[2])
    _draw(axC, pivB, ordB, colB, xs, spans, ylabB, ymax=ymaxB,
          show_xticklabels=False, logy=logyB, letter="C")
    _legend(axC, ordB, colB, "phage host genus")

    axD = fig.add_subplot(gs[3])
    _draw(axD, pivN, ordN, colN, xs, spans, ylabN, ymax=ymaxB,
          show_xticklabels=True, logy=logyB, letter="D")
    _legend(axD, ordN, colN, f"non-phage ({level})")

    # no suptitle — the per-panel group labels (A/B/C/D) would collide with it
    fig.subplots_adjust(left=0.07, right=0.82, top=0.97, bottom=0.08)
    os.makedirs(outdir, exist_ok=True)
    base = f"{outdir}/Fig4_persample_{tag}"
    for ext in ("png", "eps"):
        fig.savefig(f"{base}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"-> {base}.png , .eps   ({len(xs)} samples | A {len(ordA)} classes | "
          f"B {len(ordB)} taxa)")
