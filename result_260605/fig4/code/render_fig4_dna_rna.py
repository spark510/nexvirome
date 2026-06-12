#!/usr/bin/env python3
"""
Fig4 — how the library protocol (DNA vs RNA) shapes the detected virome
(asthma cohort, NexVirome method B = Fig2/Fig3 pipeline). Each panel is a pair of
cohort-mean stacked bars (DNA vs RNA); panels carry NO titles (described in the
caption), per request.

  (A)  Genome-type composition           — Baltimore class, mean read fraction
  (B)  Phage host-genus composition      — phage rolled to bacterial host genus
  (C)  Non-phage composition, species    — eukaryotic-virus species
  (C') Non-phage composition, genus      — same, summarised to genus

Inputs (result_260603/fig4/tables/, built by fig4_build_tables.py):
  fig4_genometype_cohort.csv      cohort, molecule, read_frac, ...
  fig4_hostgenus_cohort.csv       cohort, host_genus, frac
  fig4_nonphage_species_cohort.csv cohort, label, frac
  fig4_nonphage_genus_cohort.csv   cohort, label, frac
Output: result_260603/fig4/Fig4_dna_rna.{png,eps}
Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260603/fig4/code/render_fig4_dna_rna.py
"""
from __future__ import annotations
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

mpl.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 18,
    "axes.titlesize": 20,
    "axes.labelsize": 18,
    "xtick.labelsize": 15,
    "ytick.labelsize": 14,
    "legend.fontsize": 19.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 1.0,
})

NX = "/home/share/programs/nexvirome"
import sys as _sys
_sys.path.insert(0, f"{NX}/result_260605")
from golden_rule import keep_samples
TBL = f"{NX}/result_260605/fig4/tables"
OUT = f"{NX}/result_260605/fig4"

COHORTS = ["dna", "rna"]
XLAB = ["DNA\nlibrary", "RNA\nlibrary"]

# ---- palettes ----
# (A) Baltimore class — same as paper/figures/Fig4
CLASS_ORDER = ["dsDNA", "ssDNA", "ssRNA", "dsRNA", "RT", "Other"]
CLASS_COLOR = {"dsDNA": "#1D4E89", "ssDNA": "#00B2CA", "ssRNA": "#FF595E",
               "dsRNA": "#7DCFB6", "RT": "#FBD1A2", "Other": "#B0B0B0"}

# (B) host genus
HOST_ORDER = ["Pseudomonas", "Streptococcus", "Bacteroides", "Fusobacterium",
              "Rothia", "Haemophilus", "Staphylococcus", "Actinomyces",
              "Phlebotominae", "Homo", "Other host"]
HOST_COLOR = {
    "Pseudomonas": "#1D4E89", "Streptococcus": "#E07A5F", "Bacteroides": "#3CAEA3",
    "Fusobacterium": "#F6BD60", "Rothia": "#9B5DE5", "Haemophilus": "#00BBF9",
    "Staphylococcus": "#F15BB5", "Actinomyces": "#8AC926", "Phlebotominae": "#C9ADA7",
    "Homo": "#6D6875", "Other host": "#B0B0B0"}

# Panels C and D share a tonal register (medium-saturation, soft) but use
# DIFFERENT palettes so the two composition panels read as visually distinct.
# (C) host-genus extras fall back to this palette; (D) non-phage uses its own.
_HOST_PALETTE = ["#1D4E89", "#E07A5F", "#3CAEA3", "#F6BD60", "#9B5DE5",
                 "#00BBF9", "#F15BB5", "#8AC926", "#FF924C", "#577590",
                 "#C9ADA7", "#6D6875"]
# (D) non-phage genus — same softness, shifted hue families (teal-blue / rust /
# olive / plum / slate) so D never mirrors C's colour order.
_NONPHAGE_PALETTE = ["#3A6EA5", "#C1666B", "#6A994E", "#BC9C22", "#8E7DBE",
                     "#48A9A6", "#D08C60", "#4D7C8A", "#A24936", "#7A9E7E",
                     "#9C89B8", "#5C5470"]
WHITE_TEXT = {"#1D4E89", "#577590", "#6D6875", "#9B5DE5", "#3CAEA3",
              "#3A6EA5", "#4D7C8A", "#A24936", "#5C5470", "#9C89B8",
              "#8E7DBE", "#6A994E", "#48A9A6"}


def _stacked(ax, mat, order, colors, xlabel):
    """mat: DataFrame index=label, cols=COHORTS (fractions). Draw two HORIZONTAL
    stacked bars (DNA top, RNA bottom), label segments >=6%. No title."""
    y = np.arange(len(COHORTS))
    lefts = {c: 0.0 for c in COHORTS}
    for lab in order:
        if lab not in mat.index:
            continue
        vals = [float(mat.at[lab, c]) for c in COHORTS]
        col = colors[lab]
        ax.barh(y, vals, 0.62, left=[lefts[c] for c in COHORTS],
                color=col, edgecolor="white", linewidth=0.6, label=lab)
        for yi, c in enumerate(COHORTS):
            if vals[yi] >= 0.06:
                ax.text(lefts[c] + vals[yi] / 2, yi, f"{vals[yi]*100:.0f}%",
                        ha="center", va="center", fontsize=11,
                        color="white" if col in WHITE_TEXT else "black")
            lefts[c] += vals[yi]
    ax.set_yticks(y); ax.set_yticklabels(XLAB)
    ax.set_ylim(-0.6, len(COHORTS) - 0.4)
    ax.invert_yaxis()                       # DNA on top, RNA below
    ax.set_xlim(0, 1.0)
    ax.set_xlabel(xlabel)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              fontsize=18, handlelength=1.1, labelspacing=0.32)


def _pivot(csv, label_col, order=None):
    """Load a cohort/label/frac CSV -> DataFrame index=label, cols=COHORTS. If
    `order` is None, order labels by descending max-cohort frac with 'Other'(*) last."""
    df = pd.read_csv(csv)
    m = df.pivot_table(index=label_col, columns="cohort", values="frac",
                       fill_value=0.0)
    for c in COHORTS:
        if c not in m.columns:
            m[c] = 0.0
    if order is None:
        rank = m.max(axis=1).sort_values(ascending=False)
        order = [l for l in rank.index if not str(l).startswith("Other")]
        order += [l for l in rank.index if str(l).startswith("Other")]
    return m, order


# ---- cohort x clinical-group (4 bars) support for panels C/D ----
GROUP_ROWS = [("dna", "Control"), ("dna", "Asthma"),
              ("rna", "Control"), ("rna", "Asthma")]
ASTHMA_TBL = f"{NX}/result_260605/fig4/tables"   # cohort×group tables now live in fig4 (fig4_asthma merged in)


def _nsamp(cohort, group):
    g = pd.read_csv(f"{NX}/paper/figures/Fig3/source_data/{cohort}_asthma_groups.csv")
    # result_260605: vir17 (DNA Asthma) is excluded -> displayed n must reflect that.
    g = g[g["sample"].isin(keep_samples(g["sample"]))]
    return int((g["group"] == group).sum())


def _pivot_grp(csv, label_col):
    """cohort,group,label,frac -> DataFrame index=label, cols=(cohort,group)."""
    df = pd.read_csv(csv)
    m = df.pivot_table(index=label_col, columns=["cohort", "group"],
                       values="frac", fill_value=0.0)
    rank = m.max(axis=1).sort_values(ascending=False)
    order = [l for l in rank.index if not str(l).startswith("Other")]
    order += [l for l in rank.index if str(l).startswith("Other")]
    return m, order


def _stacked4(ax, m, order, colors, xlabel):
    """4 horizontal stacked bars (DNA-Ctrl/DNA-Asth/RNA-Ctrl/RNA-Asth)."""
    y = np.arange(len(GROUP_ROWS))
    lefts = np.zeros(len(GROUP_ROWS))
    for lab in order:
        if lab not in m.index:
            continue
        vals = np.array([float(m.at[lab, (c, g)]) if (c, g) in m.columns else 0.0
                         for c, g in GROUP_ROWS])
        col = colors[lab]
        ax.barh(y, vals, 0.66, left=lefts, color=col, edgecolor="white",
                linewidth=0.6, label=lab)
        for yi, v in enumerate(vals):
            if v >= 0.06:
                ax.text(lefts[yi] + v / 2, yi, f"{v*100:.0f}%", ha="center",
                        va="center", fontsize=10,
                        color="white" if col in WHITE_TEXT else "black")
        lefts += vals
    ax.set_yticks(y)
    ax.set_yticklabels([f"{c.upper()} · {g}\n(n={_nsamp(c, g)})" for c, g in GROUP_ROWS],
                       fontsize=11)
    ax.set_ylim(-0.6, len(GROUP_ROWS) - 0.4); ax.invert_yaxis()
    ax.set_xlim(0, 1.0); ax.set_xlabel(xlabel)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              fontsize=18, handlelength=1.1, labelspacing=0.30)


# phage/non-phage colours chosen to stand clearly apart from the genome-type
# palette (no blues/reds): muted teal-green vs warm gold.
PHAGE_KIND_COLOR = {"phage": "#4C9F70", "non-phage": "#E6B800"}


def panel_genometype(ax):
    """Panel A — 4 horizontal bars in cohort order DNA, DNA, RNA, RNA:
       row0 DNA  genome-type (TPM)     row1 DNA  phage / non-phage (TPM)
       row2 RNA  genome-type (TPM)     row3 RNA  phage / non-phage (TPM)
    Two metrics share the panel; each row is its own composition (sums to 1)."""
    gt = pd.read_csv(f"{TBL}/fig4_genometype_cohort.csv")
    pk = pd.read_csv(f"{TBL}/fig4_phage_kind_cohort.csv")
    gm = gt.pivot_table(index="molecule", columns="cohort", values="tpm_frac",
                        fill_value=0.0)
    pm = pk.pivot_table(index="kind", columns="cohort", values="tpm_frac",
                        fill_value=0.0)

    # rows, top->bottom: DNA-genometype, DNA-phagekind, RNA-genometype, RNA-phagekind
    rows = [("dna", "genometype"), ("dna", "phagekind"),
            ("rna", "genometype"), ("rna", "phagekind")]
    ylabels = ["DNA  genome type", "DNA  phage / non-phage",
               "RNA  genome type", "RNA  phage / non-phage"]
    y = np.arange(len(rows))
    for yi, (coh, metric) in enumerate(rows):
        left = 0.0
        if metric == "genometype":
            order, colors, src = [c for c in CLASS_ORDER if c in gm.index], CLASS_COLOR, gm
        else:
            order, colors, src = ["phage", "non-phage"], PHAGE_KIND_COLOR, pm
        for lab in order:
            v = float(src.at[lab, coh]) if lab in src.index else 0.0
            if v <= 0:
                continue
            col = colors[lab]
            ax.barh(yi, v, 0.7, left=left, color=col, edgecolor="white", linewidth=0.6)
            if v >= 0.06:
                ax.text(left + v / 2, yi, f"{v*100:.0f}%", ha="center", va="center",
                        fontsize=10, color="white" if col in WHITE_TEXT
                        or col == PHAGE_KIND_COLOR["phage"] else "black")
            left += v
    ax.set_yticks(y); ax.set_yticklabels(ylabels, fontsize=13)
    ax.set_ylim(-0.6, len(rows) - 0.4); ax.invert_yaxis()
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Mean TPM relative abundance")
    # two legends: genome type + phage kind
    gt_handles = [Patch(facecolor=CLASS_COLOR[c], edgecolor="white", label=c)
                  for c in CLASS_ORDER if c in gm.index]
    pk_handles = [Patch(facecolor=PHAGE_KIND_COLOR[k], edgecolor="white", label=k)
                  for k in ("phage", "non-phage")]
    leg1 = ax.legend(handles=gt_handles, loc="center left",
                     bbox_to_anchor=(1.01, 0.72), frameon=False, fontsize=16.5,
                     title="genome type", title_fontsize=16.5, handlelength=1.1,
                     labelspacing=0.3)
    ax.add_artist(leg1)
    ax.legend(handles=pk_handles, loc="center left", bbox_to_anchor=(1.01, 0.25),
              frameon=False, fontsize=16.5, title="phage kind", title_fontsize=16.5,
              handlelength=1.1, labelspacing=0.3)


def panel_hostgenus(ax):
    """C — phage host-genus composition, cohort x clinical group (4 bars),
    from the fig4_asthma cohort-group table."""
    m, _ = _pivot_grp(f"{ASTHMA_TBL}/hostgenus_cohort_group.csv", "host_genus")
    extra = [h for h in m.index if h not in HOST_ORDER and h != "Other host"]
    order = [h for h in HOST_ORDER if h in m.index and h != "Other host"]
    order += sorted(extra)
    if "Other host" in m.index:
        order += ["Other host"]
    colors = dict(HOST_COLOR)
    # assign fallback colours to extra genera, skipping any hue already used by a
    # named host genus that is present (otherwise e.g. Aggregatibacter would reuse
    # Pseudomonas' navy and the panel shows the same colour twice).
    used = {colors[h] for h in order if h in colors}
    free = [c for c in _HOST_PALETTE if c not in used]
    for i, h in enumerate(extra):
        colors.setdefault(h, free[i] if i < len(free)
                          else _HOST_PALETTE[i % len(_HOST_PALETTE)])
    def _plab(h):
        return "Other host phage" if h == "Other host" else f"{h} phage"
    m2 = m.rename(index=_plab)
    order2 = [_plab(h) for h in order]
    colors2 = {_plab(h): colors[h] for h in order}
    _stacked4(ax, m2, order2, colors2, "Mean phage TPM rel. abundance (by host genus)")


def _nonphage_colors(order):
    cols, i = {}, 0
    for lab in order:
        if str(lab).startswith("Other"):
            cols[lab] = "#B0B0B0"
        else:
            cols[lab] = _NONPHAGE_PALETTE[i % len(_NONPHAGE_PALETTE)]; i += 1
    return cols


def panel_nonphage(ax, level, xlabel):
    """D — non-phage genus composition, cohort x clinical group (4 bars),
    from the fig4_asthma cohort-group table. (level kept for API compat; genus.)"""
    m, order = _pivot_grp(f"{ASTHMA_TBL}/nonphage_genus_cohort_group.csv", "label")
    cols = _nonphage_colors(order)
    _stacked4(ax, m, order, cols, xlabel)


# phage-fraction tables now in fig4/tables (built by phage_ratio_dna_rna.py)
RATIO_TBL = f"{NX}/result_260605/fig4/tables"
# (B) DNA vs RNA — chosen to NOT clash with panel A's Baltimore palette
# (A uses navy #1D4E89 / cyan / coral-red); B uses muted violet + sea-green.
COH_RATIO_COLOR = {"dna": "#7B6CA6", "rna": "#5BA081"}


def panel_phage_ratio(ax, metric="read"):
    """Per-sample phage fraction, DNA vs RNA (vertical box + points). Default
    metric='read' (the significant one, p~0.04); TPM is n.s. (length effect)."""
    d = pd.read_csv(f"{RATIO_TBL}/phage_ratio_per_sample_DNA_RNA.csv")
    t = pd.read_csv(f"{RATIO_TBL}/phage_ratio_test_DNA_RNA.csv")
    d = d[d["metric"] == metric]
    for xi, c in enumerate(["dna", "rna"]):
        vals = d[d["cohort"] == c]["phage_frac"].values
        col = COH_RATIO_COLOR[c]
        # filled box in the cohort colour, black outline / median / whiskers
        # at the default (thin) line width
        ax.boxplot([vals], positions=[xi], widths=0.55, patch_artist=True,
                   showfliers=False, medianprops=dict(color="black"),
                   whiskerprops=dict(color="black"),
                   capprops=dict(color="black"),
                   boxprops=dict(facecolor=col, edgecolor="black", alpha=0.9))
        jit = (np.random.RandomState(xi).rand(len(vals)) - 0.5) * 0.26
        ax.scatter(np.full(len(vals), xi) + jit, vals, s=22, color=col,
                   edgecolor="black", linewidth=0.5, alpha=0.95, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["DNA", "RNA"])
    ax.set_xlim(-0.6, 1.6); ax.set_ylim(0, 1.08)
    ax.set_ylabel("Phage fraction")
    row = t[t["metric"] == metric].iloc[0]
    p = float(row["p"])
    star = "n.s." if p >= 0.05 else ("*" if p >= 0.01 else "**")
    y = 1.0
    ax.plot([0, 0, 1, 1], [y, y + 0.02, y + 0.02, y], lw=1.2, color="black")
    ax.text(0.5, y + 0.025, f"{star} (p={p:.3f})", ha="center", va="bottom", fontsize=12)


def main():
    fig = plt.figure(figsize=(15, 13))
    # top row: A (overall, wide) | B (phage ratio, narrow); then C, D (4-bar
    # cohort x group each, taller) full width.
    outer = GridSpec(3, 1, figure=fig, hspace=0.55, height_ratios=[1.7, 1.35, 1.35])
    top = outer[0].subgridspec(1, 2, width_ratios=[3.2, 1.0], wspace=0.62)

    ax_a = fig.add_subplot(top[0])
    panel_genometype(ax_a)                          # A
    ax_b = fig.add_subplot(top[1])
    panel_phage_ratio(ax_b, metric="read")          # B
    ax_c = fig.add_subplot(outer[1])
    panel_hostgenus(ax_c)                           # C
    ax_d = fig.add_subplot(outer[2])
    panel_nonphage(ax_d, "genus",                   # D
                   "Mean TPM rel. abundance (non-phage, genus)")

    for ax, lab, dx in [(ax_a, "A", -0.10), (ax_b, "B", -0.32),
                        (ax_c, "C", -0.10), (ax_d, "D", -0.10)]:
        ax.text(dx, 1.06, lab, transform=ax.transAxes,
                fontsize=24, weight="bold", va="bottom")

    fig.subplots_adjust(left=0.12, right=0.82, top=0.96, bottom=0.06)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/Fig4_dna_rna.{ext}", dpi=300, bbox_inches="tight")
    print(f"-> {OUT}/Fig4_dna_rna.png , .eps")


if __name__ == "__main__":
    main()
