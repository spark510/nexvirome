#!/usr/bin/env python3
"""
Fig4 D — TEST variant: per-sample, SPECIES-level non-phage composition.

Instead of D's four group-mean bars, this draws EVERY sample as its own stacked
bar, at SPECIES level (not genus), with the four clinical groups separated only by
a small horizontal gap (and a group label above each block):

    DNA·Control | DNA·Asthma   ||   RNA·Control | RNA·Asthma

Data lineage is identical to the production panel D (load_inputs() from
fig4_build_tables: Fig3 NexVirome species detections + asthma group manifest +
per-sample TPM weight; phage taxids removed via the pipeline's own _phage_set).
The only differences are: species (not genus) and per-sample (not group-mean).

Output: result_260603/fig4/test/Fig4D_persample_species_test.{png,eps}
Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260603/fig4/code/render_fig4d_persample_test.py
"""
from __future__ import annotations
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import fig4_build_tables as B   # reuse the production loaders / phage set / DB

mpl.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "font.size": 13, "axes.titlesize": 15, "axes.labelsize": 14,
    "xtick.labelsize": 8, "ytick.labelsize": 12, "legend.fontsize": 10,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 1.1,
})

NX = "/home/share/programs/nexvirome"
OUT = f"{NX}/result_260603/fig4/test"
TOP_N = 14   # most-abundant non-phage species shown individually; rest -> Other

# four clinical groups in display order, with a label and the gap before each
GROUPS = [("dna", "Control", "DNA · Control"),
          ("dna", "Asthma",  "DNA · Asthma"),
          ("rna", "Control", "RNA · Control"),
          ("rna", "Asthma",  "RNA · Asthma")]

# same soft tonal register as production D, its own palette order
PALETTE = ["#3A6EA5", "#C1666B", "#6A994E", "#BC9C22", "#8E7DBE",
           "#48A9A6", "#D08C60", "#4D7C8A", "#A24936", "#7A9E7E",
           "#9C89B8", "#5C5470", "#E8A87C", "#41658A"]
OTHER = "#B0B0B0"
WHITE_TEXT = {"#3A6EA5", "#4D7C8A", "#A24936", "#5C5470", "#9C89B8",
              "#8E7DBE", "#6A994E", "#48A9A6", "#41658A", "#C1666B"}


def build_persample_species():
    """Return (tpm_pivot, sample_order_per_group). tpm_pivot: index=species (incl.
    'Other'), columns=sample; each column sums to 1 over that sample's non-phage
    detections. Samples with no non-phage detection get an all-zero column."""
    nv, manifest = B.load_inputs()
    phage = B._phage_set(B.DB)
    nv = nv[~nv["taxid"].isin(phage)].copy()      # non-phage only

    # per-sample TPM renormalised over the sample's NON-PHAGE detections
    samp_tot = nv.groupby(["cohort", "sample"])["tpm"].transform("sum")
    nv["w"] = np.where(samp_tot > 0, nv["tpm"] / samp_tot, 0.0)

    # global top-N species by summed weight; everything else -> 'Other'
    top = (nv.groupby("name")["w"].sum().sort_values(ascending=False)
           .head(TOP_N).index.tolist())
    nv["label"] = np.where(nv["name"].isin(top), nv["name"], "Other")

    # sample order: per group, sort by that group's dominant detected sample first
    manifest = manifest.copy()
    order_by_group = {}
    for cohort, grp, _ in GROUPS:
        sams = manifest[(manifest.cohort == cohort) & (manifest.group == grp)]["sample"].tolist()
        # order samples within a group by their total non-phage weight (rich first)
        wt = nv[nv["sample"].isin(sams)].groupby("sample")["w"].sum()
        sams = sorted(sams, key=lambda s: -float(wt.get(s, 0.0)))
        order_by_group[(cohort, grp)] = sams

    piv = (nv.groupby(["sample", "label"])["w"].sum().unstack(fill_value=0.0))
    return piv, order_by_group, top


def main():
    piv, order_by_group, top = build_persample_species()
    species_order = top + ["Other"]
    colors = {sp: PALETTE[i % len(PALETTE)] for i, sp in enumerate(top)}
    colors["Other"] = OTHER

    # flat x layout with a gap between groups
    bar_w = 0.9
    gap = 1.6
    xs, xticklabels, group_spans = [], [], []
    x = 0.0
    for gi, (cohort, grp, glab) in enumerate(GROUPS):
        sams = order_by_group[(cohort, grp)]
        start = x
        for s in sams:
            xs.append((s, x))
            xticklabels.append(s)
            x += 1.0
        group_spans.append((glab, start, x - 1.0, len(sams)))
        x += gap   # blank gap before the next group

    fig, ax = plt.subplots(figsize=(max(13, 0.16 * len(xs) + 4), 6.4))
    for s, xpos in xs:
        col_data = piv.loc[s] if s in piv.index else None
        bottom = 0.0
        for sp in species_order:
            v = float(col_data[sp]) if (col_data is not None and sp in col_data) else 0.0
            if v <= 0:
                continue
            ax.bar(xpos, v, bar_w, bottom=bottom, color=colors[sp],
                   edgecolor="white", linewidth=0.25)
            bottom += v

    # group labels + faint separators
    for glab, x0, x1, n in group_spans:
        cx = (x0 + x1) / 2
        ax.text(cx, 1.02, f"{glab}  (n={n})", ha="center", va="bottom",
                fontsize=12, weight="bold")
        ax.axvspan(x0 - 0.6, x1 + 0.6, color="0.96", zorder=0)

    ax.set_xticks([xp for _, xp in xs])
    ax.set_xticklabels(xticklabels, rotation=90, fontsize=6)
    ax.set_xlim(-1.0, xs[-1][1] + 1.0)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel("Per-sample TPM rel. abundance\n(non-phage, species)")
    # title placed well above the group labels so they don't overlap
    fig.suptitle("Fig 4D test — per-sample, species-level non-phage composition",
                 fontsize=14, y=0.99)

    handles = [Patch(facecolor=colors[sp], edgecolor="white", label=sp)
               for sp in species_order]
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.005, 0.5),
              frameon=False, fontsize=9, handlelength=1.1, labelspacing=0.35,
              title="Non-phage species", title_fontsize=10)

    fig.subplots_adjust(left=0.07, right=0.80, top=0.86, bottom=0.20)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/Fig4D_persample_species_test.{ext}", dpi=300,
                    bbox_inches="tight")
    print(f"-> {OUT}/Fig4D_persample_species_test.png , .eps")
    print(f"   samples: {len(xs)} ; top species shown: {len(top)}")


if __name__ == "__main__":
    main()
