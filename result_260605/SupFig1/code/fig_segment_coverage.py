#!/usr/bin/env python3
"""
Supplementary figure: per-segment coverage of a multipartite virus (Reovirus 3,
Orthoreovirus mammalis, 10 segments L1-L3/M1-M3/S1-S4), in MagNA_1, with the 10
segments laid end-to-end into ONE continuous track. Reads are coloured by the STRAIN
reference they were assigned to — Dearing strain (taxid 10886, NC_077xxx) vs the
other Reovirus-3 reference set (taxid 538123, NC_013xxx). Shows the same genuine
species' reads splitting across two near-identical strain references segment by
segment (the strain-split that LCA / species-level rollup resolves), and that some
segments are filled while others are empty (partial-segment detection is normal).

  conda run -n shotgun_virome python scripts/benchmark/fig_segment_coverage.py
"""
from __future__ import annotations
import sys, os
sys.path.insert(0, "/home/share/programs/nexvirome/notebooks")
sys.path.insert(0, "/home/share/programs/nexvirome/scripts")
import numpy as np, pandas as pd, sqlite3
import matplotlib
matplotlib.use("Agg")
import paper_style  # noqa: F401  (sets journal font/style)
import matplotlib.pyplot as plt
from benchmark_utils import save_figure

NX = "/home/share/programs/nexvirome"
DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
KITD = "/home/share/programs/vshot/result_kit/newdb_unmasked_mmseqs2"
COLS = ["query", "target", "pident", "alnlen", "mismatch", "gapopen",
        "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen", "taxid"]
SAMPLE = "MagNA_1"
_c = sqlite3.connect(DB)

# segment order L1..S4; for each segment, the two strain references (Dearing 10886 /
# 538123) — coverage is computed per reference and laid on the shared segment axis.
SEG_ORDER = ["L1", "L2", "L3", "M1", "M2", "M3", "S1", "S2", "S3", "S4"]
# Greyscale contrast (no colour). Top panel (Dearing) = WHITE fill with a dark
# outline; bottom panel (the other reference set) = diagonal hatch over light
# grey. Both readable on B/W print.
STRAINS = {10886: ("Mammalian orthoreovirus 3 Dearing (taxid 10886)",
                   {"facecolor": "white", "hatch": None,
                    "edgecolor": "#1a1a1a", "linewidth": 1.0}),
           538123: ("Mammalian orthoreovirus 3 (taxid 538123)",
                   {"facecolor": "#cccccc", "hatch": "////",
                    "edgecolor": "#1a1a1a", "linewidth": 0.6})}


def seg_refs():
    """segment_name -> {strain_taxid: (accession, length)} for Reovirus 3."""
    out = {}
    rows = _c.execute(
        "SELECT accession,taxid,length,segment_name FROM refseq_sequences "
        "WHERE taxid IN (10886,538123) AND segment_name IS NOT NULL")
    for acc, tax, L, seg in rows:
        out.setdefault(seg, {})[int(tax)] = (acc.split(".")[0], int(L))
    return out


def depth(g, L):
    d = np.zeros(L + 1, dtype=np.int32)
    for s, e in zip(g.tstart, g.tend):
        s = max(1, min(int(s), L)); e = max(1, min(int(e), L))
        d[s - 1] += 1; d[e] -= 1
    return np.cumsum(d[:L])


WIN = 100   # bp window for average identity


MIN_HITS = 5   # only show identity for windows with enough reads (drops noisy tails)


def identity_windows(g, L):
    """mean pident per WIN-bp window; only windows with >=MIN_HITS reads (else NaN)."""
    nb = max(1, L // WIN)
    centers = np.arange(nb) * WIN + WIN / 2
    mid = ((g["tstart"] + g["tend"]) / 2).astype(int).clip(0, L - 1)
    wi = (mid // WIN).clip(0, nb - 1)
    s = pd.DataFrame({"w": wi.values, "pid": g["pident"].values})
    grp = s.groupby("w")["pid"]
    mean, cnt = grp.mean(), grp.size()
    vals = np.full(nb, np.nan)
    for w in mean.index:
        if cnt[w] >= MIN_HITS:
            vals[int(w)] = mean[w]
    return centers, vals


def main():
    df = pd.concat([pd.read_csv(f"{KITD}/{SAMPLE}_{r}.result", sep="\t", header=None, names=COLS)
                    for r in ("R1", "R2")], ignore_index=True)
    df = df[(df["pident"] >= 0.85) & (df["alnlen"] / df["qlen"] >= 0.5)]
    refs = seg_refs()

    # two stacked panels: top = Dearing (10886), bottom = 538123. No text, legend only.
    GAP = 200            # within-group gap (between segments of the same size class)
    GROUP_GAP = 900      # extra gap at L→M and M→S boundaries to show the three groups
    # shared segment layout + shared y-limit
    layout = []   # (seg, x0, Lmax)
    x0 = 0
    prev_group = None
    for seg in SEG_ORDER:
        if seg not in refs:
            continue
        cur_group = seg[0]               # 'L' / 'M' / 'S'
        if prev_group is not None and cur_group != prev_group:
            x0 += GROUP_GAP - GAP        # bump beyond the within-gap that was already added
        Lmax = max(L for _, L in refs[seg].values())
        layout.append((seg, x0, Lmax, cur_group))
        x0 += Lmax + GAP
        prev_group = cur_group
    xtot = x0
    # group spans (for L/M/S header bands)
    group_spans = {}
    for seg, sx, Lmax, gr in layout:
        if gr not in group_spans:
            group_spans[gr] = [sx, sx + Lmax]
        else:
            group_spans[gr][1] = sx + Lmax
    ymax = 1
    for seg, sx, Lmax, _gr in layout:
        for tax, (acc, L) in refs[seg].items():
            g = df[df["target"].astype(str).str.startswith(acc)]
            ymax = max(ymax, depth(g, L).max())

    ytop = int(ymax * 1.05)
    seg_centers = [(sx + Lmax / 2, seg) for seg, sx, Lmax, _gr in layout]
    handles, labels = [], []
    id_handle = None
    fig, axes = plt.subplots(2, 1, figsize=(15, 5.6), sharex=True)
    for ax, tax in zip(axes, (10886, 538123)):
        label, style = STRAINS[tax]
        ax2 = ax.twinx()                       # right axis = avg identity
        for seg, sx, Lmax, _gr in layout:
            if tax in refs[seg]:
                acc, L = refs[seg][tax]
                g = df[df["target"].astype(str).str.startswith(acc)]
                ax.fill_between(np.arange(L) + sx, 0, depth(g, L), **style)
                if len(g):                     # identity: open dots, B/W
                    cx, cy = identity_windows(g, L)
                    m = ~np.isnan(cy)
                    ax2.plot((cx + sx)[m], cy[m] * 100, color="#1a1a1a", lw=2.0,
                             alpha=0.85, zorder=3)
                    sc = ax2.scatter((cx + sx)[m], cy[m] * 100, s=18,
                                     facecolor="#1a1a1a", edgecolor="#1a1a1a",
                                     linewidth=0.5, zorder=4)
                    if id_handle is None:
                        id_handle = sc
            ax.axvline(sx + Lmax + GAP / 2, color="#dddddd", lw=0.6, ls="--")
        sh = ax.fill_between([], [], [], **style)
        handles.append(sh); labels.append(label)
        ax.set_xlim(-GAP, xtot); ax.set_ylim(0, ytop)
        ax.set_yticks([0, ytop])
        ax.set_ylabel("read depth", fontsize=16, labelpad=2)
        ax.tick_params(axis="y", labelsize=14, pad=2)
        ax2.set_ylim(80, 100); ax2.set_yticks([80, 85, 90, 95, 100])
        ax2.set_ylabel("% identity", fontsize=16, labelpad=3)
        ax2.tick_params(axis="y", labelsize=14)
        for s in ("top", "right", "bottom", "left"):
            ax.spines[s].set_visible(True); ax.spines[s].set_color("#333333")
            ax.spines[s].set_linewidth(1.0)
        ax.set_xticks([c for c, _ in seg_centers])
        ax.set_xticklabels([s for _, s in seg_centers], fontsize=15)

    # group header band above the top panel (L / M / S)
    group_labels = {"L": "L  (large)", "M": "M  (medium)", "S": "S  (small)"}
    ymax_top = axes[0].get_ylim()[1]
    for gr, (xL, xR) in group_spans.items():
        cx = (xL + xR) / 2
        axes[0].annotate(group_labels.get(gr, gr), xy=(cx, ymax_top), xytext=(0, 8),
                         textcoords="offset points", ha="center", va="bottom",
                         fontsize=15, fontweight="bold")
        axes[0].plot([xL, xR], [ymax_top * 1.01] * 2, color="#1a1a1a",
                     lw=1.2, clip_on=False)

    if id_handle is not None:
        handles.append(id_handle)
        labels.append("avg % identity (100 bp window, ≥5 reads)")
    axes[1].set_xlabel("Reovirus 3 segments (L → M → S)", fontsize=16)
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=14,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout(rect=[0, 0.06, 1, 0.94])
    save_figure(fig, "FigS6_segment_coverage_reovirus")
    plt.close(fig)
    print("saved -> FigS6_segment_coverage_reovirus (B/W, L|M|S grouped, identity 80-100)")


if __name__ == "__main__":
    main()
