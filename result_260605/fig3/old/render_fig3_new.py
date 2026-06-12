#!/usr/bin/env python3
"""
Fig 3 (healthy DNA-vs-RNA) COMBINED — self-contained vector render.
  A  alpha diversity (Richness / Shannon / inv.Simpson), DNA vs RNA per tool
  B  Bray-Curtis PCoA per tool, DNA vs RNA, points + 95% group ellipses (NEW)
  C  genus composition per tool: shared genera (>=3 tools) coloured + each tool's
     TOP-3 non-shared genera stacked on top (black border, named above)  (NEW)

Self-contained (does not import the broken fig3_style/panel_* paths); reads
fig3/tables/ and writes fig3/Fig3_combined.{eps,png}.
Run: conda run -n shotgun_virome python result_260605/fig3/code/render_fig3_new.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams.update({
    "font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 1.2,
})
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch, Ellipse, Rectangle

T = "/home/share/programs/nexvirome/result_260605/fig3/tables"
OUT = "/home/share/programs/nexvirome/result_260605/fig3"
TOOLS = ["NexVirome", "Ganon", "Phanta", "Kraken2", "Metabuli"]
CC = {"dna": "#E69F00", "rna": "#009E73"}            # orange=DNA, green=RNA
INDICES = [("richness", "Richness"), ("shannon", "Shannon H'"), ("invsimpson", "inv. Simpson")]
SHARED_PAL = ["#4e79a7", "#59a14f", "#9c755f", "#edc948", "#b07aa1", "#76b7b2",
              "#ff9da7", "#f28e2b", "#bab0ac", "#86bcb6", "#d37295", "#a0cbe8"]
NS3 = ["#cb181d", "#fb6a4a", "#fcae91"]              # non-shared rank1->3
MIN_TOOLS, TOPN, NS_N = 3, 12, 3
KEY, NAMECOL, VAL = "genus_taxid", "genus_name", "rel_abund"
ABBR = {"Human endogenous retrovirus K": "HERV-K"}

def short(nm):
    """Consistent label: unnamed OTU for taxid_*, else abbreviate the leading
    genus word of any multi-word name to its initial (so all 'Streptococcus
    phage X' read 'S. phage X', etc.)."""
    nm = ABBR.get(nm, nm)
    s = str(nm)
    if s.startswith("taxid_"):
        return "unnamed OTU"
    parts = s.split()
    if len(parts) >= 2 and parts[0][:1].isupper():
        return parts[0][0] + ". " + " ".join(parts[1:])
    return s

def stars(p):
    if p != p:
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""

def panel_letter(fig, lab, y):
    fig.text(0.005, y, lab, fontsize=26, weight="bold", va="top", ha="left")


# ---------------------------------------------------------------- panel A
def cell_a(ax, alpha, mw, idx):
    xpos, lab, tmax = [], [], {}
    for ti, tool in enumerate(TOOLS):
        vmax = 0.0
        for gi, coh in enumerate(["dna", "rna"]):
            vals = alpha[(alpha.tool == tool) & (alpha.cohort == coh)][idx].dropna().values
            x = ti*2.6 + gi
            if len(vals):
                ax.boxplot([vals], positions=[x], widths=0.8, patch_artist=True, showfliers=False,
                           medianprops=dict(color="white", lw=1.8),
                           boxprops=dict(facecolor=CC[coh], edgecolor="black", lw=0.8),
                           whiskerprops=dict(color="black"), capprops=dict(color="black"))
                jit = np.random.default_rng(ti*10+gi).normal(x, 0.06, len(vals))
                ax.scatter(jit, vals, s=10, color="black", alpha=0.4, edgecolor="none", zorder=3)
                vmax = max(vmax, float(np.max(vals)))
        tmax[tool] = vmax; xpos.append(ti*2.6+0.5); lab.append(tool)
    ax.set_xticks(xpos); ax.set_xticklabels(lab, rotation=90, ha="center", fontsize=18)
    ax.tick_params(axis="y", labelsize=15)   # match Panel B y-tick size
    ymax = max(tmax.values()) if tmax else 1
    for ti, tool in enumerate(TOOLS):
        row = mw[(mw.tool == tool) & (mw["index"] == idx)]
        if row.empty:
            continue
        st = stars(float(row.p.iloc[0]))
        if not st:
            continue
        x0, x1 = ti*2.6, ti*2.6+1
        y = tmax[tool] + ymax*0.06
        ax.plot([x0, x0, x1, x1], [y, y+ymax*0.02, y+ymax*0.02, y], color="black", lw=1.1)
        ax.text((x0+x1)/2, y+ymax*0.03, st, ha="center", va="bottom", fontsize=24)
    ax.set_ylim(top=ymax*1.25)


# ---------------------------------------------------------------- panel B (+ ellipses)
def _ellipse(ax, x, y, color, n_std=2.0):
    if len(x) < 3:
        return
    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)) or cov[0, 0] <= 0 or cov[1, 1] <= 0:
        return
    pear = cov[0, 1] / np.sqrt(cov[0, 0]*cov[1, 1])
    rx, ry = np.sqrt(1+pear), np.sqrt(1-pear)
    e = Ellipse((0, 0), width=rx*2, height=ry*2, facecolor="none",
                edgecolor=color, lw=1.8, zorder=1)
    tr = (mtransforms.Affine2D().rotate_deg(45)
          .scale(np.sqrt(cov[0, 0])*n_std, np.sqrt(cov[1, 1])*n_std)
          .translate(np.mean(x), np.mean(y)))
    e.set_transform(tr + ax.transData); ax.add_patch(e)

def cell_b(ax, pcoa, perm, tool):
    ax.set_box_aspect(1)
    d = pcoa[pcoa.tool == tool]
    if d.empty:
        ax.text(0.5, 0.5, "(no viral)", ha="center", va="center", transform=ax.transAxes,
                fontsize=20, color="grey")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel(tool, fontsize=21, labelpad=8); return
    for coh in ("dna", "rna"):
        g = d[d.cohort == coh]
        _ellipse(ax, g.PC1.values, g.PC2.values, CC[coh])
        ax.scatter(g.PC1, g.PC2, s=30, color=CC[coh], alpha=0.80, edgecolor="white",
                   linewidth=0.4, zorder=3, label=coh.upper())
    v1, v2 = d.var1.iloc[0]*100, d.var2.iloc[0]*100
    pr = perm[perm.tool == tool]
    ptxt = f"PERM p={pr.p.iloc[0]:.3f}" if len(pr) else ""
    ax.set_xlabel(f"PC1 ({v1:.0f}%)", fontsize=18); ax.set_ylabel(f"PC2 ({v2:.0f}%)", fontsize=18)
    ax.tick_params(labelsize=15)
    ax.text(0.04, 0.96, ptxt, transform=ax.transAxes, ha="left", va="top", fontsize=16.5)
    ax.text(0.5, -0.30, tool, transform=ax.transAxes, ha="center", va="top", fontsize=21)


# ---------------------------------------------------------------- panel C (new)
def comp_of(gn, tool, coh):
    sub = gn[(gn.tool == tool) & (gn.cohort == coh)]
    if sub.empty:
        return pd.Series(dtype=float)
    comp = sub.groupby(KEY)[VAL].sum() / max(sub["sample"].nunique(), 1)
    s = comp.sum()
    return comp / s if s else comp

def cell_c(ax, gn, tool, shared, shared_color, ns_of, name_of):
    ax.set_box_aspect(1.05)
    ns = ns_of[tool]
    for xi, coh in enumerate(["dna", "rna"]):
        comp = comp_of(gn, tool, coh)
        bottom = 0.0
        for x in shared:
            v = comp.get(x, 0.0)
            ax.bar(xi, v*100, width=0.8, bottom=bottom*100, color=shared_color[x],
                   edgecolor="white", linewidth=0.2); bottom += v
        ns_total = sum(comp.get(x, 0.0) for x in ns)
        other = max(0.0, 1.0 - bottom - ns_total)
        ax.bar(xi, other*100, width=0.8, bottom=bottom*100, color="#dddddd",
               edgecolor="white", linewidth=0.2); bottom += other
        for j in range(len(ns)-1, -1, -1):                # non-shared on top, rank1 highest
            v = comp.get(ns[j], 0.0)
            if v <= 0:
                continue
            ax.bar(xi, v*100, width=0.8, bottom=bottom*100, color=NS3[j],
                   edgecolor="black", linewidth=1.4); bottom += v
    ax.set_xlim(-0.7, 1.7); ax.set_ylim(0, 100)
    ax.tick_params(axis="y", labelsize=15)   # match Panel B y-tick size
    ax.set_xticks([0, 1]); ax.set_xticklabels(["DNA", "RNA"], fontsize=16.5)
    # tool name: match Panel B's tool label (fontsize 21, not bold)
    ax.text(0.5, -0.17, tool, transform=ax.transAxes, ha="center", va="top",
            fontsize=21)
    if ns:                                             # centred mini-legend above the panel
        hs = [Patch(facecolor=NS3[j], edgecolor="black", linewidth=1.1) for j in range(len(ns))]
        ls = [short(name_of.get(x, x)) for x in ns]
        ax.legend(hs, ls, loc="lower center", bbox_to_anchor=(0.5, 1.08),
                  frameon=False, fontsize=16.5, handlelength=1.5, handletextpad=0.5,
                  labelspacing=0.30, borderaxespad=0.0)


# ---------------------------------------------------------------- main
def main():
    alpha = pd.read_csv(f"{T}/fig3_alpha.csv")
    mw = pd.read_csv(f"{T}/fig3_alpha_dna_vs_rna.csv")
    pcoa = pd.read_csv(f"{T}/fig3_pcoa_pooled.csv")
    perm = (pd.read_csv(f"{T}/fig3_permanova_dna_vs_rna.csv")
            if os.path.exists(f"{T}/fig3_permanova_dna_vs_rna.csv")
            else pd.DataFrame(columns=["tool", "p"]))
    gn = pd.read_csv(f"{T}/fig3_genus_long.csv"); gn = gn[gn.tool.isin(TOOLS)]
    # Fig3 = Healthy/Control cohort only (group manifest); drop Asthma samples so
    # Panel C composition matches the Control alpha/PCoA panels.
    GROUPS = "/home/share/programs/nexvirome/paper/figures/Fig3/source_data"
    keep = set()
    for f in ("dna_asthma_groups.csv", "rna_asthma_groups.csv"):
        g = pd.read_csv(f"{GROUPS}/{f}")
        keep |= set(g.loc[g["group"] == "Control", "sample"])
    gn = gn[gn["sample"].isin(keep)]
    name_of = gn.drop_duplicates(KEY).set_index(KEY)[NAMECOL].to_dict()

    g = gn[gn.reads > 0]
    ntool = g.groupby(KEY)["tool"].nunique()
    tot = gn.groupby(KEY)["reads"].sum()
    shared = list(tot.loc[ntool[ntool >= MIN_TOOLS].index].sort_values(ascending=False).head(TOPN).index)
    shared_set = set(shared)
    shared_color = {x: SHARED_PAL[i % len(SHARED_PAL)] for i, x in enumerate(shared)}
    ns_of = {}
    for t in TOOLS:
        r = g[g.tool == t].groupby(KEY)["reads"].sum()
        r = r[[x for x in r.index if x not in shared_set]]
        ns_of[t] = list(r.sort_values(ascending=False).head(NS_N).index)

    fig = plt.figure(figsize=(18.0, 17.8))
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1.05, 1.0, 1.30],
                  hspace=0.52, left=0.055, right=0.99, top=0.965, bottom=0.07)

    gsa = gs[0].subgridspec(1, len(INDICES), wspace=0.28)
    for gi, (idx, ilab) in enumerate(INDICES):
        ax = fig.add_subplot(gsa[0, gi]); cell_a(ax, alpha, mw, idx); ax.set_ylabel(ilab, fontsize=18)
    fig.legend([Patch(facecolor=CC["dna"]), Patch(facecolor=CC["rna"])], ["DNA", "RNA"],
               loc="center", ncol=2, frameon=False, bbox_to_anchor=(0.5, 0.672), fontsize=18)

    gsb = gs[1].subgridspec(1, len(TOOLS), wspace=0.32)
    for ti, tool in enumerate(TOOLS):
        cell_b(fig.add_subplot(gsb[0, ti]), pcoa, perm, tool)

    gsc = gs[2].subgridspec(1, len(TOOLS), wspace=0.18)
    for ti, tool in enumerate(TOOLS):
        ax = fig.add_subplot(gsc[0, ti])
        cell_c(ax, gn, tool, shared, shared_color, ns_of, name_of)
        ax.set_ylabel("Rel. abundance (%)" if ti == 0 else "", fontsize=21)
        if ti != 0:
            ax.tick_params(labelleft=False)
    # leading entry carries the "shared genus" meaning so no separate title is
    # needed (a centred title overlapped the per-tool names above the legend).
    handles = [Patch(facecolor="none", edgecolor="none")] + \
              [Patch(facecolor=shared_color[x]) for x in shared] + \
              [Patch(facecolor=NS3[0], edgecolor="black", linewidth=1.4), Patch(facecolor="#dddddd")]
    labels = ["shared genus (≥3 tools):"] + \
             [short(name_of.get(x, x)) for x in shared] + \
             ["top-3 non-shared/tool (named above)", "Other"]
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.004), frameon=False,
               fontsize=16.5, ncol=5, columnspacing=1.2, handletextpad=0.5)

    panel_letter(fig, "A", 0.985)
    panel_letter(fig, "B", 0.635)
    panel_letter(fig, "C", 0.345)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/Fig3_combined.{ext}", dpi=300)
    print(f"-> {OUT}/Fig3_combined.png , .eps")
    for t in TOOLS:
        print(f"  {t:10s} top-3 non-shared: " + " | ".join(short(name_of.get(x, x)) for x in ns_of[t]))


if __name__ == "__main__":
    main()
