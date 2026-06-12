#!/usr/bin/env python3
"""
Fig2 composition — single plot, absolute mapped reads (linear y).

One figure, x = 5 tools each with its 4 KIT samples side by side (20 bars).
Each bar stacked = 6 GT species (Set2 pastel) + host (Homo sapiens, distinct
colour) + FP non-host (grey hatch). No 100% normalization — bar height is the
real number of reads each tool assigned, so:
  - native tools show host (~95%) dominating, GT/FP as thin bands below;
  - NexVirome (viral-only DB) has host=0, so its bar is essentially all GT.

Abundance = raw mapped reads (NOT TPM): host/bacteria/viruses are mixed, so
genome-length normalization is meaningless here.

Output: result_260605/fig2/fig2_composition_abs.{png,pdf}
Run: conda run -n shotgun_virome python result_260605/fig2/render_fig2_composition_abs.py
"""
from __future__ import annotations
import os, sys, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["pdf.fonttype"] = 42

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/scripts"); sys.path.insert(0, f"{NX}/notebooks")
sys.path.insert(0, f"{NX}/scripts/benchmark")
import three_strategy_breadth_lca as TS
from virome_classifier.alignment.filters.filter import MaskingFilter
from virome_classifier.taxonomy import TaxonomyDB
from benchmark_utils import (GROUND_TRUTH, parse_kreport, kreport_to_species_counts,
                             parse_metabuli_report, metabuli_to_species_counts)

DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
MASK = f"{NX}/result_260605/mask/mask_v3_full.bed"
NAT = f"{NX}/resources/db_20260525/external_native"
CACHE = "/tmp/hq_cache"
SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
CUT = 0.01
HOST = {9606}
OUT = f"{NX}/result_260605/fig2"

TS.DB = DB; TS._TAX_PATH = DB; TS._TLEN = TS._load_tlen()
TAX = TaxonomyDB.from_sqlite(DB)

GT_GROUPS = [(g["name"], set(g["taxids"]), g["genus"]) for g in GROUND_TRUTH]
GT_TAXIDS = set().union(*[s for _, s, _ in GT_GROUPS])
SPECIES = ["Adenovirus 40", "HHV-5 (CMV)", "Human RSV", "Influenza B", "Reovirus 3", "Zika virus"]
# Set2 pastel for the 6 GT species
SP_COLORS = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F"]
HOST_COLOR = "#7A7A7A"      # host (Homo sapiens) — solid dark grey
FP_COLOR = "#CFCFCF"        # FP non-host — light grey, hatched
TOOLS = ["NexVirome", "Ganon", "Kraken2", "Metabuli", "Phanta"]


def comp(counts):
    """{cat: reads} — 6 GT species + 'host' + 'FP' (non-host, non-GT)."""
    out = {n: sum(counts.get(t, 0) for t in tx) for n, tx, _ in GT_GROUPS}
    out["host"] = sum(counts.get(t, 0) for t in HOST)
    out["FP"] = sum(r for t, r in counts.items() if t not in GT_TAXIDS and t not in HOST)
    return out


def nexvirome():
    bed = pd.read_csv(MASK, sep="\t", header=None, usecols=[0, 1, 2], names=["target", "start", "end"])
    mf = MaskingFilter.from_dataframe(bed)
    res = {}
    for s in SAMPLES:
        df = TS._add_tlen(pd.read_parquet(f"{CACHE}/{s}.parquet"))
        best = TS._best_hit(df)
        refs = TS._refs_at(TS._breadth_by_ref(best, mf), CUT)
        res[s] = comp(TS._besthit_counts(best[best["target"].isin(refs)]))
    return res


def ganon_native(s):
    p = f"{NAT}/kit_ganon/{s}.tre"
    if not os.path.exists(p):
        return {}
    df = pd.read_csv(p, sep="\t", header=None,
                     names=["rank", "taxid", "lineage", "name", "u", "sh", "ch", "cum", "pct"])
    sp = df[df["rank"] == "species"]
    return {int(r.taxid): int(r.cum) for r in sp.itertuples() if int(r.cum) > 0}


def tool_native(tool, s):
    if tool == "Ganon":
        return ganon_native(s)
    if tool == "Kraken2":
        return kreport_to_species_counts(parse_kreport(f"{NAT}/kit_kraken2/{s}.kreport"), filter_virus=False)
    if tool == "Metabuli":
        cand = glob.glob(f"{NAT}/kit_metabuli/{s}/*report.tsv")
        return metabuli_to_species_counts(parse_metabuli_report(cand[0]), filter_virus=False) if cand else {}
    if tool == "Phanta":
        for pat in [f"{NAT}/kit_phanta/results/classification/{s}.krak.report_bracken_species.filtered",
                    f"{NAT}/kit_phanta/results/classification/{s}.krak.report.filtered"]:
            if os.path.exists(pat):
                return kreport_to_species_counts(parse_kreport(pat), filter_virus=False)
        return {}
    return {}


def main():
    data = {"NexVirome": nexvirome()}
    for tool in ["Ganon", "Kraken2", "Metabuli", "Phanta"]:
        data[tool] = {s: comp(tool_native(tool, s)) for s in SAMPLES}

    cats = SPECIES + ["host", "FP"]
    colors = SP_COLORS + [HOST_COLOR, FP_COLOR]
    hatches = [""] * 6 + ["", "//"]

    # x layout: tool groups, 4 sample bars per group
    n_s = len(SAMPLES)
    group_w = 1.0
    bar_w = group_w / (n_s + 0.6)
    fig, ax = plt.subplots(figsize=(15, 6.2))
    xticks, xticklabels = [], []
    for gi, tool in enumerate(TOOLS):
        base = gi * (group_w + 0.5)
        for si, s in enumerate(SAMPLES):
            x = base + si * bar_w
            d = data[tool].get(s, {})
            bottom = 0.0
            for cat, col, ht in zip(cats, colors, hatches):
                v = d.get(cat, 0)
                ax.bar(x, v, bar_w * 0.92, bottom=bottom, color=col, edgecolor="white",
                       linewidth=0.3, hatch=ht,
                       label=cat if (gi == 0 and si == 0) else None)
                bottom += v
            xticks.append(x)
            xticklabels.append(s.split("_")[1])  # 1/2
        # group label
        ax.text(base + (n_s - 1) * bar_w / 2, -0.07, tool, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=12, weight="bold")

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=8)
    ax.set_ylabel("Mapped reads", fontsize=12)
    ax.tick_params(axis="x", length=0)
    # legend: GT species + host + FP
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in SP_COLORS]
    handles += [plt.Rectangle((0, 0), 1, 1, color=HOST_COLOR),
                plt.Rectangle((0, 0), 1, 1, facecolor=FP_COLOR, hatch="//")]
    labels = SPECIES + ["host (Homo sapiens)", "FP (non-host)"]
    ax.legend(handles, labels, title="taxon", bbox_to_anchor=(1.01, 1),
              loc="upper left", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/fig2_composition_abs.{ext}", dpi=300, bbox_inches="tight")
    print(f"-> {OUT}/fig2_composition_abs.png , .pdf")


if __name__ == "__main__":
    main()
