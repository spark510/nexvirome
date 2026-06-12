#!/usr/bin/env python3
"""
Fig2 relative-abundance bar WITHOUT the viral-only filter — host (Homo sapiens,
9606) removed, but everything else KEPT (so non-target bacteria/fungi/other
viruses all count as FP). Shows what each classifier reports across the whole
non-host space, not just the viral subtree.

Composition = 6 GT species (color) + lumped FP-grey, raw read %, sum=100%
(raw reads, NOT TPM: mixing bacteria+viruses makes genome-length TPM meaningless).
NexVirome + 4 external tools, native DB (same-DB dropped). species level only.

Output: result_260605/fig2/fig2_relabund_bar_nohostfilter.{png,pdf}
Run: conda run -n shotgun_virome python scripts/benchmark/render_fig2_relabund_nohostfilter.py
"""
from __future__ import annotations
import os, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

matplotlib.rcParams["font.family"] = "Arial"
matplotlib.rcParams["pdf.fonttype"] = 42

NX = "/home/share/programs/nexvirome"
import sys
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
# display label only (Roche = MagNA Pure kit); internal id MagNA_1 unchanged.
SAMPLE_DISPLAY = {"MagNA_1": "Roche 1", "MagNA_2": "Roche 2",
                  "Qiagen_1": "Qiagen 1", "Qiagen_2": "Qiagen 2"}
CUT = 0.01
HOST = {9606}                          # Homo sapiens (host) — removed
OUT = f"{NX}/result_260605/fig2"

TS.DB = DB; TS._TAX_PATH = DB; TS._TLEN = TS._load_tlen()
TAX = TaxonomyDB.from_sqlite(DB)

GT_GROUPS = [(g["name"], set(g["taxids"]), g["genus"]) for g in GROUND_TRUTH]
GT_TAXIDS = set().union(*[s for _, s, _ in GT_GROUPS])
SPECIES = ["Adenovirus 40", "HHV-5 (CMV)", "Human RSV", "Influenza B", "Reovirus 3", "Zika virus"]
SP_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#CCB974"]
FP_COLOR = "#9e9e9e"

BAR_ORDER = ["NexVirome", "Ganon", "Kraken2", "Metabuli", "Phanta"]


def drop_host(c):
    return {t: r for t, r in (c or {}).items() if t not in HOST}


def composition(counts):
    """6 GT species + lumped FP as ABSOLUTE mapped-read counts (host removed).
    NOT normalized to 100% — the y-axis is the real number of reads each tool
    assigned, so per-tool detection magnitude and the true TP:FP ratio show."""
    out = {}
    for name, tx, _ in GT_GROUPS:
        out[name] = sum(counts.get(t, 0) for t in tx)
    out["FP"] = sum(r for t, r in counts.items() if t not in GT_TAXIDS)
    return out


def nexvirome():
    bed = pd.read_csv(MASK, sep="\t", header=None, usecols=[0, 1, 2], names=["target", "start", "end"])
    mf = MaskingFilter.from_dataframe(bed)
    comps = []
    for s in SAMPLES:
        df = TS._add_tlen(pd.read_parquet(f"{CACHE}/{s}.parquet"))
        best = TS._best_hit(df)
        refs = TS._refs_at(TS._breadth_by_ref(best, mf), CUT)
        c = drop_host(TS._besthit_counts(best[best["target"].isin(refs)]))
        comps.append(composition(c))
    return comps


def ganon_native(s):
    p = f"{NAT}/kit_ganon/{s}.tre"
    if not os.path.exists(p):
        return None
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
        return metabuli_to_species_counts(parse_metabuli_report(cand[0]), filter_virus=False) if cand else None
    if tool == "Phanta":
        for pat in [f"{NAT}/kit_phanta/results/classification/{s}.krak.report_bracken_species.filtered",
                    f"{NAT}/kit_phanta/results/classification/{s}.krak.report.filtered"]:
            if os.path.exists(pat):
                return kreport_to_species_counts(parse_kreport(pat), filter_virus=False)
        return None
    return None


def main():
    comp_by_tool = {"NexVirome": nexvirome()}
    for tool in ["Ganon", "Kraken2", "Metabuli", "Phanta"]:
        comps = []
        for s in SAMPLES:
            c = drop_host(tool_native(tool, s))
            if c:
                comps.append(composition(c))
        comp_by_tool[tool] = comps

    # ABSOLUTE mapped reads, per sample (depth differs 7-15x across samples, so a
    # 4-panel layout instead of averaging into one bar).
    cats = SPECIES + ["FP"]
    x = np.arange(len(BAR_ORDER))
    fig, axes = plt.subplots(1, 4, figsize=(18, 5.4), sharey=False)
    for si, sample in enumerate(SAMPLES):
        ax = axes[si]
        # value per tool for this sample
        def cell(tool, cat):
            lst = comp_by_tool[tool]
            return lst[si].get(cat, 0) if si < len(lst) else 0
        bottom = np.zeros(len(BAR_ORDER))
        for sp, col in zip(SPECIES, SP_COLORS):
            vals = np.array([cell(t, sp) for t in BAR_ORDER], dtype=float)
            ax.bar(x, vals, bottom=bottom, color=col, edgecolor="white", linewidth=0.4,
                   label=sp if si == 3 else None)
            bottom += vals
        fp = np.array([cell(t, "FP") for t in BAR_ORDER], dtype=float)
        ax.bar(x, fp, bottom=bottom, color=FP_COLOR, edgecolor="white", linewidth=0.4,
               label="FP (non-host)" if si == 3 else None, hatch="//")
        ax.axvline(0.5, color="black", lw=0.8, alpha=0.6)
        ax.set_xticks(x); ax.set_xticklabels(BAR_ORDER, fontsize=9, rotation=30, ha="right")
        ax.set_title(SAMPLE_DISPLAY.get(sample, sample), fontsize=11)
        if si == 0:
            ax.set_ylabel("Mapped reads (non-host)", fontsize=11)
    axes[3].legend(title="taxon (TP) / FP", bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=9)
    fig.tight_layout()
    for ext in ("png", "eps"):
        fig.savefig(f"{OUT}/fig2_relabund_bar_nohostfilter.{ext}", dpi=300, bbox_inches="tight")
    # digest
    pd.set_option("display.width", 220)
    print("ABSOLUTE non-host mapped reads (per sample, TP6 sum vs FP):")
    for si, sample in enumerate(SAMPLES):
        print(f"\n[{sample}]")
        for t in BAR_ORDER:
            d = comp_by_tool[t][si] if si < len(comp_by_tool[t]) else {}
            gt6 = sum(d.get(s, 0) for s in SPECIES); fp = d.get("FP", 0)
            print(f"  {t:10}: GT6={gt6:>9,}  FP={fp:>9,}")
    print(f"\n-> {OUT}/fig2_relabund_bar_nohostfilter.png , .pdf")


if __name__ == "__main__":
    main()
