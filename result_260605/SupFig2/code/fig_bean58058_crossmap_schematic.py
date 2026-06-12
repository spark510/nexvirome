#!/usr/bin/env python3
"""
Supplementary figure: NC_032111.1 (BeAn 58058 virus) megaBLAST cross-mapping
to T2T-CHM13v2.0 — viral mask candidate schematic.

Top panel    : viral genome (NC_032111.1, 163,005 bp). One short 8,262-8,658
               window (397 bp, 0.24 %) is the *only* region that ever aligns
               to the human reference. This is the proposed mask target.
Bottom panel : human T2T-CHM13v2.0 chromosomes (chr1-22, X, Y) drawn at true
               relative length, with vertical tick marks at every BLAST hit.
               Hits exist on every chromosome (~280,000 HSPs total), with
               typical alignment length ~300 bp and ~85-93 % identity — the
               classic signature of a dispersed repeat element (SINE/LINE)
               buried inside the viral reference deposit.

Source data  : pre-computed megablast XML
               resources/db/ncbi/ncbi_human_genome/refseq_viral_to_hgt2t.xml
               (one <Iteration> per viral query; this script scans for the
               NC_032111.1 iteration only.)

This is the result_260605/SupFig2 copy of the figure script. It saves directly
into result_260605/SupFig2/ as SupFig2_bean58058_crossmap.{png,tiff,eps} so the
folder is self-contained. (The original lives at
scripts/benchmark/fig_bean58058_crossmap_schematic.py and saves to paper/figures/.)

NOTE on masking: the red/black window (8,262–8,658 bp) is the exact viral-side
cross-mapping window, computed from the human-HSP query coordinates (min query-from
to max query-to; see build_supfig2_hsp_table.py). In production the whole NC_032111.1 reference is
masked (mask_v3_full.bed, `prod` tag, 0–163005) — see SupFig2_result.md.

Run:
  conda run -n shotgun_virome python \
      result_260605/SupFig2/code/fig_bean58058_crossmap_schematic.py
Outputs (result_260605/SupFig2/):
  SupFig2_bean58058_crossmap.{png,tiff,eps}
"""
from __future__ import annotations
import os, re, sys
import xml.etree.ElementTree as ET
from collections import defaultdict

sys.path.insert(0, "/home/share/programs/nexvirome/notebooks")
sys.path.insert(0, "/home/share/programs/nexvirome/scripts/benchmark")
import matplotlib
matplotlib.use("Agg")
import paper_style  # noqa: F401  (sets journal font/style)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from benchmark_utils import save_figure

NX       = "/home/share/programs/nexvirome"
XML_PATH = f"{NX}/resources/db/ncbi/ncbi_human_genome/refseq_viral_to_hgt2t.xml"
QUERY_ACC = "NC_032111.1"          # BeAn 58058 virus, RefSeq
QUERY_NAME = "BeAn 58058 virus"
FIG_NAME = "SupFig2_bean58058_crossmap"
OUT_DIR = f"{NX}/result_260605/SupFig2"

CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]


def extract_iteration(xml_path: str, accession: str):
    """Stream the multi-document XML and pull the single <Iteration> block
    whose <Iteration_query-def> starts with the requested accession.
    Returns (root_element, query_len).
    """
    buf, in_iter, found = [], False, False
    with open(xml_path) as f:
        for line in f:
            if "<Iteration>" in line:
                buf, in_iter = [line], True
                continue
            if in_iter:
                buf.append(line)
                if f"<Iteration_query-def>{accession}" in line:
                    found = True
                if "</Iteration>" in line:
                    if found:
                        root = ET.fromstring("".join(buf))
                        return root, int(root.findtext("Iteration_query-len", "0"))
                    buf, in_iter = [], False
    raise RuntimeError(f"{accession} not found in {xml_path}")


def parse_hits(iteration):
    """Return dict[chrom] -> list of HSP dicts; dict[chrom] -> chr_length."""
    per_chr = defaultdict(list)
    chr_len = {}
    for hit in iteration.findall(".//Hit"):
        hit_def = hit.findtext("Hit_def", "") or ""
        m = re.search(r"chromosome\s+([0-9XYM]+|MT)", hit_def)
        if not m:
            continue
        chrom = m.group(1)
        chr_len[chrom] = int(hit.findtext("Hit_len", "0"))
        for hsp in hit.findall(".//Hsp"):
            qf = int(hsp.findtext("Hsp_query-from", "0"))
            qt = int(hsp.findtext("Hsp_query-to", "0"))
            hf = int(hsp.findtext("Hsp_hit-from", "0"))
            ht = int(hsp.findtext("Hsp_hit-to", "0"))
            alen = int(hsp.findtext("Hsp_align-len", "0") or 1)
            nid  = int(hsp.findtext("Hsp_identity", "0"))
            pid  = 100 * nid / alen if alen else 0.0
            per_chr[chrom].append({
                "qa": min(qf, qt), "qb": max(qf, qt),
                "ha": min(hf, ht), "hb": max(hf, ht),
                "alen": alen, "pid": pid,
            })
    return per_chr, chr_len


def make_figure(per_chr, chr_len, qlen, mask_window=(8262, 8658)):
    # Grayscale palette (publication-safe; no color)
    C_BACKBONE = "#cfcfcf"   # chromosome / viral backbone fill
    C_MASK     = "#202020"   # mask candidate region (near-black)
    C_TICK     = "#000000"   # HSP tick marks
    C_HSPCNT   = "#333333"   # HSP-count label

    fig = plt.figure(figsize=(12.0, 9.0))
    gs  = fig.add_gridspec(
        2, 1, height_ratios=[1.1, 4.0], hspace=0.22, left=0.09,
        right=0.97, top=0.95, bottom=0.07,
    )

    # ── Top: viral genome ───────────────────────────────────────────
    ax_v = fig.add_subplot(gs[0])
    ax_v.barh(0, qlen, height=0.55, color=C_BACKBONE,
              edgecolor="black", linewidth=0.7)
    ms, me = mask_window
    ax_v.barh(0, me - ms, left=ms, height=0.55,
              color=C_MASK, edgecolor="black", linewidth=0.7)
    # annotation placed above the bar, arrow points to the mask region
    ax_v.annotate(
        f"Mask candidate region in viral ref\n"
        f"(cross-mapping window  {ms:,}–{me:,} bp,  397 bp, 0.24 %)",
        xy=((ms + me) / 2, 0.30),
        xytext=(qlen * 0.18, 2.3),
        ha="left", va="center", fontsize=13,
        arrowprops=dict(arrowstyle="->", lw=1.0, color="black"),
    )
    ax_v.set_xlim(0, qlen)
    ax_v.set_ylim(-1.0, 3.4)
    ax_v.text(-0.075, 1.02, "A", transform=ax_v.transAxes,
              ha="left", va="bottom", fontsize=20, fontweight="bold")
    ax_v.set_yticks([])
    ax_v.tick_params(axis="x", labelsize=11)
    ax_v.set_xlabel(f"{QUERY_ACC}  ({QUERY_NAME})  position (bp)", fontsize=13)
    for s in ("top", "right"):
        ax_v.spines[s].set_visible(False)

    # ── Bottom: human chromosomes ───────────────────────────────────
    ax_h = fig.add_subplot(gs[1])
    max_len = max(chr_len.values())
    y_step  = 1.0
    bar_h   = 0.45

    for i, c in enumerate(CHR_ORDER):
        if c not in chr_len:
            continue
        y = -i * y_step
        clen = chr_len[c]
        # chromosome backbone
        ax_h.barh(y, clen, height=bar_h, color=C_BACKBONE,
                  edgecolor="black", linewidth=0.5, zorder=1)
        # hit ticks (rasterised so 280k vertical lines don't bloat EPS)
        rows = per_chr.get(c, [])
        if rows:
            xs = [(r["ha"] + r["hb"]) / 2 for r in rows]
            ax_h.vlines(
                xs, ymin=y - bar_h / 2, ymax=y + bar_h / 2,
                color=C_TICK, linewidth=0.18, alpha=0.55,
                rasterized=True, zorder=2,
            )
        # label + HSP count
        ax_h.text(-max_len * 0.012, y, f"chr{c}",
                  ha="right", va="center", fontsize=11)
        ax_h.text(clen + max_len * 0.005, y, f"{len(rows):>6,} HSPs",
                  ha="left", va="center", fontsize=10, color=C_HSPCNT)

    ax_h.set_xlim(-max_len * 0.05, max_len * 1.12)
    ax_h.set_ylim(-(len(CHR_ORDER) - 0.4) * y_step, 0.6)
    ax_h.text(-0.075, 1.02, "B", transform=ax_h.transAxes,
              ha="left", va="bottom", fontsize=20, fontweight="bold")
    ax_h.set_yticks([])
    ax_h.tick_params(axis="x", labelsize=11)
    ax_h.set_xlabel("Human T2T-CHM13v2.0 chromosome position (bp, true relative scale)",
                    fontsize=13)
    for s in ("top", "right", "left"):
        ax_h.spines[s].set_visible(False)

    # legend (grayscale)
    leg = [
        mpatches.Patch(facecolor=C_BACKBONE, edgecolor="black", linewidth=0.5,
                       label="Human chromosome"),
        mpatches.Patch(facecolor=C_TICK, edgecolor="black", linewidth=0.5,
                       label="megaBLAST HSP (~281 bp, 85–93 % id)"),
    ]
    ax_h.legend(handles=leg, loc="lower right", fontsize=11,
                frameon=False, ncol=1)

    return fig


def main():
    print(f"[1/3] extracting <Iteration> for {QUERY_ACC} …")
    iteration, qlen = extract_iteration(XML_PATH, QUERY_ACC)
    print(f"      query length = {qlen:,} bp")

    print(f"[2/3] parsing HSPs …")
    per_chr, chr_len = parse_hits(iteration)
    total = sum(len(v) for v in per_chr.values())
    print(f"      {len(per_chr)} chromosomes, {total:,} HSPs")

    print(f"[3/3] rendering figure …")
    fig = make_figure(per_chr, chr_len, qlen)
    save_figure(fig, FIG_NAME, outdir=OUT_DIR, close=True)


if __name__ == "__main__":
    main()
