#!/usr/bin/env python3
"""
SupFig2 — per-chromosome HSP statistics behind the BeAn 58058 (NC_032111.1)
cross-mapping schematic. Re-uses the XML parser of
`fig_bean58058_crossmap_schematic.py` and writes the numbers the figure shows
(HSP count, alignment length, identity, strand balance) as a table, so the figure
has an auditable data backing.

Outputs (result_260605/SupFig2/tables/):
  supfig2_bean58058_hsp_per_chromosome.csv   chrom, chr_len_bp, n_hsp, median_alen,
                                             median_pid, plus_strand, minus_strand
  supfig2_bean58058_hsp_overall.csv          one-row global summary

  conda run -n shotgun_virome python result_260605/SupFig2/code/build_supfig2_hsp_table.py
"""
from __future__ import annotations
import os, re, sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from statistics import median

NX = "/home/share/programs/nexvirome"
XML_PATH = f"{NX}/resources/db/ncbi/ncbi_human_genome/refseq_viral_to_hgt2t.xml"
QUERY_ACC = "NC_032111.1"
OUT = f"{NX}/result_260605/SupFig2/tables"
CHR_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y"]


def extract_iteration(xml_path, accession):
    """Stream the multi-document XML; return (iteration_root, query_len) for the
    single <Iteration> whose query-def starts with `accession`."""
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
    """dict[chrom] -> list of HSP dicts (alen, pid, strand); dict[chrom] -> chr_len."""
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
            hf = int(hsp.findtext("Hsp_hit-from", "0"))
            ht = int(hsp.findtext("Hsp_hit-to", "0"))
            qf = int(hsp.findtext("Hsp_query-from", "0"))   # viral (BeAn) side
            qt = int(hsp.findtext("Hsp_query-to", "0"))
            alen = int(hsp.findtext("Hsp_align-len", "0") or 1)
            nid = int(hsp.findtext("Hsp_identity", "0"))
            pid = 100 * nid / alen if alen else 0.0
            # hit strand: plus if hit-from <= hit-to, minus otherwise
            strand = "+" if hf <= ht else "-"
            per_chr[chrom].append({"alen": alen, "pid": pid, "strand": strand,
                                   "qf": qf, "qt": qt})
    return per_chr, chr_len


def main():
    import pandas as pd
    os.makedirs(OUT, exist_ok=True)
    print(f"[1/2] extracting <Iteration> for {QUERY_ACC} …")
    it, qlen = extract_iteration(XML_PATH, QUERY_ACC)
    per_chr, chr_len = parse_hits(it)

    rows = []
    for c in CHR_ORDER:
        rows_c = per_chr.get(c, [])
        if not rows_c:
            continue
        alens = [r["alen"] for r in rows_c]
        pids = [r["pid"] for r in rows_c]
        plus = sum(1 for r in rows_c if r["strand"] == "+")
        rows.append(dict(
            chrom=f"chr{c}",
            chr_len_bp=chr_len.get(c, 0),
            n_hsp=len(rows_c),
            median_alen=int(median(alens)),
            median_pid=round(median(pids), 1),
            plus_strand=plus,
            minus_strand=len(rows_c) - plus,
        ))
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/supfig2_bean58058_hsp_per_chromosome.csv", index=False)

    all_alen = [r["alen"] for v in per_chr.values() for r in v]
    all_pid = [r["pid"] for v in per_chr.values() for r in v]
    all_plus = sum(1 for v in per_chr.values() for r in v if r["strand"] == "+")
    n = len(all_alen)
    # viral-side (BeAn) cross-mapping window, COMPUTED from the HSP query
    # coordinates (min query-from to max query-to over all human HSPs).
    all_q = [p for v in per_chr.values() for r in v for p in (r["qf"], r["qt"])]
    win_start, win_end = min(all_q), max(all_q)
    win_bp = win_end - win_start + 1
    overall = pd.DataFrame([dict(
        query_acc=QUERY_ACC,
        query_name="BeAn 58058 virus",
        viral_genome_bp=qlen,
        mask_window=f"{win_start}-{win_end}",
        mask_window_bp=win_bp,
        mask_window_pct=round(win_bp / qlen * 100, 2),
        n_chromosomes=len(df),
        total_hsp=n,
        median_alen_bp=int(median(all_alen)),
        median_pid_pct=round(median(all_pid), 1),
        plus_strand=all_plus,
        minus_strand=n - all_plus,
        strand_ratio=round(all_plus / max(n - all_plus, 1), 3),
    )])
    overall.to_csv(f"{OUT}/supfig2_bean58058_hsp_overall.csv", index=False)

    print(f"[2/2] wrote tables -> {OUT}/")
    print("\n=== per-chromosome ===")
    print(df.to_string(index=False))
    print("\n=== overall ===")
    print(overall.T.to_string(header=False))


if __name__ == "__main__":
    main()
