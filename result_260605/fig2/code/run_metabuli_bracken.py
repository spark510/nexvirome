#!/usr/bin/env python3
"""
Metabuli abundance re-estimation with Bracken (Metabuli-compatible per the
Metabuli docs: refseq_standard == Kraken2 PlusPF sequences/taxonomy, so the
Kraken2 PlusPF Bracken matrix applies).

Steps per KIT sample:
  1. Convert Metabuli native report -> Kraken-style report (rank full-name -> code,
     same column layout Bracken expects).
  2. Run Bracken est_abundance.py with PlusPF database150mers.kmer_distrib (reads ~150bp),
     pushing reads to species level (-l S).
  3. Bracken writes a _bracken species-abundance report; we keep it for Fig2.

Inputs : resources/db_20260525/external_native/kit_metabuli/{s}/{s}_report.tsv
DB     : /home/share/kraken_db/pluspf_20250714/database150mers.kmer_distrib
Outputs: result_260605/fig2/metabuli_bracken/{s}.kraken_report
         result_260605/fig2/metabuli_bracken/{s}.bracken          (Bracken table)
         result_260605/fig2/metabuli_bracken/{s}.bracken_report   (re-est kraken report)

Run: conda run -n shotgun_virome python scripts/benchmark/run_metabuli_bracken.py
"""
from __future__ import annotations
import os, subprocess, sys

NX = "/home/share/programs/nexvirome"
MET = f"{NX}/resources/db_20260525/external_native/kit_metabuli"
KMER = "/home/share/kraken_db/pluspf_20250714/database150mers.kmer_distrib"
OUT = f"{NX}/result_260605/fig2/metabuli_bracken"
SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
os.makedirs(OUT, exist_ok=True)

# Standard rank code -> Kraken tree depth (level_num). Bracken builds the tree by
# requiring each node's level_num == parent.level_num + 1, so we RE-NORMALIZE the
# indentation: assign a fixed depth to the 8 main ranks and slot every non-standard
# rank (subclass, clade, parvorder, ...) as depth(parent)+1 via a running counter.
# Metabuli's own deep/irregular indentation breaks stock Bracken; this rebuild lets
# stock est_abundance.py consume it while preserving the lineage order.
MAIN_DEPTH = {"D": 1, "K": 2, "P": 3, "C": 4, "O": 5, "F": 6, "G": 7, "S": 8}
RANK2CODE = {
    "superkingdom": "D", "domain": "D", "kingdom": "K", "phylum": "P",
    "class": "C", "order": "O", "family": "F", "genus": "G", "species": "S",
}


def convert(metabuli_report, kraken_out):
    """Metabuli report -> Kraken-style report with re-normalized indentation.
    Cols: clade_proportion, clade_count, taxon_count, rank, taxID, name
       -> pct, clade_reads, taxon_reads, rank_code, taxid, name(indent = 2*depth)."""
    n = 0
    prev_depth = 0
    with open(metabuli_report) as fi, open(kraken_out, "w") as fo:
        for ln in fi:
            if ln.startswith("#"):
                continue
            c = ln.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            pct, clade, taxon, rank, taxid = c[0], c[1], c[2], c[3].strip().lower(), c[4].strip()
            raw_name = c[5].strip()

            if taxid == "0":          # unclassified
                fo.write(f"{pct}\t{clade}\t{taxon}\tU\t0\tunclassified\n"); n += 1; continue
            if taxid == "1":          # root
                fo.write(f"{pct}\t{clade}\t{taxon}\tR\t1\troot\n"); prev_depth = 0; n += 1; continue

            if rank in RANK2CODE:
                code = RANK2CODE[rank]
                depth = MAIN_DEPTH[code]
            else:                      # non-standard rank -> child of previous node
                code = "-"
                depth = prev_depth + 1
            # Bracken's tree builder only allows level_num == parent+1 going DOWN;
            # an UP jump of +2 (e.g. a report that skips a main rank) breaks it.
            # Clamp any increase to exactly +1; decreases are fine (handled by the
            # parent walk). This preserves lineage order and the S leaf placement.
            if depth > prev_depth + 1:
                depth = prev_depth + 1
            indent = " " * (2 * depth)
            fo.write(f"{pct}\t{clade}\t{taxon}\t{code}\t{taxid}\t{indent}{raw_name}\n")
            prev_depth = depth
            n += 1
    return n


def main():
    for s in SAMPLES:
        rep = f"{MET}/{s}/{s}_report.tsv"
        if not os.path.exists(rep):
            print(f"  -- {s}: report missing, skip"); continue
        kr = f"{OUT}/{s}.kraken_report"
        nrows = convert(rep, kr)
        bk = f"{OUT}/{s}.bracken"
        bkrep = f"{OUT}/{s}.bracken_report"
        cmd = ["est_abundance.py", "-i", kr, "-k", KMER, "-o", bk,
               "-l", "S", "--out-report", bkrep, "-t", "1"]
        r = subprocess.run(cmd, capture_output=True, text=True)
        ok = os.path.exists(bk)
        print(f"  {s}: converted {nrows} rows -> bracken {'OK' if ok else 'FAIL'}")
        if not ok:
            print("    STDERR:", r.stderr.strip()[:300])
    print(f"\n-> {OUT}/  (*.bracken, *.bracken_report)")


if __name__ == "__main__":
    main()
