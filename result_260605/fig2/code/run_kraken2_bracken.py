#!/usr/bin/env python3
"""
Kraken2 abundance re-estimation with Bracken (Bracken's native use case).

Kraken2 .kreport is already standard Kraken format -> fed directly to
est_abundance.py (no conversion). Two db_modes, each with its own Bracken matrix:
  same-DB : reads classified vs OUR viral DB; Bracken DB = external/kraken2/database150mers.kmer_distrib
  native  : PlusPF DB;                        Bracken DB = pluspf_20250714/database150mers.kmer_distrib
reads ~150bp -> 150mers ; push to species (-l S).

Outputs -> result_260605/fig2/kraken2_bracken/{db_mode}/{sample}.bracken
Run: conda run -n shotgun_virome python scripts/benchmark/run_kraken2_bracken.py
"""
from __future__ import annotations
import os, subprocess

NX = "/home/share/programs/nexvirome"
SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
CFG = {
    "same-DB": dict(
        rep=f"{NX}/resources/db_20260525/external/kit_kraken2",
        kmer=f"{NX}/resources/db_20260525/external/kraken2/database150mers.kmer_distrib"),
    "native": dict(
        rep=f"{NX}/resources/db_20260525/external_native/kit_kraken2",
        kmer="/home/share/kraken_db/pluspf_20250714/database150mers.kmer_distrib"),
}
OUT = f"{NX}/result_260605/fig2/kraken2_bracken"


def main():
    for mode, cfg in CFG.items():
        odir = f"{OUT}/{mode}"
        os.makedirs(odir, exist_ok=True)
        if not os.path.exists(cfg["kmer"]):
            print(f"  {mode}: kmer_distrib MISSING ({cfg['kmer']}) -> skip"); continue
        for s in SAMPLES:
            kr = f"{cfg['rep']}/{s}.kreport"
            if not os.path.exists(kr):
                print(f"  {mode}/{s}: kreport missing, skip"); continue
            bk = f"{odir}/{s}.bracken"
            bkrep = f"{odir}/{s}.bracken_report"
            r = subprocess.run(["est_abundance.py", "-i", kr, "-k", cfg["kmer"],
                                "-o", bk, "-l", "S", "--out-report", bkrep, "-t", "1"],
                               capture_output=True, text=True)
            ok = os.path.exists(bk)
            print(f"  {mode}/{s}: bracken {'OK' if ok else 'FAIL'}")
            if not ok:
                print("    ", r.stderr.strip().splitlines()[-1] if r.stderr.strip() else "(no msg)")
    print(f"\n-> {OUT}/{{same-DB,native}}/*.bracken")


if __name__ == "__main__":
    main()
