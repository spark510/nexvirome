#!/usr/bin/env python3
"""
Reproducible SupFig3 input: the no/pre/post masking regime comparison on the KIT
mock, under the GOLDEN_RULE config (method B = best-hit + unmasked breadth ≥ 0.01
+ per-taxon read floor n ≥ 3, rel-abundance OFF; mask_v3_full; DB
tax_seq_v20260526_MSL41). Same parameters as Fig2/3/4 so the masking-regime
result is consistent with the rest of result_260605.

Every number here is produced by running the actual classifier, so the figure is
regenerable.

The three masking regimes (the only variable is HOW the mask is applied):
  No mask            : unmasked alignment + EMPTY bed         (no masking at all)
  Pre-mask (DB-level): alignment vs N-substituted DB FASTA    (mask baked into the
                       reference; classifier given an empty bed so it adds nothing)
  Post-mask (full)   : unmasked alignment + mask_v3_full coord bed (mask applied
                       only to coverage-breadth scoring — NexVirome's approach)
Single classifier = method B (the GOLDEN production method). The earlier
LCA/Coverage/EM mode sweep is dropped: masking effect is the variable of interest,
not the classification mode.

Outputs (overwrite, consistent run):
  result_260605/SupFig3/tables/masking_regimes_methodB.csv
  benchmark_runs/masking_pre_vs_post.csv  (kept for back-compat)

Run:  conda run -n shotgun_virome python result_260605/SupFig3/code/run_masking_regimes.py
"""
from __future__ import annotations
import sys, os, subprocess, tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
sys.path.insert(0, "/home/share/programs/nexvirome/notebooks")
sys.path.insert(0, "/home/share/programs/nexvirome/scripts")
import numpy as np, pandas as pd
from benchmark_utils import (GROUND_TRUTH, evaluate_sample, parse_kreport,
                             kreport_to_species_counts)

NX = "/home/share/programs/nexvirome"
DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
MASK = f"{NX}/result_260605/mask/mask_v3_full.bed"          # GOLDEN mask
UNMASKED = "/home/share/programs/vshot/result_kit/newdb_unmasked_mmseqs2"  # unmasked aln
PREMASKED = "/home/share/programs/vshot/result_kit/newdb_masked_mmseqs2"   # N-subst DB aln
OUT = "/tmp/masking_regimes_golden"
SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]

# empty BED => MaskingFilter with zero masked coordinates (no masking effect)
EMPTY_BED = "/tmp/masking_regimes_empty.bed"

# GOLDEN_RULE method B: best-hit read assignment + unmasked breadth ≥ 0.01
# + per-taxon read floor n ≥ 3, rel-abundance OFF. (= Fig2/3/4)
METHOD_B = ["--min-identity", "0.85", "--min-length", "60",
            "--read-assign", "best_hit",
            "--min-unmasked-coverage", "0.01",
            "--min-read-count", "3",
            "--min-rel-abundance", "0"]
MIN_READS = 3   # scoring floor (matches the per-taxon n≥3 above)
MODES = {
    "method B": ("lca", list(METHOD_B)),   # lca mode + read-assign best_hit = method B
}
# regime -> (alignment dir, mask file)
REGIMES = {
    "No mask":             (UNMASKED, EMPTY_BED),
    "Pre-mask (DB-level)": (PREMASKED, EMPTY_BED),
    "Post-mask (full)":    (UNMASKED, MASK),
}


def _job(args):
    """One classify run (regime x mode x sample). Module-level => picklable."""
    regime, modelabel, mode, extra, s = args
    rdir, mask = REGIMES[regime]
    out = f"{OUT}/{regime.replace(' ', '_').replace('(', '').replace(')', '')}/{modelabel}/{s}"
    cmd = [sys.executable, "-m", "virome_classifier.cli.classify", "--mode", mode,
           "--r1", f"{rdir}/{s}_R1.result", "--r2", f"{rdir}/{s}_R2.result",
           "--taxonomy", DB, "--mask", mask, "--output", out, "--sample", s] + extra
    subprocess.run(cmd, capture_output=True, cwd=f"{NX}/scripts")
    kp = f"{out}/{s}.kreport"
    return regime, modelabel, s, (kp if os.path.exists(kp) else None)


def main():
    open(EMPTY_BED, "w").close()                       # 0-byte BED = no masked coords
    jobs = [(rg, ml, mode, extra, s)
            for rg in REGIMES
            for ml, (mode, extra) in MODES.items()
            for s in SAMPLES]
    print(f"dispatching {len(jobs)} classify jobs (regime x mode x sample)...", flush=True)
    kps = {}
    with ProcessPoolExecutor(max_workers=16) as ex:
        for fu in as_completed([ex.submit(_job, j) for j in jobs]):
            rg, ml, s, kp = fu.result()
            kps[(rg, ml, s)] = kp

    rows = []
    for rg in REGIMES:
        for ml in MODES:
            F, FP = [], []
            for s in SAMPLES:
                kp = kps.get((rg, ml, s))
                if not kp:
                    continue
                c = kreport_to_species_counts(parse_kreport(kp))
                ev = evaluate_sample(c, GROUND_TRUTH, min_reads=MIN_READS)
                F.append(ev["f1"]); FP.append(ev["FP"])
            if F:
                rows.append(dict(masking=rg, mode=ml,
                                 F1=round(np.mean(F), 3), FP=round(np.mean(FP), 1)))
                print(f"  {rg:20s} {ml:9s} F1={np.mean(F):.3f} FP={np.mean(FP):.1f}", flush=True)

    df = pd.DataFrame(rows)
    out_tbl = f"{NX}/result_260605/SupFig3/tables"
    os.makedirs(out_tbl, exist_ok=True)
    os.makedirs(f"{NX}/benchmark_runs", exist_ok=True)
    df.to_csv(f"{out_tbl}/masking_regimes_methodB.csv", index=False)
    df.to_csv(f"{NX}/benchmark_runs/masking_pre_vs_post.csv", index=False)  # back-compat
    print("\n=== masking regimes (GOLDEN method B, one consistent run) ===")
    print(df.to_string(index=False))
    print(f"\nsaved -> {out_tbl}/masking_regimes_methodB.csv")


if __name__ == "__main__":
    main()
