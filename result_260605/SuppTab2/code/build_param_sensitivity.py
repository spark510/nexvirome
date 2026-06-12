#!/usr/bin/env python3
"""
Supplementary Table 2 — read-level parameter sensitivity WITHOUT the breadth gate.

Goal: show that the three read-level hit-quality filters (min identity, min query
coverage, min aligned length) each control detection ON THEIR OWN — even with the
breadth reference-filter switched off — and that the number of passing reads
swings sharply with each parameter. This isolates the read-level filter from the
breadth gate.

Design (one variable at a time; others held at default; NO downstream gate):
  - breadth OFF, per-taxon read floor OFF, rel-abundance OFF  <- pure read filter
  - best-hit assignment only; NO breadth, NO read floor, NO rel-abundance
  - sweep one of:
      identity ∈ {0.70, 0.75, 0.80, 0.85, 0.90}
      qcov     ∈ {0.30..0.90 step 0.10}
      length   ∈ {20,35,50,65,80,95,110,130,150}
  - defaults: identity 0.85, qcov 0.5, length 60.

Cohorts:
  - KIT  (4 samples, GROUND_TRUTH): species TP / FP / F1 / recall.
  - DNA  (5 representative samples) and RNA (5): NO ground truth, so we report the
    number of PASSING reads and detected species per sample — the read count
    "swings" with the parameter.
All cohorts also report passing-read count + detected-species count, so the
read-count sensitivity is shown for every cohort.

Input: raw (pre-HitQualityFilter) .result files, so values can be RELAXED below
the HQ-parquet floor (0.85/60/0.5):
  KIT  /home/share/programs/vshot/result_kit/newdb_unmasked_mmseqs2
  RNA  /home/share/programs/vshot/results_asthma/newdb_mmseqs2
  DNA  /home/share/data_processed/asthma_2024/newdb_mmseqs2

Outputs (result_260605/SuppTab2/tables/):
  param_sensitivity_per_sample.csv   cohort, param, value, sample, reads, n_species, TP, FP, F1, recall
  param_sensitivity_summary.csv      cohort, param, value, reads_mean, n_species_mean, (KIT: TP/FP/F1/recall mean)

Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260605/SuppTab2/code/build_param_sensitivity.py
"""
from __future__ import annotations
import os, sys, glob, subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/notebooks")
from benchmark_utils import GROUND_TRUTH, evaluate_sample, parse_kreport, kreport_to_species_counts
sys.path.insert(0, f"{NX}/result_260605")
from golden_rule import keep_samples   # drops EXCLUDE_SAMPLES (vir17)

DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
MASK = f"{NX}/result_260605/mask/mask_v3_full.bed"
OUT = f"{NX}/result_260605/SuppTab2/tables"

RAW_DIR = {
    "KIT": "/home/share/programs/vshot/result_kit/newdb_unmasked_mmseqs2",
    "RNA": "/home/share/programs/vshot/results_asthma/newdb_mmseqs2",
    "DNA": "/home/share/data_processed/asthma_2024/newdb_mmseqs2",
}
N_REAL = 5   # representative real-data samples per cohort

DEF = dict(identity=0.85, qcov=0.5, length=60)
# PURE read-level filter only: NO downstream gate at all — breadth OFF,
# per-taxon read floor OFF (n>=1), rel-abundance OFF. Detection = reads passing
# the swept identity/qcov/length filter (best-hit), nothing else.
MIN_READS = 1                      # no read floor (count any detected taxon)
SCRATCH = "/tmp/suptab2_paramsens_pure"   # separate cache from earlier runs
SWEEPS = {
    "identity": [0.70, 0.75, 0.80, 0.85, 0.90],
    "qcov":     [0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90],
    "length":   [20, 35, 50, 65, 80, 95, 110, 130, 150],
}
N_WORKERS = 96


def _list_samples(cohort):
    d = RAW_DIR[cohort]
    s = sorted({os.path.basename(f)[:-len("_R1.result")]
                for f in glob.glob(f"{d}/*_R1.result") if not f.endswith(".bak")})
    if cohort == "KIT":
        return ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
    # real data: pick N_REAL representatives by raw .result size (median band),
    # so we don't only sample the extremes. Exclude EXCLUDE_SAMPLES (vir17)
    # BEFORE selecting representatives so it can never be picked.
    s = keep_samples(s)
    sizes = [(x, os.path.getsize(f"{d}/{x}_R1.result")) for x in s]
    sizes.sort(key=lambda t: t[1])
    if len(sizes) <= N_REAL:
        return [x for x, _ in sizes]
    # evenly spaced indices across the size-sorted list
    idx = np.linspace(0, len(sizes) - 1, N_REAL).round().astype(int)
    return [sizes[i][0] for i in idx]


def _job(args):
    """One classify run (cohort, param, value, sample) with breadth gate OFF."""
    cohort, param, value, s = args
    p = dict(DEF); p[param] = value
    rdir = RAW_DIR[cohort]
    out = f"{SCRATCH}/{cohort}/{param}_{value}/{s}"
    os.makedirs(out, exist_ok=True)
    kp = f"{out}/{s}.kreport"
    if not os.path.exists(kp):
        cmd = [sys.executable, "-m", "virome_classifier.cli.classify", "--mode", "lca",
               "--read-assign", "best_hit",
               "--r1", f"{rdir}/{s}_R1.result", "--r2", f"{rdir}/{s}_R2.result",
               "--taxonomy", DB, "--mask", MASK,
               "--no-breadth-gate",             # breadth OFF
               "--min-identity", str(p["identity"]),
               "--min-length", str(int(p["length"])),
               "--min-query-coverage", str(p["qcov"]),
               "--min-read-count", "0",         # per-taxon read floor OFF
               "--min-rel-abundance", "0",      # rel-abundance OFF
               "--output", out, "--sample", s]
        subprocess.run(cmd, capture_output=True, cwd=f"{NX}/scripts")
    if not os.path.exists(kp):
        return cohort, param, value, s, None
    counts = kreport_to_species_counts(parse_kreport(kp))
    reads = int(sum(counts.values()))
    n_species = int(len(counts))     # any detected taxon (no read floor)
    rec = dict(cohort=cohort, param=param, value=value, sample=s,
               reads=reads, n_species=n_species, TP=None, FP=None, F1=None, recall=None)
    if cohort == "KIT":
        ev = evaluate_sample(counts, GROUND_TRUTH, min_reads=MIN_READS)
        rec.update(TP=ev["TP"], FP=ev["FP"], F1=round(ev["f1"], 4),
                   recall=round(ev["recall"], 4))
    return cohort, param, value, s, rec


def main():
    os.makedirs(OUT, exist_ok=True)
    cohort_samples = {c: _list_samples(c) for c in RAW_DIR}
    print("samples per cohort:",
          {c: len(v) for c, v in cohort_samples.items()},
          "\nreal-data picks:", {c: cohort_samples[c] for c in ("DNA", "RNA")})

    jobs = [(c, param, v, s)
            for c, samples in cohort_samples.items()
            for param, vals in SWEEPS.items()
            for v in vals
            for s in samples]
    print(f"dispatching {len(jobs)} classify jobs (breadth OFF; param sweep)...", flush=True)

    rows = []
    done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        for fu in as_completed([ex.submit(_job, j) for j in jobs]):
            c, param, value, s, rec = fu.result()
            done += 1
            if rec is not None:
                rows.append(rec)
            if done % 50 == 0:
                print(f"  progress {done}/{len(jobs)}", flush=True)
    per = pd.DataFrame(rows).sort_values(["cohort", "param", "value", "sample"])
    per.to_csv(f"{OUT}/param_sensitivity_per_sample.csv", index=False)

    agg = {"reads": ["mean", "std"], "n_species": ["mean", "std"]}
    summ = per.groupby(["cohort", "param", "value"]).agg(agg)
    summ.columns = [f"{a}_{b}" for a, b in summ.columns]
    # KIT scoring means
    kit = (per[per.cohort == "KIT"].groupby(["cohort", "param", "value"])
              .agg(TP_mean=("TP", "mean"), FP_mean=("FP", "mean"),
                   F1_mean=("F1", "mean"), recall_mean=("recall", "mean")))
    summ = summ.join(kit).round(2).reset_index()
    summ.to_csv(f"{OUT}/param_sensitivity_summary.csv", index=False)

    print(f"\n-> {OUT}/param_sensitivity_per_sample.csv ({len(per)} rows)")
    print(f"-> {OUT}/param_sensitivity_summary.csv")
    for c in ("KIT", "DNA", "RNA"):
        print(f"\n===== {c}: passing reads / species by parameter (breadth OFF) =====")
        for param in SWEEPS:
            sub = summ[(summ.cohort == c) & (summ.param == param)]
            cols = ["value", "reads_mean", "n_species_mean"]
            if c == "KIT":
                cols += ["TP_mean", "FP_mean", "F1_mean"]
            print(f"[{param}]")
            print(sub[cols].to_string(index=False))


if __name__ == "__main__":
    main()
