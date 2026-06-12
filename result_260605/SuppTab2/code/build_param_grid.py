#!/usr/bin/env python3
"""
SuppTab2 — 2-parameter GRID sweep for the read-level filters, so the *combined*
stability of (identity, qcov, length) can be shown as 2D heatmaps instead of the
one-variable-at-a-time line plots.

For each of the three parameter PAIRS, both axes are swept jointly while the
THIRD parameter is held at its production default; the KIT mock (ground truth) is
classified at every grid cell and precision = TP/(TP+FP) and FP are recorded.
(Only KIT is gridded — it is the cohort with ground truth, the 4 KIT .result
files are tiny, and precision is what defines the "stable region".)

Reuses build_param_sensitivity's classify job (_job / DEF / DB / MASK / SCRATCH)
so the kreport cache is shared with the 1D sweep.

Outputs (result_260605/SuppTab2/tables/):
  param_grid_KIT.csv   pair, x_param, x_value, y_param, y_value, sample,
                       reads, n_species, TP, FP, F1, recall   (per-sample)
  param_grid_KIT_summary.csv  pair, x_param, x_value, y_param, y_value,
                       precision_mean, FP_mean, recall_mean, F1_mean   (cell mean)

Run: /usr/local/bin/miniconda3/envs/shotgun_virome/bin/python \
       result_260605/SuppTab2/code/build_param_grid.py
"""
from __future__ import annotations
import os, sys, subprocess
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, f"{NX}/scripts"); sys.path.insert(0, f"{NX}/scripts/benchmark")
import build_param_sensitivity as B   # reuse DB/MASK/DEF/SCRATCH/RAW_DIR/_list_samples
from benchmark_utils import GROUND_TRUTH, evaluate_sample, parse_kreport, kreport_to_species_counts

OUT = f"{NX}/result_260605/SuppTab2/tables"
SCRATCH = f"{B.SCRATCH}/grid"          # separate grid cache
MIN_READS = B.MIN_READS
N_WORKERS = B.N_WORKERS

# axis values per parameter (reuse the 1D sweep grids)
VALS = {
    "identity": [0.70, 0.75, 0.80, 0.85, 0.90],
    "qcov":     [0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90],
    "length":   [20, 35, 50, 65, 80, 95, 110, 130, 150],
}
PAIRS = [("identity", "qcov"), ("identity", "length"), ("qcov", "length")]
COHORT = "KIT"


def _grid_job(args):
    """Classify one KIT sample at a (xparam=xval, yparam=yval) grid cell; the
    third parameter stays at DEF. Returns a record dict."""
    xp, xv, yp, yv, s = args
    p = dict(B.DEF); p[xp] = xv; p[yp] = yv
    rdir = B.RAW_DIR[COHORT]
    tag = f"{xp}{xv}_{yp}{yv}"
    out = f"{SCRATCH}/{tag}/{s}"
    os.makedirs(out, exist_ok=True)
    kp = f"{out}/{s}.kreport"
    if not os.path.exists(kp):
        cmd = [sys.executable, "-m", "virome_classifier.cli.classify", "--mode", "lca",
               "--read-assign", "best_hit",
               "--r1", f"{rdir}/{s}_R1.result", "--r2", f"{rdir}/{s}_R2.result",
               "--taxonomy", B.DB, "--mask", B.MASK,
               "--no-breadth-gate",
               "--min-identity", str(p["identity"]),
               "--min-length", str(int(p["length"])),
               "--min-query-coverage", str(p["qcov"]),
               "--min-read-count", "0", "--min-rel-abundance", "0",
               "--output", out, "--sample", s]
        subprocess.run(cmd, capture_output=True, cwd=f"{NX}/scripts")
    if not os.path.exists(kp):
        return None
    counts = kreport_to_species_counts(parse_kreport(kp))
    ev = evaluate_sample(counts, GROUND_TRUTH, min_reads=MIN_READS)
    tp, fp = ev["TP"], ev["FP"]
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    return dict(x_param=xp, x_value=xv, y_param=yp, y_value=yv, sample=s,
                reads=int(sum(counts.values())), n_species=len(counts),
                TP=tp, FP=fp, F1=round(ev["f1"], 4),
                recall=round(ev["recall"], 4), precision=round(prec, 4))


def main():
    os.makedirs(OUT, exist_ok=True)
    samples = B._list_samples(COHORT)             # the 4 KIT mock libraries
    jobs = []
    for xp, yp in PAIRS:
        for xv in VALS[xp]:
            for yv in VALS[yp]:
                for s in samples:
                    jobs.append((xp, xv, yp, yv, s))
    print(f"grid jobs: {len(jobs)} ({len(PAIRS)} pairs x KIT {len(samples)} samples)",
          flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        for r in ex.map(_grid_job, jobs):
            if r is not None:
                # attach which pair this cell belongs to
                r["pair"] = f"{r['x_param']}x{r['y_param']}"
                rows.append(r)

    df = pd.DataFrame(rows)
    cols = ["pair", "x_param", "x_value", "y_param", "y_value", "sample",
            "reads", "n_species", "TP", "FP", "F1", "recall", "precision"]
    df = df[cols]
    df.to_csv(f"{OUT}/param_grid_KIT.csv", index=False)

    summ = (df.groupby(["pair", "x_param", "x_value", "y_param", "y_value"])
              .agg(precision_mean=("precision", "mean"),
                   FP_mean=("FP", "mean"),
                   recall_mean=("recall", "mean"),
                   F1_mean=("F1", "mean"))
              .reset_index())
    summ = summ.round(4)
    summ.to_csv(f"{OUT}/param_grid_KIT_summary.csv", index=False)

    print(f"-> {OUT}/param_grid_KIT.csv ({len(df)} rows)")
    print(f"-> {OUT}/param_grid_KIT_summary.csv ({len(summ)} cells)")
    print("\n=== precision range per pair ===")
    for pair, g in summ.groupby("pair"):
        print(f"  {pair}: precision {g.precision_mean.min():.3f}–{g.precision_mean.max():.3f}, "
              f"FP {g.FP_mean.min():.1f}–{g.FP_mean.max():.1f}")


if __name__ == "__main__":
    main()
