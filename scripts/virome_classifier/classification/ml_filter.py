"""
ML-based real/fake filter (CellTypist-style logistic regression), as a separate
mode that REFINES coverage-mode species calls instead of replacing them.

Why a TP/FP binary model (not read->species multi-class): the KIT mock has only 6
species, so a species-classifier would overfit and never generalize to asthma's
hundreds of species. Instead we learn, from the per-species alignment/coverage
FEATURES, the probability that a detected species is a true positive vs a
cross-mapping false positive. The label is species-level (in GT or not), the
features are species-agnostic -> generalizes to new datasets.

This is the learned analogue of the rule-based unique-fraction / genus-competition
gates: the model discovers the right combination from data.

Features (per detected species, from coverage_summary): weighted_breadth,
unmasked_weighted_breadth, unique_fraction, multi_fraction, masked_fraction,
log10(total_reads), log10(unique_reads), unique_breadth fraction, segments_detected,
best/avg entropy, strain_count, log10(avg_genome_length).

Usage (train):
  python -m virome_classifier.classification.ml_filter train \
     --summaries benchmark_runs/newdb_*/coverage_summary.tsv \
     --model resources/build/ml_filter.joblib
Usage (apply): handled in classify.py --mode ml_filter (loads model, scores species).
"""
from __future__ import annotations

import argparse
import glob
import sys

import numpy as np
import pandas as pd

# KIT ground-truth species taxids
GT = {10521, 28284, 130309, 3241426, 10359, 3050295, 11250, 11259, 208893,
      3049954, 11520, 2955465, 114101, 538123, 3422298, 351073, 64320, 3048459}

FEATURES = [
    "weighted_breadth", "unmasked_weighted_breadth",
    "unique_fraction", "multi_fraction", "masked_fraction",
    "log_total_reads", "log_unique_reads", "unique_breadth_frac",
    "segments_detected", "best_entropy", "avg_entropy",
    "strain_count", "log_genome_len",
    "both_frac",  # vectorized paired-end support (TP~0.84 vs FP~0.16 on KIT)
]


def build_features(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    tot = d["total_reads"].replace(0, 1)
    d["unique_fraction"] = d["unique_reads"] / tot
    d["multi_fraction"] = d["multi_reads"] / tot
    d["masked_fraction"] = d.get("masked_reads", 0) / tot
    d["log_total_reads"] = np.log10(d["total_reads"].clip(lower=1))
    d["log_unique_reads"] = np.log10(d["unique_reads"].clip(lower=1))
    ubp = d.get("unique_breadth_bp", 0)
    allbp = (d.get("unique_breadth_bp", 0) + d.get("multi_breadth_bp", 0)).replace(0, 1)
    d["unique_breadth_frac"] = ubp / allbp
    d["log_genome_len"] = np.log10(d.get("avg_genome_length", 1).clip(lower=1))
    for c in ("best_entropy", "avg_entropy", "segments_detected", "strain_count",
              "weighted_breadth", "unmasked_weighted_breadth", "both_frac"):
        if c not in d:
            d[c] = 0.0
    return d[FEATURES].fillna(0.0)


def train(args) -> int:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.pipeline import Pipeline
    from sklearn.model_selection import GroupKFold, cross_val_predict
    from sklearn.metrics import roc_auc_score, precision_recall_fscore_support
    import joblib

    files = []
    for pat in args.summaries:
        files.extend(glob.glob(pat))
    if not files:
        sys.exit("no coverage_summary files matched")
    dfs = []
    for f in files:
        d = pd.read_csv(f, sep="\t")
        d["__src"] = f
        dfs.append(d)
    df = pd.concat(dfs, ignore_index=True)
    df["label"] = df["taxon_taxid"].isin(GT).astype(int)
    X = build_features(df)
    y = df["label"].values
    groups = df["sample"] if "sample" in df else df["__src"]
    print(f"training rows: {len(df)} | TP {y.sum()} | FP {(~y.astype(bool)).sum()}")

    pipe = Pipeline([
        ("scale", StandardScaler()),
        ("lr", LogisticRegression(max_iter=1000, class_weight="balanced")),
    ])
    # grouped CV by sample to estimate honest generalization
    metrics = {"n_rows": int(len(df)), "n_tp": int(y.sum()),
               "n_fp": int((~y.astype(bool)).sum()), "n_groups": int(groups.nunique())}
    try:
        proba = cross_val_predict(pipe, X, y, groups=groups,
                                  cv=GroupKFold(n_splits=min(4, groups.nunique())),
                                  method="predict_proba")[:, 1]
        auc = roc_auc_score(y, proba)
        pred = (proba >= 0.5).astype(int)
        p, r, f1, _ = precision_recall_fscore_support(y, pred, average="binary", zero_division=0)
        metrics.update(cv_auc=round(float(auc), 3), cv_precision=round(float(p), 3),
                       cv_recall=round(float(r), 3), cv_f1=round(float(f1), 3))
        print(f"grouped-CV: AUC={auc:.3f}  P={p:.3f} R={r:.3f} F1={f1:.3f}")
    except Exception as e:
        print(f"CV skipped: {e}")

    pipe.fit(X, y)
    # report coefficients (which features matter)
    lr = pipe.named_steps["lr"]
    coef = sorted(zip(FEATURES, lr.coef_[0]), key=lambda t: -abs(t[1]))
    metrics["coefficients"] = {name: round(float(w), 3) for name, w in coef}
    # persist metrics inside the model and as a JSON sidecar so the paper can
    # quote the CV AUC/F1 and feature weights without re-running training.
    joblib.dump({"pipe": pipe, "features": FEATURES, "metrics": metrics}, args.model)
    import json
    sidecar = args.model.rsplit(".", 1)[0] + "_metrics.json"
    with open(sidecar, "w") as fh:
        json.dump(metrics, fh, indent=2)
    print("\ntop feature weights (TP-direction):")
    for name, w in coef[:8]:
        print(f"  {name:<24} {w:+.3f}")
    print(f"\nsaved model -> {args.model}")
    print(f"saved metrics -> {sidecar}")
    return 0


def score(summary_df: pd.DataFrame, model_path: str, threshold: float = 0.5) -> pd.DataFrame:
    """Apply a trained model to a coverage_summary df; returns df with tp_proba + ml_real."""
    import joblib
    obj = joblib.load(model_path)
    X = build_features(summary_df)
    proba = obj["pipe"].predict_proba(X)[:, 1]
    out = summary_df.copy()
    out["tp_proba"] = proba
    out["ml_real"] = proba >= threshold
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description="ML TP/FP filter for virome species calls.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    tr = sub.add_parser("train")
    tr.add_argument("--summaries", nargs="+", required=True)
    tr.add_argument("--model", required=True)
    args = ap.parse_args()
    if args.cmd == "train":
        return train(args)
    return 1


if __name__ == "__main__":
    sys.exit(main())
