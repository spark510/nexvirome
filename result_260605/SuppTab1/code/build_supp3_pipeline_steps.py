#!/usr/bin/env python3
"""
Supplementary Table S3 — read counts at each pipeline step, per sample / cohort.

Five steps (all in READS; R1 and R2 are processed as independent queries by the
classifier, so they are counted separately throughout):

  raw        raw reads                         fastp before_filtering.total_reads
  qc_pass    reads passing fastp QC            fastp after_filtering.total_reads
  nonhost    host-removed (viral candidates)   bowtie2 unmapped pairs * 2
  viral_aln  aligned to the viral DB (mmseqs2) distinct R1+R2 queries with >=1 hit
  confident  confidently classified            GOLDEN method B survivors (see below)

`raw / qc_pass / nonhost / viral_aln` are taken verbatim from the existing
pipeline-step table (paper/tables/pipeline_qc_steps_per_sample.csv): those steps
are independent of the classifier and `viral_aln` (= a read produced any DB hit)
does not change when the downstream confidence rule changes.

`confident` is RE-COMPUTED here under the locked GOLDEN_RULE method B (it used to
be the raw HitQuality gate; the confidence definition changed):
    best-hit  ->  reference unmasked-breadth >= 0.01  ->  per-species read floor n>=3
applied with the SAME helpers the Fig2/Fig3/Fig4 detection uses (TS._best_hit /
_breadth_by_ref / _refs_at). `confident` is then counted as the number of
DISTINCT READS (the `query` column) assigned to a kept species, to match the
read-level unit of `viral_aln` and keep the table monotonic. (Note: the parquet
`read_id` is per-alignment-segment, ~2x the read count, so it must NOT be used as
the read unit here — `query` is the read identity, as in the upstream mmseqs
`.result` distinct-query count of `viral_aln`.)
Input = HQ parquet /tmp/hq_cache/{sample}.parquet (HitQualityFilter pre-applied:
fident>=0.85, alnlen>=60, qcov>=0.5, e<=1e-3).

Outputs (result_260605/SuppTab1/tables/):
  suppTbl1_pipeline_steps_per_sample.csv
  suppTbl1_pipeline_steps_cohort_summary.csv

  conda run -n shotgun_virome python result_260605/SuppTab1/code/build_supp3_pipeline_steps.py
"""
from __future__ import annotations
import os, sys
import pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/scripts")
sys.path.insert(0, f"{NX}/scripts/benchmark")
sys.path.insert(0, f"{NX}/result_260605")
import three_strategy_breadth_lca as TS
from virome_classifier.alignment.filters.filter import MaskingFilter
from golden_rule import keep_samples   # drops EXCLUDE_SAMPLES (vir17)

# --- GOLDEN_RULE (result_260605/GOLDEN_RULE.md): method B confidence rule ---
DB    = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
MASK  = f"{NX}/result_260605/mask/mask_v3_full.bed"
CACHE = "/tmp/hq_cache"
CUT = 0.01            # unmasked breadth (strict)
MIN_TAXON_READS = 3   # per-species read floor, applied BEFORE any roll-up
TS.DB = DB; TS._TAX_PATH = DB; TS._TLEN = TS._load_tlen()

# upstream steps come from the existing classifier-independent table
PREV = f"{NX}/paper/tables/pipeline_qc_steps_per_sample.csv"
OUT  = f"{NX}/result_260605/SuppTab1/tables"


def confident_reads(sample, mf):
    """Distinct reads (query) surviving GOLDEN method B
    (best-hit -> breadth>=CUT -> per-species n>=3). Returns None if no HQ parquet.

    The per-species floor n>=3 is evaluated on the best-hit read count per taxid
    (same as nexvirome_counts); reads of taxids that clear breadth AND the floor
    are then counted as distinct `query` IDs to stay on the read unit."""
    p = f"{CACHE}/{sample}.parquet"
    if not os.path.exists(p):
        return None
    df = TS._add_tlen(pd.read_parquet(p))
    best = TS._best_hit(df)                                   # one row per read_id
    refs = TS._refs_at(TS._breadth_by_ref(best, mf), CUT)     # breadth-passing refs
    sub = best[best["target"].isin(refs)].copy()
    if sub.empty:
        return 0
    # per-species floor on best-hit read count (taxid units, pre-rollup, n>=3)
    per_tax = sub["taxid"].value_counts()
    kept_tax = per_tax.index[per_tax.values >= MIN_TAXON_READS]
    sub = sub[sub["taxid"].isin(kept_tax)]
    # report as DISTINCT READS (query), matching viral_aln's read unit
    return int(sub["query"].nunique())


def main():
    prev = pd.read_csv(PREV)
    # result_260605: drop EXCLUDE_SAMPLES (vir17) before anything else, so the
    # per-sample table and the recomputed cohort aggregates exclude it everywhere.
    prev = prev[prev["sample"].isin(keep_samples(prev["sample"]))].copy()
    # columns we keep verbatim (classifier-independent / DB-hit step)
    keep = prev[["cohort", "sample",
                 "raw_reads", "post_fastp_reads",
                 "host_filtered_reads", "viral_aligned_reads"]].copy()
    keep = keep.rename(columns={
        "raw_reads": "raw", "post_fastp_reads": "qc_pass",
        "host_filtered_reads": "nonhost", "viral_aligned_reads": "viral_aln"})

    bed = pd.read_csv(MASK, sep="\t", header=None, usecols=[0, 1, 2],
                      names=["target", "start", "end"])
    mf = MaskingFilter.from_dataframe(bed)

    conf = []
    for _, r in keep.iterrows():
        conf.append(confident_reads(r["sample"], mf))
        print(f"  {r['cohort']:4s} {r['sample']:12s} confident={conf[-1]}", flush=True)
    keep["confident"] = conf

    cols = ["cohort", "sample", "raw", "qc_pass", "nonhost", "viral_aln", "confident"]
    df = keep[cols]
    df.to_csv(f"{OUT}/suppTbl1_pipeline_steps_per_sample.csv", index=False)

    # cohort summary: mean count per step + mean % of raw
    s = df.copy()
    for c in ["qc_pass", "nonhost", "viral_aln", "confident"]:
        s[f"{c}_pct"] = s[c] / s["raw"] * 100
    summ = (s.groupby("cohort")
            .agg(n=("sample", "nunique"),
                 raw_mean=("raw", "mean"),
                 qc_pass_pct=("qc_pass_pct", "mean"),
                 nonhost_pct=("nonhost_pct", "mean"),
                 viral_aln_pct=("viral_aln_pct", "mean"),
                 confident_pct=("confident_pct", "mean"))
            .reset_index())
    summ["raw_mean"] = summ["raw_mean"].round(0).astype("int64")
    for c in ["qc_pass_pct", "nonhost_pct", "viral_aln_pct", "confident_pct"]:
        summ[c] = summ[c].round(4)
    summ.to_csv(f"{OUT}/suppTbl1_pipeline_steps_cohort_summary.csv", index=False)

    print("\n=== per-sample (head) ===")
    print(df.head(8).to_string(index=False))
    print("\n=== cohort summary (mean, % of raw) ===")
    print(summ.to_string(index=False))
    print(f"\nsaved -> {OUT}/suppTbl1_pipeline_steps_per_sample.csv + _cohort_summary.csv")


if __name__ == "__main__":
    main()
