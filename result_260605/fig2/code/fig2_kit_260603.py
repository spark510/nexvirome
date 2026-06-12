#!/usr/bin/env python3
"""
Fig2 data (FINALIZED 2026-06-03): NexVirome (new mask, method B, breadth 0.01)
vs 4 external tools (Kraken2 / Ganon / Metabuli / Phanta), each run TWICE:
  - same-DB : external tool given OUR viral-only DB  (external/kit_*)      -> algorithm comparison
  - native  : external tool on its OWN integrated DB (external_native/kit_*) -> real-world

For every (tool x db_mode) we score, at SPECIES and GENUS level vs GROUND_TRUTH:
  TP / FP / FN / precision / recall / F1   +   VAE (TPM, genome-length normalized)
and a per-GT-species TPM relative-abundance table (every tool gives read counts).

VAE = mean |observed_TPM_frac - 1/6| over the 6 GT taxa.
TPM: per matched taxid RPK = reads/(genome_len/1000); summed within a GT group;
normalized over the 6 GT groups. genus VAE == species VAE here (1 GT species = 1 genus).

NexVirome = our pipeline on HQ parquet (method B = best-hit + breadth>=0.01,
mask_v3_full.bed). It has only one "db_mode" (our DB) -> reported once as same-DB.

Outputs -> result_260605/fig2/:
  fig2_kit_scores.csv              tool x db_mode x (per-sample + MEAN): sp/gn TP/FP/F1 + VAE
  fig2_kit_rel_abundance_TPM.csv   tool x db_mode x sample x GT-species: reads, raw%, TPM%
Run: conda run -n shotgun_virome python scripts/benchmark/fig2_kit_260603.py
"""
from __future__ import annotations
import os, sys, glob
import numpy as np
import pandas as pd

NX = "/home/share/programs/nexvirome"
sys.path.insert(0, f"{NX}/scripts"); sys.path.insert(0, f"{NX}/notebooks")
sys.path.insert(0, f"{NX}/scripts/benchmark")
import three_strategy_breadth_lca as TS
from virome_classifier.alignment.filters.filter import MaskingFilter
from virome_classifier.taxonomy import TaxonomyDB
from benchmark_utils import (GROUND_TRUTH, parse_kreport, kreport_to_species_counts,
                             parse_metabuli_report, metabuli_to_species_counts)

DB    = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
MASK  = f"{NX}/result_260605/mask/mask_v3_full.bed"
GLEN  = f"{NX}/paper/figures/Fig5_extra/tables/species_genome_length.csv"
EXT   = f"{NX}/resources/db_20260525/external"          # same-DB
NAT   = f"{NX}/resources/db_20260525/external_native"   # native
CACHE = "/tmp/hq_cache"
SAMPLES = ["MagNA_1", "MagNA_2", "Qiagen_1", "Qiagen_2"]
# GOLDEN RULE (result_260605/GOLDEN_RULE.md): method B = best-hit + breadth>=CUT
# + per-taxon read floor n>=MIN_TAXON_READS (applied at SPECIES taxid level,
# BEFORE any genus roll-up), rel-abund OFF.
CUT = 0.01            # unmasked breadth (strict; was 0.005)
MIN_TAXON_READS = 3   # per-taxon read floor (species-level, pre-rollup)
EXPECTED = 1.0 / 6
OUT = f"{NX}/result_260605/fig2"
os.makedirs(OUT, exist_ok=True)

GT_GROUPS = [(g["name"], set(g["taxids"]), g["genus"]) for g in GROUND_TRUTH]
GT_ORDER = [g[0] for g in GT_GROUPS]
GT_TAXIDS = set().union(*[s for _, s, _ in GT_GROUPS])
GT_GENERA = {gid for _, _, gid in GT_GROUPS}

# ---- shared taxonomy / genome length ----
TS.DB = DB; TS._TAX_PATH = DB; TS._TLEN = TS._load_tlen()
TAX = TaxonomyDB.from_sqlite(DB)
GLEN_MAP = dict(pd.read_csv(GLEN).itertuples(index=False))

# Representative genome length per GT GROUP. External tools (same-DB) classify GT
# spike-ins under NEWER species taxids (e.g. CMV 3050295, Adeno 3241426) that are
# absent from species_genome_length.csv (built on the older taxids NexVirome uses).
# Since every member taxid of a GT group is the SAME virus, use the group's known
# member length(s) as the per-group genome length, so RPK works for any matched
# taxid regardless of which taxid version the tool reported.
GT_GROUP_LEN = {}
for _name, _tx, _ in GT_GROUPS:
    cands = [GLEN_MAP[t] for t in _tx if t in GLEN_MAP]
    GT_GROUP_LEN[_name] = max(cands) if cands else 0


def _glen_for(taxid, group_name):
    """Genome length for a matched taxid: its own length, else the GT group's."""
    return GLEN_MAP.get(taxid) or GT_GROUP_LEN.get(group_name, 0)


# ============ scoring (species + genus, TPM-VAE) ============
def score_counts(counts):
    """counts = {taxid: reads}. Return per-sample dict with sp_/gn_ TP/FP/F1 + VAE."""
    if not counts:
        return None
    # species TP/FP/FN
    sp_tp = sum(1 for _, tx, _ in GT_GROUPS if sum(counts.get(t, 0) for t in tx) >= 1)
    sp_fp = sum(1 for t, c in counts.items() if t not in GT_TAXIDS and c >= 1)
    sp_fn = len(GT_GROUPS) - sp_tp
    # genus roll-up
    gcounts = TS._genus_counts(counts, TAX)
    gn_tp = sum(1 for gid in GT_GENERA if gcounts.get(gid, 0) >= 1)
    gn_fp = sum(1 for g, c in gcounts.items() if g not in GT_GENERA and c >= 1)
    gn_fn = len(GT_GENERA) - gn_tp

    def prf(tp, fp, fn):
        p = tp / (tp + fp) if tp + fp else 0.0
        r = tp / (tp + fn) if tp + fn else 0.0
        f = 2 * p * r / (p + r) if p + r else 0.0
        return p, r, f
    sp_p, sp_r, sp_f = prf(sp_tp, sp_fp, sp_fn)
    gn_p, gn_r, gn_f = prf(gn_tp, gn_fp, gn_fn)

    # TPM VAE over GT groups (per-taxid RPK with group-fallback length, summed per group)
    rpk = {}
    for name, tx, _ in GT_GROUPS:
        tot = 0.0
        for t in tx:
            r = counts.get(t, 0)
            if r <= 0:
                continue
            L = _glen_for(t, name)
            if L > 0:
                tot += r / (L / 1000)
        rpk[name] = tot
    ksum = sum(rpk.values())
    vae = np.mean([abs(rpk[n] / ksum - EXPECTED) for n in GT_ORDER]) if ksum else float("nan")
    return dict(sp_TP=sp_tp, sp_FP=sp_fp, sp_FN=sp_fn, sp_F1=round(sp_f, 4),
                gn_TP=gn_tp, gn_FP=gn_fp, gn_FN=gn_fn, gn_F1=round(gn_f, 4),
                VAE=round(float(vae), 4))


def tpm_rel_abundance(counts):
    """Per-GT-species reads, raw%, TPM% (GT-normalized)."""
    rows = []
    raw = {n: sum(counts.get(t, 0) for t in tx) for n, tx, _ in GT_GROUPS}
    rpk = {}
    for n, tx, _ in GT_GROUPS:
        tot = 0.0
        for t in tx:
            r = counts.get(t, 0)
            if r > 0:
                L = _glen_for(t, n)
                if L > 0:
                    tot += r / (L / 1000)
        rpk[n] = tot
    rsum, ksum = sum(raw.values()), sum(rpk.values())
    for n in GT_ORDER:
        rows.append(dict(gt_taxon=n, reads=raw[n],
                         raw_pct=round(raw[n] / rsum * 100, 2) if rsum else 0.0,
                         TPM_pct=round(rpk[n] / ksum * 100, 2) if ksum else 0.0))
    return rows


def _taxid_len(taxid):
    """Genome length for any taxid (FP or TP): own length, else the GT-group's
    representative length if it belongs to one, else 0 (skipped from TPM)."""
    if taxid in GLEN_MAP:
        return GLEN_MAP[taxid]
    for name, tx, _ in GT_GROUPS:
        if taxid in tx:
            return GT_GROUP_LEN.get(name, 0)
    return 0


def tpm_with_fp(counts, level="species"):
    """TPM composition over the 6 GT taxa PLUS one lumped 'FP' category, summing
    to 100%. level='species' uses taxid counts; level='genus' rolls counts up to
    genus first (GT genera vs non-GT genera). FP TPM = sum of per-taxon RPK of all
    non-GT taxa (each by its own genome length; unknown-length taxa skipped)."""
    if level == "genus":
        c = TS._genus_counts(counts, TAX)
        gt_keysets = [(name, {gid}) for name, _, gid in GT_GROUPS]
        gt_all = GT_GENERA
    else:
        c = counts
        gt_keysets = [(name, tx) for name, tx, _ in GT_GROUPS]
        gt_all = GT_TAXIDS

    rpk = {}
    for name, keys in gt_keysets:
        tot = 0.0
        for k in keys:
            r = c.get(k, 0)
            if r > 0:
                L = _taxid_len(k) if level == "species" else (GT_GROUP_LEN.get(name, 0) or _taxid_len(k))
                if L > 0:
                    tot += r / (L / 1000)
        rpk[name] = tot
    # FP: every taxon not in GT
    fp = 0.0
    for k, r in c.items():
        if k in gt_all or r <= 0:
            continue
        L = _taxid_len(k)
        if L > 0:
            fp += r / (L / 1000)
    rpk["FP"] = fp
    tot = sum(rpk.values())
    return {k: (v / tot * 100 if tot else 0.0) for k, v in rpk.items()}


# ====== NexVirome counts (GOLDEN method B: best-hit + breadth>=CUT + n>=3) ======
def nexvirome_counts():
    bed = pd.read_csv(MASK, sep="\t", header=None, usecols=[0, 1, 2], names=["target", "start", "end"])
    mf = MaskingFilter.from_dataframe(bed)
    out = {}
    for s in SAMPLES:
        df = TS._add_tlen(pd.read_parquet(f"{CACHE}/{s}.parquet"))
        best = TS._best_hit(df)
        refs = TS._refs_at(TS._breadth_by_ref(best, mf), CUT)
        c = TS._besthit_counts(best[best["target"].isin(refs)])
        out[s] = {int(t): int(n) for t, n in c.items() if n >= MIN_TAXON_READS}  # GOLDEN: n>=3
    return out


# ============ external tool counts (same-DB / native) ============
def ganon_tre_counts(path, viral_only_db=False):
    """Parse a ganon .tre. For a NATIVE (integrated) DB keep only species whose
    lineage passes the Viruses node (10239). For a SAME-DB run the DB is already
    viral-only and its lineage has NO 10239 node (root '1||phylum|...'), so keep
    all species."""
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["rank", "taxid", "lineage", "name", "uniq", "shared", "children", "cum", "cumpct"])
    sp = df[df["rank"] == "species"]
    if not viral_only_db:
        sp = sp[sp["lineage"].astype(str).str.contains(r"(?:^|\|)10239(?:\||$)")]
    return {int(r.taxid): int(r.cum) for r in sp.itertuples() if int(r.cum) > 0}


def kraken_counts(path):
    return kreport_to_species_counts(parse_kreport(path), filter_virus=True) if os.path.exists(path) else None


def metabuli_counts(folder_glob):
    cand = glob.glob(folder_glob)
    if not cand:
        return None
    return metabuli_to_species_counts(parse_metabuli_report(cand[0]), filter_virus=True)


# Metabuli-Bracken species counts (abundance re-estimated; new_est_reads).
# Metabuli docs: refseq_standard == Kraken2 PlusPF, so the PlusPF Bracken matrix
# applies. Built by run_metabuli_bracken.py -> result_260605/fig2/metabuli_bracken/.
MB_DIR = f"{NX}/result_260605/fig2/metabuli_bracken"


K2B_DIR = f"{NX}/result_260605/fig2/kraken2_bracken"


def _bracken_viral_counts(path):
    """Read a Bracken table (.bracken), keep viral species (under 10239),
    return {taxid: new_est_reads}."""
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, sep="\t")
    out = {}
    for _, r in df.iterrows():
        tid = int(r["taxonomy_id"]); reads = int(r["new_est_reads"])
        if reads <= 0:
            continue
        try:
            if TAX.get_taxid_at_rank(tid, "superkingdom") == 10239:
                out[tid] = reads
        except Exception:
            pass
    return out or None


def metabuli_bracken_counts(sample):
    return _bracken_viral_counts(f"{MB_DIR}/{sample}.bracken")


def kraken2_bracken_counts(sample, db_mode):
    return _bracken_viral_counts(f"{K2B_DIR}/{db_mode}/{sample}.bracken")


def external_counts(tool, db_mode, sample):
    base = EXT if db_mode == "same-DB" else NAT
    if tool == "Ganon":
        return ganon_tre_counts(f"{base}/kit_ganon/{sample}.tre",
                                viral_only_db=(db_mode == "same-DB"))
    if tool == "Kraken2":
        return kraken_counts(f"{base}/kit_kraken2/{sample}.kreport")
    if tool == "Metabuli":
        return metabuli_counts(f"{base}/kit_metabuli/{sample}/*report.tsv") or \
               metabuli_counts(f"{base}/kit_metabuli/{sample}*report.tsv")
    if tool == "Phanta":
        for pat in [f"{base}/kit_phanta/results/classification/{sample}.krak.report_bracken_species.filtered",
                    f"{base}/kit_phanta/results/classification/{sample}.krak.report.filtered"]:
            if os.path.exists(pat):
                return kreport_to_species_counts(parse_kreport(pat), filter_virus=True)
        return None
    return None


# ============ main ============
def main():
    score_rows, ab_rows, comp_rows = [], [], []

    def add(tool, db_mode, counts_by_sample, abund_by_sample=None):
        """Detection (TP/FP/F1) uses counts_by_sample. Relative abundance + VAE use
        abund_by_sample when given (e.g. Metabuli native -> Bracken-re-estimated
        counts); otherwise the same detection counts. VAE in the score row is taken
        from the abundance counts so it reflects the abundance estimator used."""
        persamp = []
        for s in SAMPLES:
            c = counts_by_sample.get(s)
            if not c:
                continue
            ac = (abund_by_sample or {}).get(s) or c
            sc = score_counts(c)                      # detection from classifier
            if sc:
                sc["VAE"] = round(score_counts(ac)["VAE"], 4)  # abundance VAE
                persamp.append(dict(tool=tool, db_mode=db_mode, sample=s, **sc))
            for r in tpm_rel_abundance(ac):           # GT-only rel-abundance
                ab_rows.append(dict(tool=tool, db_mode=db_mode, sample=s, **r))
            # FP-inclusive composition (6 GT taxa + lumped FP, sum=100), sp & gn
            for level in ("species", "genus"):
                comp = tpm_with_fp(ac, level=level)
                for cat, pct in comp.items():
                    comp_rows.append(dict(tool=tool, db_mode=db_mode, sample=s,
                                          level=level, category=cat, TPM_pct=round(pct, 3)))
        if not persamp:
            score_rows.append(dict(tool=tool, db_mode=db_mode, sample="MEAN", status="no_output"))
            return
        score_rows.extend(persamp)
        dfp = pd.DataFrame(persamp)
        mean = dict(tool=tool, db_mode=db_mode, sample="MEAN")
        for k in ["sp_TP", "sp_FP", "sp_FN", "sp_F1", "gn_TP", "gn_FP", "gn_FN", "gn_F1", "VAE"]:
            mean[k] = round(dfp[k].mean(), 4)
        score_rows.append(mean)

    # NexVirome (our DB; reported as same-DB)
    print("scoring NexVirome (method B, breadth 0.01, mask_v3_full)...", flush=True)
    add("NexVirome", "same-DB", nexvirome_counts())

    # external tools x {same-DB, native}
    for tool in ["Ganon", "Kraken2", "Metabuli", "Phanta"]:
        for mode in ["same-DB", "native"]:
            cbs = {s: external_counts(tool, mode, s) for s in SAMPLES}
            # Abundance re-estimation with Bracken where supported. Detection
            # (TP/FP/F1) always stays from each tool's own classification.
            #   Kraken2 (both modes): Bracken is Kraken2's native abundance estimator.
            #   Metabuli native: PlusPF Bracken matrix (Metabuli docs).
            # Ganon/Phanta provide their own abundance (ganon reads; phanta bracken-built),
            # so they use raw counts directly.
            abund = None
            if tool == "Kraken2":
                abund = {s: kraken2_bracken_counts(s, mode) for s in SAMPLES}
            elif tool == "Metabuli" and mode == "native":
                abund = {s: metabuli_bracken_counts(s) for s in SAMPLES}
            if abund is not None:
                nb = sum(1 for v in abund.values() if v)
                print(f"  {tool:9s} {mode:8s}: Bracken abundance {nb}/4", flush=True)
            add(tool, mode, cbs, abund_by_sample=abund)
            n_ok = sum(1 for v in cbs.values() if v)
            print(f"  {tool:9s} {mode:8s}: {n_ok}/4 samples parsed", flush=True)

    sc = pd.DataFrame(score_rows)
    ab = pd.DataFrame(ab_rows)
    comp = pd.DataFrame(comp_rows)
    sc.to_csv(f"{OUT}/fig2_kit_scores.csv", index=False)
    ab.to_csv(f"{OUT}/fig2_kit_rel_abundance_TPM.csv", index=False)
    comp.to_csv(f"{OUT}/fig2_kit_composition_with_FP.csv", index=False)

    pd.set_option("display.width", 220)
    mean = sc[sc["sample"] == "MEAN"].copy()
    print("\n===== Fig2 KIT MEAN — species TP/FP/F1/VAE | genus TP/FP/F1/VAE =====")
    cols = ["tool", "db_mode", "sp_TP", "sp_FP", "sp_F1", "VAE", "gn_TP", "gn_FP", "gn_F1"]
    print(mean.reindex(columns=cols).to_string(index=False))
    print(f"\n-> {OUT}/fig2_kit_scores.csv , fig2_kit_rel_abundance_TPM.csv")


if __name__ == "__main__":
    main()
