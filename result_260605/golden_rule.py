"""
GOLDEN RULE — single source of truth for result_260603 NexVirome parameters.

Every figure/table builder in result_260603 imports these constants instead of
hard-coding values, so the analysis can never drift. See GOLDEN_RULE.md for the
rationale. Method B (Strategy B) = best-hit + breadth>=0.01 + per-taxon n>=3,
rel-abundance OFF, on mask_v3_full / tax_seq_v20260526_MSL41 / HQ parquet.

Usage:
    import sys; sys.path.insert(0, "/home/share/programs/nexvirome/result_260603")
    from golden_rule import DB, MASK, BREADTH_CUT, MIN_TAXON_READS, apply_method_b

`apply_method_b(df, mask_filter, TS)` runs the full Strategy-B detection on one
sample's HQ-parquet DataFrame and returns {taxid: reads} already filtered by
breadth>=BREADTH_CUT and n>=MIN_TAXON_READS — the canonical detection call.
"""
from __future__ import annotations

NX = "/home/share/programs/nexvirome"

# ---- paths ----
# result_260605 = result_260603 re-run with vir17 EXCLUDED (see EXCLUDE_SAMPLES).
# The mask is self-contained inside this folder.
RESULT_DIR = f"{NX}/result_260605"
DB = f"{NX}/resources/db/custom/tax_seq_v20260526_MSL41.db"
MASK = f"{RESULT_DIR}/mask/mask_v3_full.bed"
GENOME_LENGTH_CSV = f"{NX}/paper/figures/Fig5_extra/tables/species_genome_length.csv"
HQ_CACHE = "/tmp/hq_cache"   # {sample}.parquet

# ---- Method B (Strategy B) detection parameters ----
BREADTH_CUT = 0.01        # unmasked breadth >= this for a reference to pass (strict; was 0.005)
MIN_TAXON_READS = 3       # per-taxon read floor: drop taxa with < this many reads
REL_ABUNDANCE_GATE = None  # OFF for method B
READ_ASSIGN = "best_hit"   # best-hit per read (no LCA)

# ---- excluded samples (result_260605) ----
# vir17 (DNA asthma) carries ~130k reads — an order of magnitude above any other
# sample — so its phage detections dominate every read-weighted aggregate and bias
# the cohort means. result_260605 drops it. Every sample-list builder must filter
# through `keep_samples()` so this exclusion propagates everywhere.
EXCLUDE_SAMPLES = {"vir17"}


def keep_samples(samples):
    """Filter an iterable of sample names, dropping EXCLUDE_SAMPLES (vir17)."""
    return [s for s in samples if s not in EXCLUDE_SAMPLES]

# ---- HitQualityFilter (HQ parquet build; LOCKED 2026-05-29) ----
HQ_MIN_IDENTITY = 0.85
HQ_MIN_ALN_LEN = 60
HQ_MIN_QCOV = 0.5
HQ_MAX_EVALUE = 1e-3

# ---- cohort manifests ----
DNA_GROUPS = f"{NX}/paper/figures/Fig3/source_data/dna_asthma_groups.csv"
RNA_GROUPS = f"{NX}/paper/figures/Fig3/source_data/rna_asthma_groups.csv"


def apply_method_b(df, mask_filter, TS, breadth_cut=BREADTH_CUT,
                   min_taxon_reads=MIN_TAXON_READS):
    """Canonical Strategy-B detection on one sample's HQ-parquet DataFrame.

    df          : HQ parquet for a sample (TS._add_tlen already applied, or apply here)
    mask_filter : MaskingFilter built from MASK
    TS          : the three_strategy_breadth_lca module (provides _best_hit etc.)

    Returns {taxid: reads}, filtered by unmasked breadth>=breadth_cut AND
    per-taxon reads>=min_taxon_reads. This is THE detection call all result_260603
    NexVirome panels should use (replaces the bare TS._besthit_counts that omitted
    the n>=3 floor)."""
    best = TS._best_hit(df)
    refs = TS._refs_at(TS._breadth_by_ref(best, mask_filter), breadth_cut)
    counts = TS._besthit_counts(best[best["target"].isin(refs)])
    return {int(t): int(n) for t, n in counts.items() if n >= min_taxon_reads}
