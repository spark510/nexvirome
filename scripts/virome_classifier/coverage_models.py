"""
Data models and pure helpers for the coverage-based classifier.

Extracted verbatim from coverage_based_classifier3.py to keep that module focused
on the classification pipeline (no behaviour change). Holds:
  - calculate_breadth_fast / compute_depth_entropy  (pure numeric helpers)
  - SegmentCoverage / StrainCoverage / TaxonCoverage (per-level coverage records)
  - CoverageThresholds                               (real/fake decision thresholds)
  - HitQualityFilter                                 (identity/length/qcov filter)
"""
from typing import Dict, Optional
from dataclasses import dataclass, field
import numpy as np
import pandas as pd

try:
    from numba import jit
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False


# =============================================================================
# Helper
# =============================================================================
def calculate_breadth_fast(intervals):
    if len(intervals) == 0:
        return 0

    covered_len = 0
    curr_s, curr_e = intervals[0]

    for s, e in intervals[1:]:
        if s <= curr_e:
            curr_e = max(curr_e, e)
        else:
            covered_len += (curr_e - curr_s)
            curr_s, curr_e = s, e

    covered_len += (curr_e - curr_s)
    return covered_len


if HAS_NUMBA:
    calculate_breadth_fast = jit(nopython=True)(calculate_breadth_fast)


def compute_depth_entropy(depth_array: np.ndarray) -> float:
    total = depth_array.sum()
    if total == 0:
        return 0.0

    p = depth_array / total
    p = p[p > 0]
    H = -np.sum(p * np.log(p))
    return float(H / np.log(len(depth_array)))


# =============================================================================
# Data Classes
# =============================================================================
@dataclass
class SegmentCoverage:
    segment_id: str
    length: int
    breadth_bp: int
    breadth_ratio: float
    avg_depth: float
    hit_count: int
    depth_entropy: Optional[float] = None
    # Read type breakdown
    unique_breadth_bp: int = 0
    multi_breadth_bp: int = 0
    masked_breadth_bp: int = 0


@dataclass
class StrainCoverage:
    strain_taxid: int
    virus_id: str
    total_genome_length: int
    segment_count: int
    expected_segment_count: int = 0
    segment_coverages: Dict[str, SegmentCoverage] = field(default_factory=dict)
    weighted_breadth: float = 0.0
    unmasked_weighted_breadth: float = 0.0
    segments_detected: int = 0
    unique_read_count: int = 0
    multi_mapping_read_count: int = 0
    total_read_count: int = 0
    masked_read_count: int = 0
    is_real: bool = False


@dataclass
class TaxonCoverage:
    taxon_taxid: int
    taxon_rank: str
    strains: Dict[int, StrainCoverage]
    avg_genome_length: float = 0.0
    weighted_breadth: float = 0.0
    unmasked_weighted_breadth: float = 0.0
    is_real: bool = False
    taxon_name: str = ""


@dataclass
class CoverageThresholds:
    min_weighted_breadth: float = 0.01
    min_segment_fraction: float = 0.1
    min_segment_breadth: float = 0.0
    min_segment_covered_bp: int = 200
    min_unmasked_coverage: float = 0.05  # ~15 spread reads on a median viral genome (44kb expected) — empirically as precise as 0.08 on KIT GT (recall 1.0, ref precision 78%) while recovering A1024O0001-like commensal phage signal that 0.08 truncates; see diag_breadth_vs_spread.csv / final_gate_comparison.csv
    require_unique_reads: bool = True
    min_depth_entropy: float = 0.20
    min_read_count: int = 10
    min_unmasked_weighted_breadth: float = 0.01
    # --- FP-leak fix (off by default => preserves current behavior) ---
    # NOTE: genus_competition was removed — the relative-abundance cut
    # (apply_fp_postfilter min_rel_abundance, default 0.05%) supersedes it.
    min_unique_fraction: float = 0.0


@dataclass
class HitQualityFilter:
    min_identity: float = 0.85
    min_aligned_length: int = 60
    min_query_coverage: float = 0.5

    def filter_df(self, df: pd.DataFrame) -> pd.DataFrame:
        fident = df['fident'] if 'fident' in df else df['pident'] / 100.0
        qcov = df.get('qcov', 1.0)

        mask = (
            (fident >= self.min_identity) &
            (df['alnlen'] >= self.min_aligned_length) &
            (qcov >= self.min_query_coverage)
        )
        return df[mask].copy()
