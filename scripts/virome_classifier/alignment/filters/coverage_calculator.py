"""
Coverage calculation engine.

This module contains all coverage calculation logic, separated from
filtering and data loading concerns.
"""
from __future__ import annotations

from typing import Dict, Tuple
import numpy as np
import pandas as pd

from .models import CoverageStats, CoverageMetrics
from .mask_loader import MaskedRegion


class CoverageCalculator:
    """
    Calculator for genome coverage statistics.

    This class encapsulates all coverage calculation logic:
    - Merged interval coverage (handles overlaps)
    - Breadth and depth calculation
    - Masked vs unmasked hit classification

    All calculations are vectorized using NumPy for performance.
    """

    def __init__(self, overlap_threshold: float = 0.5):
        """
        Initialize calculator.

        Args:
            overlap_threshold: Minimum overlap ratio to classify hit as "masked"
                              (default=0.5 means 50%+ overlap required)
        """
        if not 0.0 <= overlap_threshold <= 1.0:
            raise ValueError(
                f"overlap_threshold must be in [0.0, 1.0], got {overlap_threshold}"
            )
        self.overlap_threshold = overlap_threshold

    def calculate_merged_coverage(
        self,
        intervals: np.ndarray,
        target_length: int
    ) -> float:
        """
        Calculate coverage ratio from potentially overlapping intervals.

        This merges overlapping intervals before calculating coverage to avoid
        double-counting overlapping regions.

        Args:
            intervals: Array of [[start1, end1], [start2, end2], ...]
                      (0-based, half-open coordinates)
            target_length: Total target sequence length in bp

        Returns:
            Coverage ratio (0.0 - 1.0), capped at 1.0
        """
        if len(intervals) == 0:
            return 0.0

        # Vectorised covered-length via a difference array (replaces a per-interval
        # Python merge loop that dominated the masking filter, ~4.5s on a 4M-hit
        # sample). Covered bp = number of positions with coverage >= 1, which is
        # exactly what merging overlapping [start,end) intervals and summing their
        # lengths yields. NOTE: the previous loop did NOT clamp to target_length —
        # it summed merged lengths as-is and only capped the final ratio at 1.0 —
        # so to stay byte-identical we size the array to the max end (no clamp) and
        # cap the ratio the same way.
        iv = np.asarray(intervals)
        s = iv[:, 0].astype(np.int64)
        e = iv[:, 1].astype(np.int64)
        valid = e > s
        if not valid.any():
            return 0.0
        s = s[valid]; e = e[valid]
        # shift so the smallest start is index 0 (handles negative starts safely)
        lo = int(s.min())
        size = int(e.max()) - lo + 1
        diff = np.zeros(size + 1, dtype=np.int64)
        np.add.at(diff, s - lo, 1)
        np.add.at(diff, e - lo, -1)
        total_covered = int((diff[:-1].cumsum() > 0).sum())

        # Cap at 1.0 (defensive: handles edge cases where coordinates exceed target length)
        return min(total_covered / target_length, 1.0)

    def calculate_breadth_and_depth(
        self,
        intervals: np.ndarray,
        target_length: int
    ) -> Tuple[int, float]:
        """
        Calculate breadth of coverage (bp covered) and average depth.

        Breadth = number of unique positions covered (≥1x coverage)
        Depth = average coverage depth across covered positions

        Args:
            intervals: Array of [[start1, end1], [start2, end2], ...]
            target_length: Total target sequence length in bp

        Returns:
            Tuple of (breadth_bp: int, avg_depth: float)
        """
        if len(intervals) == 0:
            return 0, 0.0

        # Vectorised difference-array (replaces a per-interval Python loop that was
        # the dominant masking-filter cost, ~23s on a 4M-hit sample): mark +1 at
        # each clamped start and -1 at each clamped end, then cumsum to get the
        # exact same per-position coverage as `coverage[s:e] += 1` would.
        iv = np.asarray(intervals)
        starts = np.clip(iv[:, 0].astype(np.int64), 0, target_length - 1)
        ends = np.clip(iv[:, 1].astype(np.int64), 0, target_length)
        valid = starts < ends
        starts, ends = starts[valid], ends[valid]
        if starts.size == 0:
            return 0, 0.0
        diff = np.zeros(target_length + 1, dtype=np.int64)
        np.add.at(diff, starts, 1)
        np.add.at(diff, ends, -1)
        coverage = np.cumsum(diff[:-1])

        breadth_bp = int(np.count_nonzero(coverage))
        avg_depth = float(coverage[coverage > 0].mean()) if breadth_bp > 0 else 0.0
        return breadth_bp, avg_depth

    def _calculate_metrics(
        self,
        intervals: np.ndarray,
        target_length: int
    ) -> CoverageMetrics:
        """
        Calculate all coverage metrics for a set of intervals.

        Args:
            intervals: Array of intervals
            target_length: Target length

        Returns:
            CoverageMetrics object
        """
        if len(intervals) == 0:
            return CoverageMetrics.empty()

        coverage_ratio = self.calculate_merged_coverage(intervals, target_length)
        breadth_bp, avg_depth = self.calculate_breadth_and_depth(intervals, target_length)

        return CoverageMetrics(
            hit_count=len(intervals),
            breadth_bp=breadth_bp,
            coverage_ratio=coverage_ratio,
            avg_depth=avg_depth
        )

    def identify_masked_hits(
        self,
        hit_starts: np.ndarray,
        hit_ends: np.ndarray,
        masked_region: MaskedRegion
    ) -> np.ndarray:
        """
        Identify which hits overlap with masked regions (vectorized).

        A hit is considered "masked" if overlap_threshold or more of its length
        overlaps with any masked region.

        Args:
            hit_starts: Array of hit start positions (n_hits,)
            hit_ends: Array of hit end positions (n_hits,)
            masked_region: MaskedRegion object for this target

        Returns:
            Boolean array indicating which hits are masked (True = masked)
        """
        n_hits = len(hit_starts)

        if n_hits == 0 or masked_region.n_regions == 0:
            return np.zeros(n_hits, dtype=bool)

        hit_lengths = hit_ends - hit_starts  # 0-based half-open (n_hits,)

        # Broadcast for vectorized overlap calculation
        hit_starts_2d = hit_starts[:, None]  # (n_hits, 1)
        hit_ends_2d = hit_ends[:, None]      # (n_hits, 1)

        # Calculate overlap length for each hit-mask pair (vectorized, half-open intervals)
        overlap_starts = np.maximum(hit_starts_2d, masked_region.starts)  # (n_hits, n_masks)
        overlap_ends = np.minimum(hit_ends_2d, masked_region.ends)        # (n_hits, n_masks)
        overlap_lengths = np.maximum(0, overlap_ends - overlap_starts)    # (n_hits, n_masks)

        # Calculate overlap ratio for each hit-mask pair
        overlap_ratios = overlap_lengths / hit_lengths[:, None]  # (n_hits, n_masks)

        # A hit is masked if ANY mask overlaps >= threshold
        is_masked = (overlap_ratios >= self.overlap_threshold).any(axis=1)  # (n_hits,)

        return is_masked

    def calculate_stats(
        self,
        target_hits: pd.DataFrame,
        mask_dict: Dict[str, MaskedRegion]
    ) -> CoverageStats:
        """
        Calculate comprehensive coverage statistics for a single target.

        This is the main entry point for coverage calculation.

        Args:
            target_hits: DataFrame of hits for a single target
                        Required columns: ['target', 'tstart', 'tend', 'tlen']
                        Optional columns: ['taxid']
            mask_dict: Dictionary mapping target name to MaskedRegion

        Returns:
            CoverageStats object with comprehensive statistics

        Raises:
            ValueError: If required columns are missing
        """
        required_cols = {"target", "tstart", "tend", "tlen"}
        if not required_cols.issubset(target_hits.columns):
            raise ValueError(
                f"DataFrame must have columns: {required_cols}, "
                f"got: {set(target_hits.columns)}"
            )

        if target_hits.empty:
            raise ValueError("Cannot calculate stats for empty DataFrame")

        # Extract target info
        target_name = target_hits.iloc[0]["target"]
        target_length = int(target_hits.iloc[0]["tlen"])
        taxid = target_hits.iloc[0].get("taxid", None)

        # Get all hit intervals
        hit_intervals = target_hits[["tstart", "tend"]].values

        # === TOTAL METRICS (all hits) ===
        total_metrics = self._calculate_metrics(hit_intervals, target_length)

        # === MASKED vs UNMASKED CLASSIFICATION ===
        masked_region = mask_dict.get(target_name)
        has_mask = masked_region is not None

        if not has_mask:
            # No mask data → all hits are unmasked
            masked_metrics = CoverageMetrics.empty()
            unmasked_metrics = total_metrics
        else:
            # Classify hits as masked or unmasked
            hit_starts = target_hits["tstart"].values
            hit_ends = target_hits["tend"].values
            is_masked = self.identify_masked_hits(hit_starts, hit_ends, masked_region)

            # Calculate metrics separately
            masked_intervals = hit_intervals[is_masked]
            unmasked_intervals = hit_intervals[~is_masked]

            masked_metrics = self._calculate_metrics(masked_intervals, target_length)
            unmasked_metrics = self._calculate_metrics(unmasked_intervals, target_length)

        return CoverageStats(
            target_name=target_name,
            target_length=target_length,
            taxid=taxid,
            has_mask=has_mask,
            total=total_metrics,
            masked=masked_metrics,
            unmasked=unmasked_metrics
        )

    def calculate_stats_for_all_targets(
        self,
        df: pd.DataFrame,
        mask_dict: Dict[str, MaskedRegion]
    ) -> pd.DataFrame:
        """
        Calculate coverage statistics for all targets in DataFrame.

        Args:
            df: DataFrame with hits for multiple targets
            mask_dict: Dictionary of masked regions

        Returns:
            DataFrame with one row per target containing all statistics
        """
        if df.empty:
            return pd.DataFrame()

        stats_list = []
        for target, group in df.groupby("target"):
            try:
                stats = self.calculate_stats(group, mask_dict)
                stats_list.append(stats.to_series())
            except Exception as e:
                # Log error but continue with other targets
                print(f"Warning: Failed to calculate stats for {target}: {e}")
                continue

        if not stats_list:
            return pd.DataFrame()

        return pd.DataFrame(stats_list)

    def per_hit_masked_labels(
        self,
        df: pd.DataFrame,
        mask_dict: Dict[str, MaskedRegion]
    ) -> "pd.Series":
        """
        Per-hit masked/unmasked label for ALL targets in one pass.

        A hit is 'masked' when >= overlap_threshold of its aligned length lies
        inside any mask region of its target (same rule as identify_masked_hits,
        which this delegates to per target). Hits on targets without a mask are
        unmasked. Returns a boolean Series aligned to df.index.

        This is the single source of truth for per-hit masking — callers that
        need labels without the breadth gate (e.g. all-masked-read dropping)
        use this instead of re-implementing the overlap math.
        """
        labels = pd.Series(False, index=df.index)
        if df.empty or not mask_dict:
            return labels
        starts = df["tstart"].to_numpy()
        ends = df["tend"].to_numpy()
        for target, idx in df.groupby("target").groups.items():
            region = mask_dict.get(target)
            if region is None or region.n_regions == 0:
                continue
            pos = df.index.get_indexer(idx)  # positional, robust to non-reset index
            is_masked = self.identify_masked_hits(starts[pos], ends[pos], region)
            labels.loc[idx] = is_masked
        return labels
