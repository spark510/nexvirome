"""
Main masking filter class.

This module provides the high-level API for filtering viral genome alignments
based on masked region coverage.
"""
from __future__ import annotations

from typing import Dict, Optional
import pandas as pd
from pathlib import Path

from .models import FilterResult, CoverageStats
from .mask_loader import MaskedRegion, MaskLoader
from .coverage_calculator import CoverageCalculator


# Default thresholds
DEFAULT_MIN_UNMASKED_COV = 0.25  # 25%
DEFAULT_MIN_TOTAL_COV = 0.5      # 50%
DEFAULT_OVERLAP_THRESHOLD = 0.5  # 50%


class MaskingFilter:
    """
    High-level API for masking-aware filtering of alignment hits.

    This class encapsulates masked region data and provides filtering
    methods based on coverage statistics.

    Example:
        >>> # Create filter from BED file
        >>> filter = MaskingFilter.from_bed_file("masked_regions.bed")
        >>>
        >>> # Filter by unmasked coverage
        >>> result = filter.filter_by_unmasked_coverage(
        ...     hits_df,
        ...     min_coverage=0.25
        ... )
        >>>
        >>> print(result.summary())
        >>> passed_hits = result.passed
    """

    def __init__(
        self,
        mask_dict: Optional[Dict[str, MaskedRegion]] = None,
        overlap_threshold: float = DEFAULT_OVERLAP_THRESHOLD
    ):
        """
        Initialize masking filter.

        Args:
            mask_dict: Dictionary mapping target name to MaskedRegion
            overlap_threshold: Minimum overlap ratio to classify hit as masked
                             (default=0.5 means 50%+ overlap required)
        """
        self._mask_dict = mask_dict or {}
        self._calculator = CoverageCalculator(overlap_threshold=overlap_threshold)
        self._overlap_threshold = overlap_threshold

    @classmethod
    def from_bed_file(
        cls,
        bed_file: str,
        overlap_threshold: float = DEFAULT_OVERLAP_THRESHOLD
    ) -> 'MaskingFilter':
        """
        Factory method: Create filter from BED file.

        Args:
            bed_file: Path to BED file with masked regions
            overlap_threshold: Overlap threshold for masked hit classification

        Returns:
            MaskingFilter instance

        Raises:
            FileNotFoundError: If BED file doesn't exist
        """
        mask_dict = MaskLoader.from_bed_file(bed_file)
        print(f"📋 Loaded {len(mask_dict)} targets with masked regions from {bed_file}")
        return cls(mask_dict, overlap_threshold)

    @classmethod
    def from_dataframe(
        cls,
        df: pd.DataFrame,
        overlap_threshold: float = DEFAULT_OVERLAP_THRESHOLD
    ) -> 'MaskingFilter':
        """
        Factory method: Create filter from DataFrame.

        Args:
            df: DataFrame with columns ['target', 'start', 'end']
            overlap_threshold: Overlap threshold for masked hit classification

        Returns:
            MaskingFilter instance
        """
        mask_dict = MaskLoader.from_dataframe(df)
        return cls(mask_dict, overlap_threshold)

    @property
    def n_targets_with_masks(self) -> int:
        """Number of targets that have masked region data."""
        return len(self._mask_dict)

    @property
    def has_masks(self) -> bool:
        """Whether any masked region data is available."""
        return len(self._mask_dict) > 0

    def calculate_stats(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Calculate coverage statistics for all targets.

        Args:
            df: DataFrame with alignment hits

        Returns:
            DataFrame with per-target statistics
        """
        return self._calculator.calculate_stats_for_all_targets(df, self._mask_dict)

    def classify_hits(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Label each hit masked / unmasked WITHOUT applying the breadth gate.

        Returns a copy of df with a boolean `is_masked` column. This is the
        classification concern split out from filter_by_unmasked_coverage so
        masking can be inspected or acted on (e.g. dropping all-masked reads)
        independently of the breadth gate.
        """
        if df is None or df.empty:
            return df
        out = df.copy()
        out["is_masked"] = self._calculator.per_hit_masked_labels(out, self._mask_dict)
        return out

    def drop_all_masked_reads(
        self,
        df: pd.DataFrame,
        read_col: str = "query",
    ) -> pd.DataFrame:
        """
        Remove reads whose every hit is masked (all_masked reads = host-mimic /
        vector(UniVec) / rRNA / MAG contamination), and drop the masked hits of
        surviving reads so they don't reach LCA/count.

        Reads are grouped by `read_col` (default 'query'; R1/R2 already merged
        via normalize_headers), matching how LCA groups reads. Reads keeping at
        least one unmasked hit survive (their masked hits are dropped); reads
        with no unmasked hit are removed entirely.

        Validated 2026-05-30 (project_allmasked_drop_validated): no KIT TP loss,
        FP reduction. See also classify_hits for the labelling step.
        """
        if df is None or df.empty or not self.has_masks:
            return df
        labelled = self.classify_hits(df)
        # a read is all-masked iff min(is_masked) over its hits is True
        read_all_masked = labelled.groupby(read_col)["is_masked"].transform("min")
        keep = (~read_all_masked.astype(bool)) & (~labelled["is_masked"])
        return labelled[keep].drop(columns=["is_masked"])

    def _split_by_targets(
        self,
        df: pd.DataFrame,
        valid_targets: pd.Index
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Split DataFrame into passed and failed based on target list.

        Args:
            df: Original DataFrame
            valid_targets: Index of targets that passed

        Returns:
            Tuple of (passed_df, failed_df)
        """
        passed = df[df["target"].isin(valid_targets)].copy()
        failed = df[~df["target"].isin(valid_targets)].copy()
        return passed, failed

    def filter_by_unmasked_coverage(
        self,
        df: pd.DataFrame,
        min_coverage: float = DEFAULT_MIN_UNMASKED_COV
    ) -> FilterResult:
        """
        Filter targets based on unmasked region coverage.

        Targets pass if unmasked regions cover at least min_coverage of
        total target length. Coverage from masked regions is ignored.

        Args:
            df: DataFrame with alignment hits
                Required columns: ['target', 'tstart', 'tend', 'tlen']
            min_coverage: Minimum unmasked coverage ratio (default=0.25 = 25%)

        Returns:
            FilterResult object with passed/failed hits and statistics

        Example:
            >>> result = filter.filter_by_unmasked_coverage(hits_df, min_coverage=0.25)
            >>> print(f"Passed: {result.n_passed_targets} targets")
        """
        if df is None or df.empty:
            return FilterResult(
                passed=df or pd.DataFrame(),
                failed=pd.DataFrame(),
                stats=pd.DataFrame(),
                filter_name="unmasked_coverage"
            )

        if not self.has_masks:
            print("⚠️  No masking data - all targets pass")
            return FilterResult(
                passed=df,
                failed=pd.DataFrame(),
                stats=pd.DataFrame(),
                filter_name="unmasked_coverage"
            )

        print(f"🔍 Filtering by unmasked coverage (min={min_coverage:.1%})")

        # Calculate stats
        stats_df = self.calculate_stats(df)

        # Filter by unmasked coverage
        valid_mask = stats_df["unmasked_coverage_ratio"] >= min_coverage
        valid_targets = stats_df[valid_mask].index

        # Split hits
        passed, failed = self._split_by_targets(df, valid_targets)

        print(f"✅ Passed: {len(valid_targets)} targets ({len(passed)} hits)")
        print(f"❌ Failed: {len(stats_df) - len(valid_targets)} targets ({len(failed)} hits)")

        return FilterResult(
            passed=passed,
            failed=failed,
            stats=stats_df,
            filter_name="unmasked_coverage"
        )

    def filter_by_total_coverage(
        self,
        df: pd.DataFrame,
        min_coverage: float = DEFAULT_MIN_TOTAL_COV
    ) -> FilterResult:
        """
        Filter targets based on total coverage (including masked regions).

        Targets pass if total coverage (masked + unmasked) is at least min_coverage.
        This counts all hits regardless of whether they overlap masked regions.

        Args:
            df: DataFrame with alignment hits
            min_coverage: Minimum total coverage ratio (default=0.5 = 50%)

        Returns:
            FilterResult object with passed/failed hits and statistics
        """
        if df is None or df.empty:
            return FilterResult(
                passed=df or pd.DataFrame(),
                failed=pd.DataFrame(),
                stats=pd.DataFrame(),
                filter_name="total_coverage"
            )

        print(f"🔍 Filtering by total coverage (min={min_coverage:.1%})")

        # Calculate stats
        stats_df = self.calculate_stats(df)

        # Filter by total coverage
        valid_mask = stats_df["total_coverage_ratio"] >= min_coverage
        valid_targets = stats_df[valid_mask].index

        # Split hits
        passed, failed = self._split_by_targets(df, valid_targets)

        print(f"✅ Passed: {len(valid_targets)} targets ({len(passed)} hits)")
        print(f"❌ Failed: {len(stats_df) - len(valid_targets)} targets ({len(failed)} hits)")

        return FilterResult(
            passed=passed,
            failed=failed,
            stats=stats_df,
            filter_name="total_coverage"
        )

    def filter_by_hybrid_coverage(
        self,
        df: pd.DataFrame,
        min_total_cov: float = DEFAULT_MIN_TOTAL_COV,
        min_unmasked_cov: float = DEFAULT_MIN_UNMASKED_COV
    ) -> FilterResult:
        """
        Hybrid filter: targets pass if BOTH total AND unmasked coverage meet thresholds.

        This is a strict filter requiring both conditions:
        - Total coverage (including masked) >= min_total_cov
        - AND unmasked coverage >= min_unmasked_cov

        Args:
            df: DataFrame with alignment hits
            min_total_cov: Minimum total coverage threshold (default=0.5)
            min_unmasked_cov: Minimum unmasked coverage threshold (default=0.25)

        Returns:
            FilterResult object with passed/failed hits and statistics
        """
        if df is None or df.empty:
            return FilterResult(
                passed=df or pd.DataFrame(),
                failed=pd.DataFrame(),
                stats=pd.DataFrame(),
                filter_name="hybrid_coverage"
            )

        print(
            f"🔍 Filtering by hybrid coverage: "
            f"total≥{min_total_cov:.1%} AND unmasked≥{min_unmasked_cov:.1%}"
        )

        # Calculate stats
        stats_df = self.calculate_stats(df)

        # Hybrid filtering: BOTH conditions required (AND logic)
        total_pass = stats_df["total_coverage_ratio"] >= min_total_cov
        unmasked_pass = stats_df["unmasked_coverage_ratio"] >= min_unmasked_cov
        valid_mask = total_pass & unmasked_pass
        valid_targets = stats_df[valid_mask].index

        # Add pass reason to stats for analysis
        stats_df = stats_df.copy()
        stats_df["total_cov_pass"] = total_pass
        stats_df["unmasked_cov_pass"] = unmasked_pass
        stats_df["passed"] = valid_mask

        # Split hits
        passed, failed = self._split_by_targets(df, valid_targets)

        print(f"✅ Passed: {len(valid_targets)} targets ({len(passed)} hits)")
        print(f"❌ Failed: {len(stats_df) - len(valid_targets)} targets ({len(failed)} hits)")

        # Show breakdown
        both_pass = (total_pass & unmasked_pass).sum()
        total_only = (total_pass & ~unmasked_pass).sum()
        unmasked_only = (~total_pass & unmasked_pass).sum()
        neither = (~total_pass & ~unmasked_pass).sum()

        print(f"   📊 Breakdown: both={both_pass}, total_only={total_only}, "
              f"unmasked_only={unmasked_only}, neither={neither}")

        return FilterResult(
            passed=passed,
            failed=failed,
            stats=stats_df,
            filter_name="hybrid_coverage"
        )

    def __repr__(self) -> str:
        return (
            f"MaskingFilter(n_targets={self.n_targets_with_masks}, "
            f"overlap_threshold={self._overlap_threshold})"
        )
