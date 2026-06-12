"""
Data models for masking filter.

This module contains value objects and data classes for coverage statistics
and filter results.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import pandas as pd


@dataclass(frozen=True)
class CoverageMetrics:
    """
    Immutable metrics for a coverage category (total/masked/unmasked).

    Attributes:
        hit_count: Number of alignment hits
        breadth_bp: Breadth of coverage in base pairs (positions covered)
        coverage_ratio: Fraction of target covered (merged intervals)
        avg_depth: Average sequencing depth across covered positions
    """
    hit_count: int
    breadth_bp: int
    coverage_ratio: float
    avg_depth: float

    @property
    def breadth_ratio(self) -> float:
        """Alias for coverage_ratio (for clarity in some contexts)."""
        return self.coverage_ratio

    @classmethod
    def empty(cls, reason: str = "no_hits") -> 'CoverageMetrics':
        """Create empty metrics (no coverage)."""
        return cls(
            hit_count=0,
            breadth_bp=0,
            coverage_ratio=0.0,
            avg_depth=0.0
        )

    def to_dict(self, prefix: str = "") -> dict:
        """
        Convert to dictionary for DataFrame creation.

        Args:
            prefix: Prefix for keys (e.g., "total_", "masked_")

        Returns:
            Dictionary with prefixed keys
        """
        return {
            f"{prefix}hit_count": self.hit_count,
            f"{prefix}breadth_bp": self.breadth_bp,
            f"{prefix}coverage_ratio": self.coverage_ratio,
            f"{prefix}breadth_ratio": self.coverage_ratio,  # Same as coverage_ratio
            f"{prefix}avg_depth": self.avg_depth,
        }


@dataclass(frozen=True)
class CoverageStats:
    """
    Comprehensive coverage statistics for a single target sequence.

    This is an immutable value object containing all coverage metrics
    for total, masked, and unmasked regions.

    Attributes:
        target_name: Target sequence accession (e.g., "NC_001806.2")
        target_length: Total length of target sequence in bp
        taxid: NCBI taxonomy ID (optional)
        has_mask: Whether this target has masked region data
        total: Metrics for all hits (masked + unmasked)
        masked: Metrics for hits in masked regions only
        unmasked: Metrics for hits NOT in masked regions
    """
    target_name: str
    target_length: int
    taxid: Optional[int]
    has_mask: bool
    total: CoverageMetrics
    masked: CoverageMetrics
    unmasked: CoverageMetrics

    def to_series(self) -> pd.Series:
        """
        Convert to pandas Series for backwards compatibility.

        Returns:
            Series with all metrics as a flat structure
        """
        data = {
            "target_length": self.target_length,
            "has_mask": self.has_mask,
            "taxid": self.taxid,
        }

        # Add metrics with prefixes
        data.update(self.total.to_dict(prefix="total_"))
        data.update(self.masked.to_dict(prefix="masked_"))
        data.update(self.unmasked.to_dict(prefix="unmasked_"))

        # Legacy aliases for backwards compatibility
        data["total_coverage"] = self.total.coverage_ratio
        data["unmasked_coverage"] = self.unmasked.coverage_ratio
        data["hit_count"] = self.total.hit_count
        data["masked_hit_count"] = self.masked.hit_count

        return pd.Series(data, name=self.target_name)

    @classmethod
    def from_series(cls, series: pd.Series, target_name: str) -> 'CoverageStats':
        """
        Create from pandas Series (for backwards compatibility).

        Args:
            series: Series with coverage metrics
            target_name: Name of the target sequence

        Returns:
            CoverageStats object
        """
        return cls(
            target_name=target_name,
            target_length=int(series["target_length"]),
            taxid=series.get("taxid"),
            has_mask=bool(series["has_mask"]),
            total=CoverageMetrics(
                hit_count=int(series["total_hit_count"]),
                breadth_bp=int(series["total_breadth_bp"]),
                coverage_ratio=float(series["total_coverage_ratio"]),
                avg_depth=float(series["total_avg_depth"]),
            ),
            masked=CoverageMetrics(
                hit_count=int(series["masked_hit_count"]),
                breadth_bp=int(series["masked_breadth_bp"]),
                coverage_ratio=float(series["masked_coverage_ratio"]),
                avg_depth=float(series["masked_avg_depth"]),
            ),
            unmasked=CoverageMetrics(
                hit_count=int(series["unmasked_hit_count"]),
                breadth_bp=int(series["unmasked_breadth_bp"]),
                coverage_ratio=float(series["unmasked_coverage_ratio"]),
                avg_depth=float(series["unmasked_avg_depth"]),
            ),
        )


@dataclass(frozen=True)
class FilterResult:
    """
    Result of a filtering operation.

    Attributes:
        passed: DataFrame of hits that passed the filter
        failed: DataFrame of hits that failed the filter
        stats: DataFrame with per-target statistics
        filter_name: Name of the filter applied
    """
    passed: pd.DataFrame
    failed: pd.DataFrame
    stats: pd.DataFrame
    filter_name: str

    @property
    def n_passed_targets(self) -> int:
        """Number of unique targets that passed."""
        return self.passed["target"].nunique() if not self.passed.empty else 0

    @property
    def n_failed_targets(self) -> int:
        """Number of unique targets that failed."""
        return self.failed["target"].nunique() if not self.failed.empty else 0

    @property
    def n_passed_hits(self) -> int:
        """Number of hits that passed."""
        return len(self.passed)

    @property
    def n_failed_hits(self) -> int:
        """Number of hits that failed."""
        return len(self.failed)

    def summary(self) -> str:
        """
        Generate a human-readable summary of the filter results.

        Returns:
            Multi-line summary string
        """
        lines = [
            f"{'='*60}",
            f"Filter: {self.filter_name}",
            f"{'='*60}",
            f"✅ Passed: {self.n_passed_targets} targets ({self.n_passed_hits} hits)",
            f"❌ Failed: {self.n_failed_targets} targets ({self.n_failed_hits} hits)",
        ]

        if not self.stats.empty:
            # Add statistics if available
            if "total_coverage_ratio" in self.stats.columns:
                avg_cov = self.stats["total_coverage_ratio"].mean()
                lines.append(f"📊 Average coverage: {avg_cov:.1%}")

        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"FilterResult(filter='{self.filter_name}', "
            f"passed={self.n_passed_targets} targets, "
            f"failed={self.n_failed_targets} targets)"
        )
