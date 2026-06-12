"""
Masking filter for viral genome alignments.

This package provides OOP-based filtering of alignment hits based on
masked (repetitive/low-complexity) regions.

Quick Start:
    >>> from masking_filter import MaskingFilter
    >>>
    >>> # Create filter from BED file
    >>> filter = MaskingFilter.from_bed_file("masked_regions.bed")
    >>>
    >>> # Filter by unmasked coverage
    >>> result = filter.filter_by_unmasked_coverage(hits_df, min_coverage=0.25)
    >>> print(result.summary())

Main Classes:
    MaskingFilter: Main API for filtering
    FilterResult: Result object with passed/failed DataFrames
    CoverageStats: Coverage statistics for a single target
    MaskedRegion: Masked region data for a target

Utility Classes:
    CoverageCalculator: Low-level coverage calculations
    MaskLoader: Loading masked regions from files
"""

# Public API
from .filter import MaskingFilter, DEFAULT_MIN_UNMASKED_COV, DEFAULT_MIN_TOTAL_COV
from .models import FilterResult, CoverageStats, CoverageMetrics
from .mask_loader import MaskedRegion, MaskLoader
from .coverage_calculator import CoverageCalculator

__version__ = "2.0.0"

__all__ = [
    # Main API
    "MaskingFilter",

    # Data models
    "FilterResult",
    "CoverageStats",
    "CoverageMetrics",
    "MaskedRegion",

    # Utilities
    "CoverageCalculator",
    "MaskLoader",

    # Constants
    "DEFAULT_MIN_UNMASKED_COV",
    "DEFAULT_MIN_TOTAL_COV",
]
