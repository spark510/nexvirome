"""
Legacy function wrappers for backwards compatibility.

This module provides wrapper functions that match the old API,
allowing existing code to work without changes.

DEPRECATED: These functions are provided for backwards compatibility only.
New code should use the MaskingFilter class directly.
"""
from __future__ import annotations

import warnings
from typing import Dict, Tuple
import pandas as pd
import numpy as np

from .filter import MaskingFilter
from .mask_loader import MaskLoader


def load_masked_bed(bed_file: str) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
    """
    Load masked regions from BED file (legacy format).

    DEPRECATED: Use MaskingFilter.from_bed_file() instead.

    Args:
        bed_file: Path to BED file

    Returns:
        Dictionary mapping accession to (starts_array, ends_array)
    """
    warnings.warn(
        "load_masked_bed() is deprecated. Use MaskingFilter.from_bed_file() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    mask_dict = MaskLoader.from_bed_file(bed_file)
    return MaskLoader.to_legacy_format(mask_dict)


def apply_unmasked_cov_filter(
    df: pd.DataFrame,
    mask_dict: Dict[str, Tuple[np.ndarray, np.ndarray]],
    min_unmasked_cov: float = 0.25,
    coverage_col: str = "qcov"
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Filter targets based on unmasked region coverage (legacy wrapper).

    DEPRECATED: Use MaskingFilter.filter_by_unmasked_coverage() instead.

    Args:
        df: MMseqs hits DataFrame
        mask_dict: Dictionary of masked regions (legacy format)
        min_unmasked_cov: Minimum unmasked coverage ratio
        coverage_col: Column name for coverage (unused, for compatibility)

    Returns:
        Tuple of (passed_df, failed_df, target_stats)
    """
    warnings.warn(
        "apply_unmasked_cov_filter() is deprecated. "
        "Use MaskingFilter.filter_by_unmasked_coverage() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    # Convert legacy format to new format
    mask_dict_new = MaskLoader.from_legacy_format(mask_dict)

    # Use new API
    filter_obj = MaskingFilter(mask_dict_new)
    result = filter_obj.filter_by_unmasked_coverage(df, min_coverage=min_unmasked_cov)

    return result.passed, result.failed, result.stats


def apply_total_cov_filter(
    df: pd.DataFrame,
    mask_dict: Dict[str, Tuple[np.ndarray, np.ndarray]],
    min_total_cov: float = 0.5,
    coverage_col: str = "qcov"
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Filter targets based on total coverage (legacy wrapper).

    DEPRECATED: Use MaskingFilter.filter_by_total_coverage() instead.

    Args:
        df: MMseqs hits DataFrame
        mask_dict: Dictionary of masked regions (legacy format)
        min_total_cov: Minimum total coverage ratio
        coverage_col: Column name for coverage (unused, for compatibility)

    Returns:
        Tuple of (passed_df, failed_df, target_stats)
    """
    warnings.warn(
        "apply_total_cov_filter() is deprecated. "
        "Use MaskingFilter.filter_by_total_coverage() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    # Convert legacy format to new format
    mask_dict_new = MaskLoader.from_legacy_format(mask_dict)

    # Use new API
    filter_obj = MaskingFilter(mask_dict_new)
    result = filter_obj.filter_by_total_coverage(df, min_coverage=min_total_cov)

    return result.passed, result.failed, result.stats


def apply_hybrid_cov_filter(
    df: pd.DataFrame,
    mask_dict: Dict[str, Tuple[np.ndarray, np.ndarray]],
    min_total_cov: float = 0.5,
    min_unmasked_cov: float = 0.2,
    coverage_col: str = "qcov"
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Hybrid filter: targets pass if BOTH conditions met (legacy wrapper).

    DEPRECATED: Use MaskingFilter.filter_by_hybrid_coverage() instead.

    Args:
        df: MMseqs hits DataFrame
        mask_dict: Dictionary of masked regions (legacy format)
        min_total_cov: Minimum total coverage threshold
        min_unmasked_cov: Minimum unmasked coverage threshold
        coverage_col: Column name for coverage (unused, for compatibility)

    Returns:
        Tuple of (passed_df, failed_df, target_stats)
    """
    warnings.warn(
        "apply_hybrid_cov_filter() is deprecated. "
        "Use MaskingFilter.filter_by_hybrid_coverage() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    # Convert legacy format to new format
    mask_dict_new = MaskLoader.from_legacy_format(mask_dict)

    # Use new API
    filter_obj = MaskingFilter(mask_dict_new)
    result = filter_obj.filter_by_hybrid_coverage(
        df,
        min_total_cov=min_total_cov,
        min_unmasked_cov=min_unmasked_cov
    )

    return result.passed, result.failed, result.stats


def calculate_coverage_stats(
    group: pd.DataFrame,
    mask_dict: Dict[str, Tuple[np.ndarray, np.ndarray]],
    overlap_threshold: float = 0.5
) -> pd.Series:
    """
    Calculate comprehensive coverage statistics (legacy wrapper).

    DEPRECATED: Use CoverageCalculator.calculate_stats() instead.

    Args:
        group: DataFrame of hits for a single target
        mask_dict: Dictionary of masked regions (legacy format)
        overlap_threshold: Minimum overlap ratio to consider hit as masked

    Returns:
        Series with comprehensive coverage statistics
    """
    warnings.warn(
        "calculate_coverage_stats() is deprecated. "
        "Use CoverageCalculator.calculate_stats() instead.",
        DeprecationWarning,
        stacklevel=2
    )

    from .coverage_calculator import CoverageCalculator

    # Convert legacy format
    mask_dict_new = MaskLoader.from_legacy_format(mask_dict)

    # Use new API
    calculator = CoverageCalculator(overlap_threshold=overlap_threshold)
    stats = calculator.calculate_stats(group, mask_dict_new)

    return stats.to_series()


# Alias for backwards compatibility
apply_masking_cov_filter = apply_unmasked_cov_filter
