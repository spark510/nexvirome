"""
Alignment parsing and filtering.

This module provides tools for:
- Parsing alignment files (MMseqs2, BLAST, etc.)
- Normalizing read headers
- Filtering based on quality criteria
- Masking filter for low-complexity regions
"""

from .parser import AlignmentParser, BatchAlignmentParser
from .header import ReadHeaderNormalizer

# Import all filtering tools from filters package
from .filters import (
    MaskingFilter,
    FilterResult,
    CoverageStats,
    CoverageMetrics,
    MaskedRegion,
    MaskLoader,
    CoverageCalculator,
    DEFAULT_MIN_UNMASKED_COV,
    DEFAULT_MIN_TOTAL_COV,
)

__all__ = [
    # Parsing
    "AlignmentParser",
    "BatchAlignmentParser",

    # Header normalization
    "ReadHeaderNormalizer",

    # Filtering - Main API
    "MaskingFilter",
    "FilterResult",
    "CoverageStats",
    "CoverageMetrics",
    "MaskedRegion",

    # Filtering - Utilities
    "MaskLoader",
    "CoverageCalculator",

    # Constants
    "DEFAULT_MIN_UNMASKED_COV",
    "DEFAULT_MIN_TOTAL_COV",
]
