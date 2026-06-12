"""
Virome Classifier - Advanced viral metagenomics classification pipeline.

A modern, OOP-based refactoring of the virome classification system with:
- Clean separation of concerns (core, taxonomy, alignment, classification, reporting)
- Type-safe interfaces with comprehensive type hints
- Immutable data classes for reliability
- Comprehensive documentation and examples
- Full test coverage

Quick Start:
    >>> from virome_classifier import (
    ...     TaxonomyDB, AlignmentParser, MaskingFilter,
    ...     FilterCriteria, set_verbose
    ... )
    >>>
    >>> # Enable verbose logging
    >>> set_verbose(True)
    >>>
    >>> # Load taxonomy database
    >>> tax = TaxonomyDB.from_sqlite("taxonomy.db")
    >>>
    >>> # Parse and filter alignments
    >>> parser = AlignmentParser()
    >>> df = parser.parse("alignments.tsv")
    >>> criteria = FilterCriteria(min_identity=0.8, min_alignment_length=30)
    >>> filtered = parser.filter(df, criteria)
    >>>
    >>> # Apply masking filter
    >>> masking = MaskingFilter.from_bed_file("masked.bed")
    >>> result = masking.filter_by_unmasked_coverage(filtered, min_coverage=0.25)
    >>> print(result.summary())
"""

from .__version__ import __version__, __author__, __description__

# ========== Core Utilities ==========
from .core import (
    # Logging
    ViromeLogger,
    get_logger,
    set_verbose,
    log_info,
    log_verbose,

    # Exceptions
    ViromeClassifierError,
    TaxonomyError,
    AlignmentError,
    ClassificationError,

    # Types
    FilterCriteria,
)

# ========== Taxonomy ==========
from .taxonomy import TaxonomyDB, SegmentInfo

# ========== Alignment ==========
from .alignment import (
    # Parsers
    AlignmentParser,
    BatchAlignmentParser,

    # Filters
    MaskingFilter,
    FilterResult,
    CoverageStats,
    CoverageMetrics,
)

# ========== Reporting ==========
from .reporting import (
    write_all_outputs,
    write_kraken_output,
    generate_kraken_report,
    write_abundance_table,
)

# # ========== Classification ==========
# from .coverage_based_classifier2 import (
#     CoverageBasedClassifier2,
    
#     HitQualityFilter,
#     ReadCandidate,
#     SegmentCoverage,
#     SpeciesCoverage,
#     CoverageThresholds,
# )

from .coverage_based_classifier3 import (
    CoverageBasedClassifier3,
    HitQualityFilter,
    SegmentCoverage,
    StrainCoverage,
    TaxonCoverage,
    CoverageThresholds,
)

# from .coverage_based_classifier4 import (
#     CoverageBasedClassifier4,
#     HitQualityFilter,
#     SegmentCoverage,
#     TaxonCoverage,
#     CoverageThresholds,
# )

# ========== Public API ==========
__all__ = [
    # Version info
    "__version__",
    "__author__",
    "__description__",

    # ===== Core =====
    # Logging
    "ViromeLogger",
    "get_logger",
    "set_verbose",
    "log_info",
    "log_verbose",

    # Exceptions
    "ViromeClassifierError",
    "TaxonomyError",
    "AlignmentError",
    "ClassificationError",

    # Types
    "FilterCriteria",

    # ===== Main Classes =====
    # Taxonomy
    "TaxonomyDB",
    "SegmentInfo",

    # Alignment
    "AlignmentParser",
    "BatchAlignmentParser",

    # Filters
    "MaskingFilter",
    "FilterResult",
    "CoverageStats",
    "CoverageMetrics",

    # ===== Reporting =====
    "write_all_outputs",
    "write_kraken_output",
    "generate_kraken_report",
    "write_abundance_table",

    # ===== Classification =====
    # "CoverageBasedClassifier2",
    "CoverageBasedClassifier3",

    "HitQualityFilter",
    "ReadCandidate",
    "SegmentCoverage",
    # "SpeciesCoverage",
    "TaxonCoverage",
    "StrainCoverage",
    "CoverageThresholds",
]
