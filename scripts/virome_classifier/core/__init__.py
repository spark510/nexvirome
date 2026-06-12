"""Core utilities and base classes for virome_classifier."""

from .logger import (
    ViromeLogger,
    get_logger,
    set_verbose,
    log_info,
    log_verbose,
    log_debug,
    log_warning,
    log_error,
    log_success,
)

from .exceptions import (
    ViromeClassifierError,
    TaxonomyError,
    AlignmentError,
    ClassificationError,
    InvalidInputError,
    ConfigurationError,
)

from .types import (
    TaxID,
    QueryName,
    TargetName,
    AlignmentFormat,
    PairedMode,
    WeightMode,
    AlignmentHit,
    FilterCriteria,
    ClassificationResult,
    TaxonInfo,
)

__all__ = [
    # Logger
    "ViromeLogger",
    "get_logger",
    "set_verbose",
    "log_info",
    "log_verbose",
    "log_debug",
    "log_warning",
    "log_error",
    "log_success",

    # Exceptions
    "ViromeClassifierError",
    "TaxonomyError",
    "AlignmentError",
    "ClassificationError",
    "InvalidInputError",
    "ConfigurationError",

    # Types
    "TaxID",
    "QueryName",
    "TargetName",
    "AlignmentFormat",
    "PairedMode",
    "WeightMode",
    "AlignmentHit",
    "FilterCriteria",
    "ClassificationResult",
    "TaxonInfo",
]
