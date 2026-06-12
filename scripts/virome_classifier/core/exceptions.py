"""Custom exceptions for virome_classifier."""


class ViromeClassifierError(Exception):
    """Base exception for all virome_classifier errors."""
    pass


class TaxonomyError(ViromeClassifierError):
    """Raised when taxonomy operations fail."""
    pass


class AlignmentError(ViromeClassifierError):
    """Raised when alignment parsing/filtering fails."""
    pass


class ClassificationError(ViromeClassifierError):
    """Raised when classification fails."""
    pass


class InvalidInputError(ViromeClassifierError):
    """Raised when input data is invalid."""
    pass


class ConfigurationError(ViromeClassifierError):
    """Raised when configuration is invalid."""
    pass
