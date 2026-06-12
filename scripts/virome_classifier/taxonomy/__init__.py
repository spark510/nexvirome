"""Taxonomy database and rank utilities."""

from .taxonomy import Taxonomy
from .taxonomy_db import TaxonomyDB, SegmentInfo
from .ranks import (
    MAJOR_RANKS,
    STANDARD_RANKS,
    normalize_rank,
    is_major_rank,
    get_rank_level,
    get_rank_code,
    compare_ranks,
)

__all__ = [
    # Core classes
    "Taxonomy",
    "TaxonomyDB",
    "SegmentInfo",

    # Rank utilities
    "MAJOR_RANKS",
    "STANDARD_RANKS",
    "normalize_rank",
    "is_major_rank",
    "get_rank_level",
    "get_rank_code",
    "compare_ranks",
]
