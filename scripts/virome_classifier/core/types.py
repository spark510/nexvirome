"""Common type definitions for virome_classifier."""

from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass
from enum import Enum
import numpy as np


# Type aliases
TaxID = int
QueryName = str
TargetName = str


class AlignmentFormat(Enum):
    """Supported alignment output formats."""
    MMSEQS = "mmseqs"
    BLAST = "blast"
    BLAST6 = "blast6"


class PairedMode(Enum):
    """Paired-end resolution modes."""
    CONSERVATIVE = "conservative"  # Only pairs that agree
    COMPREHENSIVE = "comprehensive"  # Union of all hits


class WeightMode(Enum):
    """Weight calculation modes for LCA voting."""
    EVALUE = "evalue"
    BITSCORE = "bitscore"
    UNIFORM = "uniform"


@dataclass(frozen=True)
class AlignmentHit:
    """
    Single alignment hit (immutable).

    Represents one alignment between a query and target sequence.
    """
    query: str
    target: str
    identity: float  # 0.0 - 1.0
    alignment_length: int
    evalue: float
    bitscore: float
    query_length: int
    target_length: int
    taxid: int

    @property
    def query_coverage(self) -> float:
        """Query coverage ratio."""
        return self.alignment_length / self.query_length if self.query_length > 0 else 0.0

    @property
    def target_coverage(self) -> float:
        """Target coverage ratio."""
        return self.alignment_length / self.target_length if self.target_length > 0 else 0.0


@dataclass(frozen=True)
class FilterCriteria:
    """
    Criteria for filtering alignment hits (immutable).
    """
    min_identity: float = 0.8  # 80%
    min_alignment_length: int = 30
    max_evalue: float = 1e-3
    min_query_coverage: float = 0.5  # 50%

    def passes(self, hit: AlignmentHit) -> bool:
        """Check if hit passes all criteria."""
        return (
            hit.identity >= self.min_identity and
            hit.alignment_length >= self.min_alignment_length and
            hit.evalue <= self.max_evalue and
            hit.query_coverage >= self.min_query_coverage
        )


@dataclass(frozen=True)
class ClassificationResult:
    """
    Result of LCA classification for a single query (immutable).
    """
    query: str
    taxid: TaxID
    read_count: int  # 1 for single-end, 2 for paired-end
    confidence: float = 1.0
    method: str = "lca"

    @property
    def is_classified(self) -> bool:
        """Check if query was successfully classified."""
        return self.taxid > 0


@dataclass
class TaxonInfo:
    """
    Taxonomic information for a taxon (mutable for caching).
    """
    taxid: TaxID
    name: str
    rank: str
    parent_taxid: Optional[TaxID] = None
    lineage: Optional[List[TaxID]] = None

    def __hash__(self):
        """Make hashable for use in sets/dicts."""
        return hash(self.taxid)

    def __eq__(self, other):
        """Equality based on taxid."""
        if not isinstance(other, TaxonInfo):
            return False
        return self.taxid == other.taxid
