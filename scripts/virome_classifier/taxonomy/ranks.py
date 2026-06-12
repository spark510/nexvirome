"""
Taxonomic rank definitions and utilities.
"""

from typing import Set, Dict


# Standard taxonomic ranks
STANDARD_RANKS: Set[str] = {
    "acellular root",
    "realm",
    "kingdom",
    "phylum",
    "subphylum",
    "class",
    "order",
    "suborder",
    "family",
    "subfamily",
    "genus",
    "subgenus",
    "species",
}

# Major ranks (commonly reported)
MAJOR_RANKS: Set[str] = {
    "acellular root",
    "realm",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
}

# Rank levels (for hierarchy)
RANK_LEVELS: Dict[str, int] = {
    "acellular root": 0,
    "realm": 1,
    "kingdom": 2,
    "phylum": 3,
    "subphylum": 4,
    "class": 5,
    "subclass": 6,
    "order": 7,
    "suborder": 8,
    "family": 9,
    "subfamily": 10,
    "genus": 11,
    "subgenus": 12,
    "species": 13,
    "subspecies": 14,
    "strain": 15,
    "no rank": 99,
}

# Canonical rank codes (for Kraken format)
RANK_CODES: Dict[str, str] = {
    "superkingdom": "D",
    "realm": "D",  # Treat realm as domain for viruses
    "kingdom": "K",
    "phylum": "P",
    "class": "C",
    "order": "O",
    "family": "F",
    "genus": "G",
    "species": "S",
    "subspecies": "SS",
    "strain": "S1",
    "no rank": "-",
    "acellular root": "-",
}

# Normalize rank names (NCBI virus variations)
RANK_NORMALIZATION: Dict[str, str] = {
    "clade": "no rank",
    "genotype": "no rank",
    "isolate": "no rank",
    "serogroup": "no rank",
    "serotype": "no rank",
    "unranked": "no rank",
}


def normalize_rank(rank: str) -> str:
    """
    Normalize rank name to standard form.

    Args:
        rank: Raw rank name

    Returns:
        Normalized rank name
    """
    if not rank:
        return "no rank"

    rank_lower = rank.lower().strip()
    return RANK_NORMALIZATION.get(rank_lower, rank_lower)


def is_major_rank(rank: str) -> bool:
    """
    Check if rank is a major rank.

    Args:
        rank: Rank name

    Returns:
        True if major rank
    """
    normalized = normalize_rank(rank)
    return normalized in MAJOR_RANKS


def get_rank_level(rank: str) -> int:
    """
    Get numeric level for rank (lower = higher in hierarchy).

    Args:
        rank: Rank name

    Returns:
        Numeric level (0 = root, 99 = unranked)
    """
    normalized = normalize_rank(rank)
    return RANK_LEVELS.get(normalized, 99)


def get_rank_code(rank: str) -> str:
    """
    Get Kraken-style rank code.

    Args:
        rank: Rank name

    Returns:
        Single-letter rank code (e.g., 'G' for genus)
    """
    normalized = normalize_rank(rank)
    return RANK_CODES.get(normalized, "-")


def compare_ranks(rank1: str, rank2: str) -> int:
    """
    Compare two ranks.

    Args:
        rank1: First rank
        rank2: Second rank

    Returns:
        -1 if rank1 < rank2 (rank1 higher in hierarchy)
        0 if equal
        1 if rank1 > rank2 (rank1 lower in hierarchy)
    """
    level1 = get_rank_level(rank1)
    level2 = get_rank_level(rank2)

    if level1 < level2:
        return -1
    elif level1 > level2:
        return 1
    else:
        return 0
