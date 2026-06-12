"""
Rank code utilities for Kraken report generation.
"""
from __future__ import annotations


def get_rank_code(rank: str) -> str:
    """
    Convert taxonomic rank to Kraken rank code.

    Kraken uses single-letter codes for standard ranks.
    Note: 'realm' is mapped to 'D' (Domain) for virus taxonomy.

    Args:
        rank: Taxonomic rank name (e.g., 'genus', 'species')

    Returns:
        Single-letter rank code or '-' for non-standard ranks

    Examples:
        >>> get_rank_code('genus')
        'G'
        >>> get_rank_code('species')
        'S'
        >>> get_rank_code('no rank')
        '-'
    """
    r = (rank or "no rank").lower()
    rank_codes = {
        "superkingdom": "D",
        "realm": "D",  # Virus taxonomy: treat realm as domain
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
    return rank_codes.get(r, "-")
