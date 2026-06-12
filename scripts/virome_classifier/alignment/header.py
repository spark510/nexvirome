"""
Read header normalization for paired-end matching.

Handles Illumina sequencing headers and other common formats.
"""

import re
from typing import Optional
import pandas as pd


class ReadHeaderNormalizer:
    """
    Normalizes read headers for paired-end matching.

    Removes paired-end indicators (/1, /2, _R1, _R2, etc.) to enable
    matching of forward and reverse reads.
    """

    # Pattern to match and remove paired-end indicators
    PAIRED_END_PATTERN = re.compile(
        r"^@|"  # Leading '@'
        r"\s+[12]:[NnYy]:\d+:.*$|"  # Casava 1.8+ format
        r"#0/[12]$|"  # Old Illumina format
        r"[/._][12]$"  # Common suffixes: /1, .1, _1
    )

    @classmethod
    def normalize(cls, header: str) -> str:
        """
        Normalize a single header string.

        Args:
            header: Raw header string

        Returns:
            Normalized header (paired-end indicators removed)

        Example:
            >>> ReadHeaderNormalizer.normalize("read1/1")
            'read1'
            >>> ReadHeaderNormalizer.normalize("@M00967:1:1101:1000:1000 1:N:0:1")
            'M00967:1:1101:1000:1000'
        """
        if not header:
            return header

        return re.sub(cls.PAIRED_END_PATTERN, '', header)

    @classmethod
    def normalize_series(cls, headers: pd.Series) -> pd.Series:
        """
        Normalize a pandas Series of headers (vectorized).

        Args:
            headers: Series of header strings

        Returns:
            Series with normalized headers
        """
        return headers.str.replace(cls.PAIRED_END_PATTERN, '', regex=True)

    @classmethod
    def extract_sample_name(cls, filename: str) -> str:
        """
        Extract sample name from filename.

        Removes common extensions and paired-end indicators.

        Args:
            filename: Input filename

        Returns:
            Sample name

        Example:
            >>> ReadHeaderNormalizer.extract_sample_name("sample1_R1.fastq.gz")
            'sample1'
        """
        import os

        basename = os.path.basename(filename) if filename else ""

        # Remove common extensions
        for ext in ['.mmseqs', '.tsv', '.txt', '.gz', '.fastq', '.fasta', '.blast6']:
            basename = basename.replace(ext, '')

        # Remove paired-end indicators
        for indicator in ['_R1', '_R2', '.R1', '.R2', '_1', '_2', '/1', '/2']:
            basename = basename.replace(indicator, '')

        return basename
