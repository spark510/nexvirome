"""
Alignment result parser for MMseqs2 and BLAST outputs.
"""
from __future__ import annotations
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np

from ..core import AlignmentFormat, FilterCriteria, AlignmentError, log_info, log_verbose
from .header import ReadHeaderNormalizer


class AlignmentParser:
    """
    Parser for sequence alignment results.

    Supports:
    - MMseqs2 output format
    - BLAST tabular format (outfmt 6)
    - Automatic format detection
    - Vectorized filtering
    - Memory-efficient processing
    """

    # Column schemas for different formats
    # New format (no mismatch/gapopen, has taxname/taxlineage)
    MMSEQS_COLUMNS = [
        'query', 'target', 'fident', 'alnlen',
        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits',
        'qlen', 'tlen', 'taxid', 'taxname', 'taxlineage'
    ]

    MMSEQS_DTYPES = {
        'query': 'string',
        'target': 'string',
        'fident': 'float32',
        'alnlen': 'int32',
        'qstart': 'int32',
        'qend': 'int32',
        'tstart': 'int32',
        'tend': 'int32',
        'evalue': 'float64',
        'bits': 'float32',
        'qlen': 'int32',
        'tlen': 'int32',
        'taxid': 'int32',
        'taxname': 'string',
        'taxlineage': 'string'
    }

    # Legacy format (has mismatch/gapopen, no taxname/taxlineage)
    MMSEQS_LEGACY_COLUMNS = [
        'query', 'target', 'fident', 'alnlen',
        'mismatch', 'gapopen',
        'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits',
        'qlen', 'tlen', 'taxid'
    ]

    MMSEQS_LEGACY_DTYPES = {
        'query': 'string',
        'target': 'string',
        'fident': 'float32',
        'alnlen': 'int32',
        'mismatch': 'int32',
        'gapopen': 'int32',
        'qstart': 'int32',
        'qend': 'int32',
        'tstart': 'int32',
        'tend': 'int32',
        'evalue': 'float64',
        'bits': 'float32',
        'qlen': 'int32',
        'tlen': 'int32',
        'taxid': 'int32',
    }

    BLAST_COLUMNS = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore',
        'sgi', 'saccver', 'slen', 'qlen', 'staxids', 'sscinames',
        'stitle', 'qcovs', 'qcovhsp'
    ]

    def __init__(self, normalize_headers: bool = True):
        """
        Initialize parser.

        Args:
            normalize_headers: Normalize read headers for paired-end matching
        """
        self.normalize_headers = normalize_headers

    def parse(
        self,
        file_path: str | Path,
        format_type: Optional[AlignmentFormat] = None,
    ) -> pd.DataFrame:
        """
        Parse alignment file.

        Args:
            file_path: Path to alignment file
            format_type: Format type (auto-detected if None)

        Returns:
            DataFrame with alignment hits

        Raises:
            AlignmentError: If parsing fails
        """
        file_path = Path(file_path)

        if not file_path.exists():
            raise AlignmentError(f"File not found: {file_path}")

        # Parquet cache: same path with .parquet suffix
        parquet_cache = file_path.with_suffix('.parquet')
        if (parquet_cache.exists() and
                parquet_cache.stat().st_mtime >= file_path.stat().st_mtime):
            log_info(f"⚡ Loading cached parquet: {parquet_cache.name}")
            df = pd.read_parquet(parquet_cache)
            log_info(f"✅ Loaded {len(df):,} hits from cache")
            return df

        # Auto-detect format
        if format_type is None:
            format_type = self._detect_format(file_path)

        log_verbose(f"📄 Parsing {file_path.name} as {format_type.value} format...")

        try:
            if format_type == AlignmentFormat.MMSEQS:
                df = self._parse_mmseqs(file_path)
            elif format_type in (AlignmentFormat.BLAST, AlignmentFormat.BLAST6):
                df = self._parse_blast(file_path)
            else:
                raise AlignmentError(f"Unsupported format: {format_type}")

            # Calculate query coverage if not present
            if 'qcov' not in df.columns:
                df['qcov'] = (df['alnlen'] / df['qlen']).astype('float32')

            # Normalize headers
            if self.normalize_headers and 'query' in df.columns:
                df['query'] = ReadHeaderNormalizer.normalize_series(df['query'])

            # Save parquet cache
            try:
                df.to_parquet(parquet_cache, compression='zstd', index=False)
                log_verbose(f"💾 Cached to {parquet_cache.name}")
            except Exception:
                pass  # Cache write failure is non-fatal

            log_info(f"✅ Loaded {len(df):,} hits from {file_path.name}")
            return df

        except Exception as e:
            raise AlignmentError(f"Failed to parse {file_path}: {e}") from e

    def filter(
        self,
        df: pd.DataFrame,
        criteria: Optional[FilterCriteria] = None,
        **kwargs
    ) -> pd.DataFrame:
        """
        Filter alignment hits.

        Args:
            df: DataFrame with hits
            criteria: FilterCriteria object
            **kwargs: Individual filter parameters (override criteria)

        Returns:
            Filtered DataFrame
        """
        if df is None or df.empty:
            return pd.DataFrame()

        # Use provided criteria or create from kwargs
        if criteria is None:
            criteria = FilterCriteria(
                min_identity=kwargs.get('min_identity', 0.8),
                min_alignment_length=kwargs.get('min_alignment_length', 30),
                max_evalue=kwargs.get('max_evalue', 1e-3),
                min_query_coverage=kwargs.get('min_query_coverage', 0.5),
            )

        original_count = len(df)

        # Vectorized filtering
        mask = (
            (df['fident'] >= criteria.min_identity) &
            (df['alnlen'] >= criteria.min_alignment_length) &
            (df['evalue'] <= criteria.max_evalue) &
            (df['qcov'] >= criteria.min_query_coverage)
        )

        filtered = df[mask].copy()

        retention = len(filtered) / original_count * 100 if original_count > 0 else 0
        log_verbose(
            f"🔍 Filtered: {len(filtered):,}/{original_count:,} hits "
            f"({retention:.1f}% retention)"
        )

        return filtered

    # ========== Private Methods ==========

    def _detect_format(self, file_path: Path) -> AlignmentFormat:
        """Auto-detect file format from extension."""
        suffix = file_path.suffix.lower()

        if 'blast' in suffix or 'blast6' in file_path.name.lower():
            return AlignmentFormat.BLAST
        else:
            return AlignmentFormat.MMSEQS  # Default

    def _parse_mmseqs(self, file_path: Path) -> pd.DataFrame:
        """Parse MMseqs2 format. Auto-detects header and column layout."""
        with open(file_path, 'r') as f:
            lines = f.readlines()

        # Empty / header-only result (a sample with 0 viral hits): return an empty
        # frame with the standard columns instead of crashing on lines[1]. The
        # downstream classifier handles an empty hit set and writes empty outputs,
        # so a 0-read sample no longer breaks the run.
        first_line = lines[0].strip() if lines else ""
        has_header = first_line.startswith('query\t')
        data_lines = lines[1:] if has_header else lines
        if not any(ln.strip() for ln in data_lines):
            log_info(f"  ⚠️  {file_path.name}: no alignment hits (0 viral reads) — empty result")
            return pd.DataFrame(columns=self.MMSEQS_COLUMNS)

        data_line = next(ln.strip() for ln in data_lines if ln.strip())
        n_cols = len(data_line.split('\t'))

        # Detect format: legacy (mismatch/gapopen, 15 cols with int-like col5)
        # vs new (taxname/taxlineage, 15 cols with evalue-like col9)
        is_legacy = False
        if n_cols == 15 and not has_header:
            parts = data_line.split('\t')
            # In legacy format, col10 (index 10) is evalue (contains 'E' or very small float)
            # In new format, col8 (index 8) is evalue
            # Quick check: col4 (index 4) is mismatch (small int) in legacy, qstart (int) in new
            # Most reliable: col10 (index 10) — in legacy it's evalue (has 'E'), in new it's qlen (pure int)
            try:
                col10 = parts[10]
                if 'E' in col10.upper() or 'e' in col10 or '.' in col10:
                    is_legacy = True
            except (IndexError, ValueError):
                pass

        if is_legacy:
            columns = self.MMSEQS_LEGACY_COLUMNS
            dtypes = self.MMSEQS_LEGACY_DTYPES
            log_info(f"  Detected legacy MMseqs2 format (with mismatch/gapopen)")
        else:
            columns = self.MMSEQS_COLUMNS
            dtypes = self.MMSEQS_DTYPES

        df = pd.read_csv(
            file_path,
            sep='\t',
            header=0 if has_header else None,
            names=columns,
            skiprows=1 if has_header else 0,
            dtype=dtypes,
            engine='c',
            low_memory=False,
        )

        # If legacy format, drop mismatch/gapopen columns (not needed downstream)
        if is_legacy:
            df = df.drop(columns=['mismatch', 'gapopen'], errors='ignore')

        # Normalize identity to a 0-1 fraction. MMseqs2 'pident' is reported on a
        # 0-100 percent scale (e.g. 96.9), but the quality filter compares against
        # min_identity as a fraction (0.85). Without this, fident values of 68-100
        # always pass a 0.85 cut and the identity filter is silently a no-op.
        if 'fident' in df.columns and len(df) and df['fident'].max() > 1.0:
            df['fident'] = (df['fident'] / 100.0).astype('float32')

        return df

    def _parse_blast(self, file_path: Path) -> pd.DataFrame:
        """Parse BLAST format and normalize to MMseqs2 schema."""
        df = pd.read_csv(
            file_path,
            sep='\t',
            header=None,
            comment='#',
            names=self.BLAST_COLUMNS,
        )

        # Rename columns to match MMseqs2 schema
        df = df.rename(columns={
            'qseqid': 'query',
            'sseqid': 'target',
            'pident': 'fident',
            'length': 'alnlen',
            'sstart': 'tstart',
            'send': 'tend',
            'bitscore': 'bits',
            'slen': 'tlen',
            'staxids': 'taxid',
        })

        # Convert percent identity to fraction
        df['fident'] = (df['fident'] / 100.0).astype('float32')

        return df


class BatchAlignmentParser:
    """
    Parse multiple alignment files efficiently.

    Useful for paired-end data or batch processing.
    """

    def __init__(self, normalize_headers: bool = True):
        """
        Initialize batch parser.

        Args:
            normalize_headers: Normalize read headers
        """
        self.parser = AlignmentParser(normalize_headers=normalize_headers)

    def parse_paired(
        self,
        forward_file: str | Path,
        reverse_file: Optional[str | Path] = None,
        filter_criteria: Optional[FilterCriteria] = None,
    ) -> tuple[pd.DataFrame, Optional[pd.DataFrame]]:
        """
        Parse paired-end alignment files.

        Args:
            forward_file: Forward read alignments
            reverse_file: Reverse read alignments (optional)
            filter_criteria: Filtering criteria

        Returns:
            Tuple of (forward_df, reverse_df)
        """
        # Parse forward
        forward_df = self.parser.parse(forward_file)

        if filter_criteria:
            forward_df = self.parser.filter(forward_df, filter_criteria)

        # Parse reverse if provided
        reverse_df = None
        if reverse_file:
            reverse_df = self.parser.parse(reverse_file)
            if filter_criteria:
                reverse_df = self.parser.filter(reverse_df, filter_criteria)

        return forward_df, reverse_df
