"""
Masked region loading and management.

This module handles loading masked regions from BED files and provides
efficient overlap checking.
"""
from __future__ import annotations

import io
from typing import Dict, Tuple
import numpy as np
import pandas as pd
from pathlib import Path

from . import mask_crypt


class MaskedRegion:
    """
    Represents masked (repetitive/low-complexity) regions for a single target sequence.

    Masked regions are stored as parallel numpy arrays for efficient vectorized
    overlap calculations.

    Attributes:
        target_name: Target sequence accession
        starts: Array of region start positions (0-based, inclusive)
        ends: Array of region end positions (0-based, exclusive / half-open)
    """

    def __init__(self, target_name: str, starts: np.ndarray, ends: np.ndarray):
        """
        Initialize masked regions for a target.

        Args:
            target_name: Target sequence accession
            starts: Array of start positions
            ends: Array of end positions

        Raises:
            ValueError: If starts and ends have different lengths
        """
        if len(starts) != len(ends):
            raise ValueError(
                f"Starts and ends must have same length: "
                f"got {len(starts)} vs {len(ends)}"
            )

        self.target_name = target_name
        self.starts = np.asarray(starts, dtype=np.int32)
        self.ends = np.asarray(ends, dtype=np.int32)

    @property
    def n_regions(self) -> int:
        """Number of masked regions."""
        return len(self.starts)

    def calculate_overlap_ratio(
        self,
        hit_start: int,
        hit_end: int
    ) -> float:
        """
        Calculate maximum overlap ratio between a hit and any masked region.

        Args:
            hit_start: Hit start position (0-based, inclusive)
            hit_end: Hit end position (0-based, exclusive)

        Returns:
            Maximum fraction of hit that overlaps with any masked region (0.0-1.0)
        """
        if self.n_regions == 0:
            return 0.0

        hit_length = hit_end - hit_start  # 0-based half-open

        # Calculate overlap with each region (vectorized, half-open intervals)
        overlap_starts = np.maximum(hit_start, self.starts)
        overlap_ends = np.minimum(hit_end, self.ends)
        overlap_lengths = np.maximum(0, overlap_ends - overlap_starts)

        # Return max overlap ratio
        max_overlap = overlap_lengths.max()
        return float(max_overlap / hit_length) if hit_length > 0 else 0.0

    def overlaps_with(
        self,
        hit_start: int,
        hit_end: int,
        threshold: float = 0.5
    ) -> bool:
        """
        Check if a hit overlaps with any masked region above threshold.

        Args:
            hit_start: Hit start position
            hit_end: Hit end position
            threshold: Minimum overlap ratio to consider as "masked" (default=0.5)

        Returns:
            True if hit overlaps >= threshold with any masked region
        """
        overlap_ratio = self.calculate_overlap_ratio(hit_start, hit_end)
        return overlap_ratio >= threshold

    def __repr__(self) -> str:
        return f"MaskedRegion(target='{self.target_name}', n_regions={self.n_regions})"


class MaskLoader:
    """
    Factory class for loading masked regions from BED files.

    This class provides static methods to load masked region data from
    various file formats.
    """

    @staticmethod
    def from_bed_file(bed_file: str) -> Dict[str, MaskedRegion]:
        """
        Load masked regions from BED file.

        BED format (tab-separated):
        Column 0: accession (target name)
        Column 1: start position (0-based)
        Column 2: end position (0-based, exclusive)

        Note: BED uses 0-based half-open intervals [start, end).
        We keep coordinates as-is (0-based half-open) to match
        MMseqs2 output and Python array slicing conventions.

        The file may be either a plaintext BED or an encrypted mask container
        (see mask_crypt). Encrypted files are auto-detected by their magic header
        and decrypted transparently in memory, so the curated mask is always
        applied with the exact same call; the rest of the pipeline is unaware of
        the distinction. Development can use plaintext, distribution the .enc.

        Args:
            bed_file: Path to BED file (plaintext or encrypted)

        Returns:
            Dictionary mapping target name to MaskedRegion object

        Raises:
            FileNotFoundError: If BED file doesn't exist
            ValueError: If BED file has invalid format (or wrong decryption key)
        """
        bed_path = Path(bed_file)
        if not bed_path.exists():
            raise FileNotFoundError(f"BED file not found: {bed_file}")

        # Encrypted mask -> decrypt to text in memory; plaintext -> read as-is.
        # pd.read_csv accepts a StringIO just like a path, so the parse body below
        # is identical for both branches.
        if mask_crypt.is_encrypted(bed_path):
            bed_source: object = io.StringIO(mask_crypt.decrypt_to_bed_text(bed_path))
        else:
            bed_source = bed_file

        try:
            # Read BED (only first 3 columns needed)
            bed_df = pd.read_csv(
                bed_source,
                sep="\t",
                header=None,
                names=["acc", "start", "end", "depth", "desc", "other_acc", "other_desc"],
                usecols=[0, 1, 2],  # acc, start, end
                dtype={"acc": "string", "start": "int32", "end": "int32"}
            )

            if bed_df.empty:
                return {}

            # Group by accession and create MaskedRegion objects
            mask_dict = {}
            for acc, group in bed_df.groupby("acc"):
                # Sort by start position for consistency
                intervals = group[["start", "end"]].sort_values("start").to_numpy()

                # Keep BED 0-based half-open coordinates as-is
                # This matches MMseqs2 tstart/tend (0-based) and Python slicing
                starts = intervals[:, 0]
                ends = intervals[:, 1]

                mask_dict[acc] = MaskedRegion(
                    target_name=acc,
                    starts=starts,
                    ends=ends
                )

            return mask_dict

        except Exception as e:
            raise ValueError(f"Error parsing BED file: {e}") from e

    @staticmethod
    def from_dataframe(df: pd.DataFrame) -> Dict[str, MaskedRegion]:
        """
        Create masked regions from a DataFrame.

        Args:
            df: DataFrame with columns ['target', 'start', 'end']
               (0-based, half-open coordinates)

        Returns:
            Dictionary mapping target name to MaskedRegion object
        """
        required_cols = {"target", "start", "end"}
        if not required_cols.issubset(df.columns):
            raise ValueError(
                f"DataFrame must have columns: {required_cols}, "
                f"got: {set(df.columns)}"
            )

        mask_dict = {}
        for target, group in df.groupby("target"):
            intervals = group[["start", "end"]].sort_values("start").to_numpy()
            mask_dict[target] = MaskedRegion(
                target_name=target,
                starts=intervals[:, 0],
                ends=intervals[:, 1]
            )

        return mask_dict

    @staticmethod
    def to_legacy_format(
        mask_dict: Dict[str, MaskedRegion]
    ) -> Dict[str, Tuple[np.ndarray, np.ndarray]]:
        """
        Convert to legacy tuple format for backwards compatibility.

        Args:
            mask_dict: Dictionary of MaskedRegion objects

        Returns:
            Dictionary mapping target to (starts_array, ends_array)
        """
        return {
            target: (region.starts, region.ends)
            for target, region in mask_dict.items()
        }

    @staticmethod
    def from_legacy_format(
        legacy_dict: Dict[str, Tuple[np.ndarray, np.ndarray]]
    ) -> Dict[str, MaskedRegion]:
        """
        Create from legacy tuple format.

        Args:
            legacy_dict: Dictionary mapping target to (starts, ends) tuples

        Returns:
            Dictionary of MaskedRegion objects
        """
        return {
            target: MaskedRegion(target, starts, ends)
            for target, (starts, ends) in legacy_dict.items()
        }
