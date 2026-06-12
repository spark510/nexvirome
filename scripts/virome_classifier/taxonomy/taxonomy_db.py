"""
OOP wrapper for taxonomy database.

This wraps the low-level Taxonomy class from the original codebase
with a clean, type-safe interface.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple
from pathlib import Path
from dataclasses import dataclass
import sqlite3

from ..core import TaxID, TaxonInfo, TaxonomyError, log_info, log_verbose
from .ranks import normalize_rank, is_major_rank
from .taxonomy import Taxonomy


@dataclass
class SegmentInfo:
    """Information about a viral reference sequence (segment)."""
    segment_id: str           # segment_name from DB (e.g., "RNA 1", "L", "S")
    segment_length: int       # sequence length
    species_taxid: int        # NCBI taxid
    virus_id: str             # Groups segments of the same virus (taxid as string)
    is_segmented: bool = False  # True if virus has multiple segments

class TaxonomyDB:   
    """
    High-level taxonomy database interface.

    Wraps the low-level Taxonomy class with a clean OOP API.

    Example:
        >>> tax = TaxonomyDB.from_sqlite("taxonomy.db")
        >>> name = tax.get_name(10239)  # Viruses
        >>> lineage = tax.get_lineage(10359)  # Human betaherpesvirus 5
    """

    def __init__(self, taxonomy: Taxonomy, db_path: Optional[Path] = None):
        """
        Initialize from existing Taxonomy instance.

        Args:
            taxonomy: Original Taxonomy instance
            db_path: Path to SQLite database (for additional queries)
        """
        self._tax = taxonomy
        self._db_path = db_path
        self._cache: dict = {}  # Simple cache for frequently accessed data

    @classmethod
    def from_sqlite(
        cls,
        db_path: str | Path,
        root_taxid: int = 10239,  # Viruses
        table: str = "ncbi_taxonomy",
        build_major_taxon: bool = False,
    ) -> "TaxonomyDB":
        """
        Load taxonomy from SQLite database.

        Args:
            db_path: Path to SQLite database
            root_taxid: Root taxid to load (default: 10239 = Viruses)
            table: Table name in database
            build_major_taxon: Build major taxon tree only

        Returns:
            TaxonomyDB instance

        Raises:
            TaxonomyError: If loading fails
        """
        db_path = Path(db_path)
        if not db_path.exists():
            raise TaxonomyError(f"Taxonomy database not found: {db_path}")

        # In-process cache: loading the tree (~282k taxa) takes ~3.8s; reuse it
        # when the same DB is loaded again in one process (in-process batch over
        # many samples, ML eval loops, etc.). Does NOT help subprocess-per-sample
        # runs — those each get a fresh interpreter.
        cache = getattr(cls, "_load_cache", None)
        if cache is None:
            cache = cls._load_cache = {}
        cache_key = (str(db_path.resolve()), int(root_taxid), table, bool(build_major_taxon))
        if cache_key in cache:
            return cache[cache_key]

        try:
            log_info(f"📚 Loading taxonomy from {db_path} (root={root_taxid})...")

            tax = Taxonomy(
                str(db_path),
                table=table,
                root_taxid=root_taxid,
                build_nodes=False,  # Don't build slow object graph
                build_major_taxon=build_major_taxon,
            )

            log_info(f"✅ Loaded {len(tax.taxids):,} taxa")
            instance = cls(tax, db_path=db_path)
            cache[cache_key] = instance
            return instance

        except Exception as e:
            raise TaxonomyError(f"Failed to load taxonomy: {e}") from e

    # ========== Basic Accessors ==========

    def exists(self, taxid: TaxID) -> bool:
        """Check if taxid exists in database."""
        return self._tax.exists(taxid)

    def get_name(self, taxid: TaxID) -> Optional[str]:
        """Get scientific name for taxid."""
        if not self.exists(taxid):
            return None
        return self._tax.name(taxid)

    def get_rank(self, taxid: TaxID) -> Optional[str]:
        """Get rank for taxid."""
        if not self.exists(taxid):
            return None
        rank = self._tax.rank(taxid)
        return normalize_rank(rank) if rank else "no rank"

    def get_parent(self, taxid: TaxID) -> Optional[TaxID]:
        """Get parent taxid."""
        if not self.exists(taxid):
            return None
        return self._tax.parent(taxid)

    def get_lineage(self, taxid: TaxID, major_only: bool = False) -> List[TaxID]:
        """
        Get full lineage from root to taxid.

        Args:
            taxid: Taxon ID
            major_only: Return only major ranks

        Returns:
            List of taxids from root to taxid
        """
        if not self.exists(taxid):
            return []

        lineage = self._tax.lineage(taxid)

        if major_only:
            # Filter to major ranks only
            lineage = [
                tid for tid in lineage
                if is_major_rank(self.get_rank(tid) or "no rank")
            ]

        return lineage

    def get_taxid_at_rank(self, taxid: int, target_rank: str) -> int:
        """
        Return the taxid at a chosen higher rank (genus, family, order, etc.)

        Walks up the lineage starting from the input taxid until it finds
        an ancestor with rank == target_rank.

        If not found, returns the original taxid.

        Args:
            taxid: Starting taxonomic ID (strain, species, genus, etc.)
            target_rank: Desired rank, e.g. 'genus', 'family', 'order'

        Returns:
            The taxid of the ancestor at the given rank, or the original taxid.

        Example:
            >>> tax.get_taxid_at_rank(10886, "genus")
            10880   # Mammalian orthoreovirus (genus)
        """

        # Invalid taxid
        if not self.exists(taxid):
            return taxid

        # Normalize rank string
        target_rank = target_rank.lower()

        # Cache key
        cache_key = f"{target_rank}_{taxid}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        current = taxid
        max_iterations = 20  # safety

        for _ in range(max_iterations):
            rank = self.get_rank(current)

            # Found matching rank
            if rank == target_rank:
                self._cache[cache_key] = current
                return current

            # Move to parent
            parent = self.get_parent(current)
            if parent is None or parent == current:
                break

            current = parent

        # No match, return original taxid
        self._cache[cache_key] = taxid
        return taxid


    def get_info(self, taxid: TaxID) -> Optional[TaxonInfo]:
        """
        Get complete taxonomic information.

        Args:
            taxid: Taxon ID

        Returns:
            TaxonInfo object or None if not found
        """
        if not self.exists(taxid):
            return None

        # Check cache
        cache_key = f"info_{taxid}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        # Build TaxonInfo
        info = TaxonInfo(
            taxid=taxid,
            name=self.get_name(taxid) or str(taxid),
            rank=self.get_rank(taxid) or "no rank",
            parent_taxid=self.get_parent(taxid),
            lineage=self.get_lineage(taxid),
        )

        # Cache it
        self._cache[cache_key] = info
        return info

    # ========== LCA Operations ==========

    def compute_lca(self, taxids: List[TaxID]) -> Optional[TaxID]:
        """
        Compute Lowest Common Ancestor (LCA) of multiple taxids.

        Args:
            taxids: List of taxon IDs

        Returns:
            LCA taxid, or None if no valid LCA exists
        """
        if not taxids:
            return None

        # Filter out invalid taxids
        valid_taxids = [tid for tid in taxids if tid > 0 and self.exists(tid)]

        if not valid_taxids:
            return None

        if len(valid_taxids) == 1:
            return valid_taxids[0]

        # Use fast LCA from original implementation
        try:
            lca_taxid = self._tax.lca_k_taxids(valid_taxids)
            return lca_taxid if lca_taxid > 0 else None
        except Exception as e:
            log_verbose(f"Warning: LCA computation failed: {e}")
            return None

    def compute_lca_with_rank(
        self,
        taxids: List[TaxID],
        min_rank: Optional[str] = None,
        major_only: bool = True,
    ) -> Optional[TaxID]:
        """
        Compute LCA and optionally snap to a minimum rank.

        Args:
            taxids: List of taxon IDs
            min_rank: Minimum rank to snap to (e.g., "genus")
            major_only: Snap to major ranks only

        Returns:
            LCA taxid at appropriate rank
        """
        lca = self.compute_lca(taxids)
        if lca is None:
            return None

        # If major_only, find nearest major rank ancestor
        if major_only:
            current = lca
            while current is not None:
                rank = self.get_rank(current)
                if rank and is_major_rank(rank):
                    return current
                current = self.get_parent(current)
            return lca  # Fallback to original LCA

        return lca

    # ========== Utility Methods ==========

    def is_ancestor(self, ancestor: TaxID, descendant: TaxID) -> bool:
        """
        Check if ancestor is an ancestor of descendant.

        Args:
            ancestor: Potential ancestor taxid
            descendant: Descendant taxid

        Returns:
            True if ancestor is in lineage of descendant
        """
        if not self.exists(ancestor) or not self.exists(descendant):
            return False

        lineage = self.get_lineage(descendant)
        return ancestor in lineage

    def get_all_taxa(self) -> Set[TaxID]:
        """Get all taxids in database."""
        return set(self._tax.taxids)

    def __len__(self) -> int:
        """Number of taxa in database."""
        return len(self._tax.taxids)

    def __repr__(self) -> str:
        return f"TaxonomyDB(taxa={len(self):,}, root={self._tax.root_taxid})"

    # ========== LCA Cache Statistics ==========

    def get_lca_cache_stats(self) -> dict:
        """
        Get LCA cache statistics.

        Returns:
            Dictionary with cache hits, misses, and hit rate
        """
        return self._tax.get_lca_cache_stats()

    def clear_lca_cache(self) -> None:
        """Clear the LCA cache and reset statistics."""
        self._tax.clear_lca_cache()

    # ========== Segment Info ==========

    def get_segment_info(
        self,
        table: str = "refseq_sequences",
    ) -> Dict[str, SegmentInfo]:
        """
        Load segment information from refseq_sequences table.

        This method queries the database for reference sequence metadata
        including segment names and lengths, which is useful for coverage-based
        classification of segmented viruses.

        Args:
            table: Table name containing reference sequences (default: "refseq_sequences")

        Returns:
            Dict mapping accession -> SegmentInfo

        Raises:
            TaxonomyError: If database path not available or query fails

        Example:
            >>> tax = TaxonomyDB.from_sqlite("tax_seq.db")
            >>> segment_info = tax.get_segment_info()
            >>> print(f"Loaded {len(segment_info)} references")
            >>> # Use with CoverageBasedClassifier
            >>> classifier = CoverageBasedClassifier(tax, segment_info=segment_info)
        """
        if self._db_path is None:
            raise TaxonomyError(
                "Database path not available. "
                "Use TaxonomyDB.from_sqlite() to load with segment info support."
            )

        # Check cache
        cache_key = f"segment_info_{table}"
        if cache_key in self._cache:
            return self._cache[cache_key]

        try:
            log_info(f"📋 Loading segment info from {table}...")

            conn = sqlite3.connect(str(self._db_path))
            cursor = conn.cursor()

            # Query refseq_sequences table. The old `is_cds=0` filter (excluding
            # title-flagged single-CDS records) was removed: is_cds was a title-keyword
            # guess that mis-classified many complete genomes, so all references are
            # now included. (n_cds replaces is_cds but is informational, not a filter.)
            query = f"""
                SELECT accession, taxid, length, segment_name, is_segmented
                FROM {table}
            """
            cursor.execute(query)
            rows = cursor.fetchall()
            conn.close()

            # Build segment_info dict
            segment_info: Dict[str, SegmentInfo] = {}

            for accession, taxid, length, segment_name, is_segmented in rows:
                # Normalize to species-level taxid
                raw_taxid = int(taxid) if taxid else 0
                species_taxid = self.get_taxid_at_rank(raw_taxid, "species") if raw_taxid > 0 else 0

                # Use species taxid as virus_id to group segments of the same virus
                virus_id = str(species_taxid)

                # segment_name can be None for non-segmented viruses
                segment_id = segment_name if segment_name else accession

                segment_info[accession] = SegmentInfo(
                    segment_id=segment_id,
                    segment_length=int(length) if length else 0,
                    species_taxid=species_taxid,
                    virus_id=virus_id,
                    is_segmented=bool(is_segmented),
                )

            # Cache the result
            self._cache[cache_key] = segment_info

            # Count segmented viruses
            segmented_count = sum(1 for s in segment_info.values() if s.is_segmented)
            log_info(f"✅ Loaded {len(segment_info):,} references "
                    f"({segmented_count:,} segmented)")

            return segment_info

        except sqlite3.Error as e:
            raise TaxonomyError(f"Failed to load segment info: {e}") from e

    def search_by_name(self, query: str) -> List[Tuple[int, str]]:
        """
        Search for taxa by name (case-insensitive substring match).

        Args:
            query: Name substring to search for

        Returns:
            List of (taxid, name) tuples matching the query
        """
        query = query.lower()
        results = []

        # Use internal arrays for speed if available
        if hasattr(self._tax, 'taxids') and hasattr(self._tax, '_names'):
            for taxid, name in zip(self._tax.taxids, self._tax._names):
                if query in name.lower():
                    results.append((taxid, name))
        else:
            for taxid in self.get_all_taxa():
                name = self.get_name(taxid)
                if name and query in name.lower():
                    results.append((taxid, name))

        return sorted(results, key=lambda x: x[1])
