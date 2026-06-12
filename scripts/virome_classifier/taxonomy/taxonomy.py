"""
Low-level taxonomy implementation with fast LCA computation.

This module provides the core Taxonomy class that loads taxonomic data from SQLite
and builds efficient data structures for:
- O(1) Lowest Common Ancestor (LCA) queries using Euler Tour + RMQ
- Fast lineage traversal
- Rank filtering and normalization
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple, Sequence
from collections import deque
from functools import lru_cache
import sqlite3

from ..core import log_verbose
from .ranks import MAJOR_RANKS, STANDARD_RANKS, normalize_rank


# -----------------------------
# Optional object graph (debug/navigation; NOT used on LCA hot path)
# -----------------------------
@dataclass
class TaxonNode:
    taxid: int
    parent_taxid: int
    rank: str
    name: str
    parent: Optional["TaxonNode"] = None

# -------------------- Types & constants --------------------

Row = Tuple[int, int, str, str]  # (taxid, parent_taxid, rank, name)


# -------------------- Main class --------------------

class Taxonomy:
    """
    Loads an NCBI or ICTV taxonomy subtree from SQLite (default: viruses root taxid=10239),
    builds a fast array-based LCA engine using Euler Tour + RMQ (Sparse Table),
    and optionally exposes a TaxonNode graph for navigation/debugging.
    """

    def __init__(self, sqlite_db_path: str, *,
                table: str = "ncbi_taxonomy",
                root_taxid: int = 10239,
                build_nodes: bool = False,
                pragmas: bool = True,
                build_major_taxon: bool = False
                ):

        self.db_path = sqlite_db_path  # Store for reproducibility
        self.nodes: List[TaxonNode] = []
        self.taxid_to_node: Dict[int, TaxonNode] = {}
        self._build_major_taxon = build_major_taxon
        self._major_default = build_major_taxon
        self.root_taxid = root_taxid

        rows = self._load_rows(sqlite_db_path, table, root_taxid, pragmas)
        # Normalize ranks to lowercase and map to standard names
        rows = [(t, p, normalize_rank(r or ""), n)
                for (t, p, r, n) in rows]

        if build_nodes:
            self._build_nodes(rows)
        self._build_lca_engine(rows)
        
        log_verbose(f"✅ Loaded taxonomy with {len(rows):,} nodes from '{table}' (root taxid={root_taxid})")


    # ---------- SQLite loading ----------
    def _load_rows(
        self,
        sqlite_db_path: str,
        table: str,
        root_taxid: int,
        pragmas: bool,
    ) -> List[Row]:
        sql = f"""
        WITH RECURSIVE subtree(taxid, parent_taxid, rank, scientific_name) AS (
          SELECT taxid, parent_taxid, rank, scientific_name
          FROM {table}
          WHERE taxid = {root_taxid}
          UNION ALL
          SELECT t.taxid, t.parent_taxid, t.rank, t.scientific_name
          FROM {table} t
          JOIN subtree v ON t.parent_taxid = v.taxid
        )
        SELECT taxid,
               parent_taxid,
               COALESCE(rank, 'no rank')           AS rank,
               COALESCE(scientific_name, '')       AS name
        FROM subtree;
        """
        con = sqlite3.connect(sqlite_db_path)
        try:
            if pragmas:
                con.execute("PRAGMA journal_mode=WAL;")
                con.execute("PRAGMA synchronous=NORMAL;")
                con.execute("PRAGMA temp_store=MEMORY;")
                con.execute(f"CREATE INDEX IF NOT EXISTS idx_{table}_taxid ON {table}(taxid);")
                con.execute(f"CREATE INDEX IF NOT EXISTS idx_{table}_parent_taxid ON {table}(parent_taxid);")
            rows = con.execute(sql).fetchall()
        finally:
            con.close()
        return rows

    # ---------- Optional: object graph (navigation/debug only) ----------
    def _build_nodes(self, rows: List[Row]) -> None:
        taxid_to_node: Dict[int, TaxonNode] = {}
        for t, p, r, n in rows:
            taxid_to_node[t] = TaxonNode(taxid=t, parent_taxid=p, rank=r, name=n)
        for node in taxid_to_node.values():
            node.parent = taxid_to_node.get(node.parent_taxid)
        self.taxid_to_node = taxid_to_node
        self.nodes = list(taxid_to_node.values())

    # ---------- Fast array-based LCA engine (Euler + RMQ) ----------
    def _build_lca_engine(self, rows: List[Row]) -> None:
        taxids = [t for t, _, _, _ in rows]
        idx_of_taxid = {t: i for i, t in enumerate(taxids)}
        N = len(rows)

        children: List[List[int]] = [[] for _ in range(N)]
        parent_idx: List[int] = [-1] * N
        root: Optional[int] = None
        for i, (t, p, _, _) in enumerate(rows):
            if p == t or p == 0 or p not in idx_of_taxid:
                root = i
                parent_idx[i] = i
            else:
                pi = idx_of_taxid[p]
                parent_idx[i] = pi
                children[pi].append(i)

        # Store ranks and names for nearest major ancestor computation
        ranks = [r for (_, _, r, _) in rows]
        names = [n for (_, _, _, n) in rows]

        # Store ranks and names for public access
        self._ranks = ranks
        self._names = names

        nearest_major = [0] * N
        if root is not None:
            q = deque([root])
            # Initialize root
            nearest_major[root] = root if ranks[root] in MAJOR_RANKS else root
            while q:
                v = q.popleft()
                base = v if ranks[v] in MAJOR_RANKS else nearest_major[parent_idx[v]]
                for w in children[v]:
                    nearest_major[w] = w if ranks[w] in MAJOR_RANKS else base
                    q.append(w)

        # Build major taxon tree if requested (filter to major ranks only)
        if self._build_major_taxon:
            # Extract major rank nodes only
            major_idxs = [i for i, rk in enumerate(ranks) if rk in MAJOR_RANKS]
            # Reconstruct rows with (taxid, parent_major_taxid, rank, name)
            rows_major: List[Row] = []
            for i in major_idxs:
                t = taxids[i]
                p_i = parent_idx[i]
                p_major_i = nearest_major[p_i]
                pt = t if p_major_i == i else taxids[p_major_i]
                rows_major.append((t, pt, ranks[i], names[i]))

            # Rebuild tree with major ranks only
            rows = rows_major
            taxids = [t for (t, _, _, _) in rows]
            idx_of_taxid = {t: i for i, t in enumerate(taxids)}
            N = len(rows)
            children = [[] for _ in range(N)]
            parent_idx = [-1] * N
            root = None
            for i, (t, p, _, _) in enumerate(rows):
                if p == t:
                    root = i
                    parent_idx[i] = i
                else:
                    pi = idx_of_taxid[p]
                    parent_idx[i] = pi
                    children[pi].append(i)
            ranks = [r for (_, _, r, _) in rows]
            names = [n for (_, _, _, n) in rows]

            # Update stored ranks and names after major taxon filtering
            self._ranks = ranks
            self._names = names


        # Euler tour (iterative), record first occurrence and depths
        euler: List[int] = []
        depth_e: List[int] = []
        first = [-1] * N
        depth = [0] * N

        if root is not None:
            # DFS to build Euler tour
            def dfs(v: int, d: int):
                if first[v] == -1:
                    first[v] = len(euler)
                euler.append(v)
                depth_e.append(d)

                for child in children[v]:
                    depth[child] = d + 1
                    dfs(child, d + 1)
                    # Revisit parent when returning from child
                    euler.append(v)
                    depth_e.append(d)

            depth[root] = 0
            dfs(root, 0)
                       

        # Sparse table over positions in euler
        M = len(depth_e)
        LOG = (M.bit_length())
        st: List[List[int]] = [[0] * M for _ in range(LOG)]
        for i in range(M):
            st[0][i] = i
        k = 1
        while (1 << k) <= M:
            half = 1 << (k - 1)
            prev, cur = st[k - 1], st[k]
            for i in range(M - (1 << k) + 1):
                a = prev[i]; b = prev[i + half]
                cur[i] = a if depth_e[a] <= depth_e[b] else b
            k += 1

        # Store computed fields
        self.idx_of_taxid = idx_of_taxid
        self.taxids = taxids
        self.first = first
        self.euler = euler
        self._depth_euler = depth_e
        self._st = st
        self._parent_idx = parent_idx

        # LCA caching: {frozenset(taxids): lca_result}
        self._lca_cache: Dict[frozenset, int] = {}
        self._lca_cache_hits = 0
        self._lca_cache_misses = 0

    # ---------- Public accessor methods ----------
    def name(self, taxid: int) -> Optional[str]:
        """Get the name of a taxid"""
        if taxid not in self.idx_of_taxid:
            return None
        idx = self.idx_of_taxid[taxid]
        return self._names[idx]
    
    def rank(self, taxid: int) -> Optional[str]:
        """Get the rank of a taxid"""
        if taxid not in self.idx_of_taxid:
            return None
        idx = self.idx_of_taxid[taxid]
        return self._ranks[idx]
    
    def parent(self, taxid: int) -> Optional[int]:
        """Get the parent taxid of a taxid"""
        if taxid not in self.idx_of_taxid:
            return None
        idx = self.idx_of_taxid[taxid]
        parent_idx = self._parent_idx[idx]
        if parent_idx == idx:  # root node
            return taxid
        return self.taxids[parent_idx]
    
    def lineage(self, taxid: int) -> List[int]:
        """Get the full lineage from root to taxid"""
        if taxid not in self.idx_of_taxid:
            return []
        
        lineage = []
        current_idx = self.idx_of_taxid[taxid]
        
        while True:
            lineage.append(self.taxids[current_idx])
            parent_idx = self._parent_idx[current_idx]
            if parent_idx == current_idx:  # reached root
                break
            current_idx = parent_idx
        
        return lineage[::-1]  # reverse to get root-to-taxid order
    
    def exists(self, taxid: int) -> bool:
        """Check if a taxid exists in the taxonomy"""
        return taxid in self.idx_of_taxid

    # ---------- Internal RMQ & 2-node LCA ----------
    def _rmq_pos(self, l: int, r: int) -> int:
        if l > r:
            l, r = r, l
        k = (r - l + 1).bit_length() - 1
        a = self._st[k][l]
        b = self._st[k][r - (1 << k) + 1]
        return a if self._depth_euler[a] <= self._depth_euler[b] else b

    def _lca_idx(self, u: int, v: int) -> int:
        pos = self._rmq_pos(self.first[u], self.first[v])
        return self.euler[pos]

    # ---------- Public LCA APIs ----------
    def lca_k_taxids(self, tids: Sequence[int]) -> int:
        """
        Return LCA taxid of a list of taxids (fold O(k)).

        Uses caching to avoid redundant calculations for common taxid combinations.
        """
        valid_tids = [t for t in tids if t > 0 and t in self.idx_of_taxid]

        if len(valid_tids) == 0:
            return 0
        if len(valid_tids) == 1:
            return valid_tids[0]

        # Create cache key from sorted taxids (order doesn't matter for LCA)
        cache_key = frozenset(valid_tids)

        # Check cache first
        if cache_key in self._lca_cache:
            self._lca_cache_hits += 1
            return self._lca_cache[cache_key]

        # Cache miss - compute LCA
        self._lca_cache_misses += 1

        it = iter(valid_tids)
        a = self.idx_of_taxid[next(it)]
        for t in it:
            a = self._lca_idx(a, self.idx_of_taxid[t])

        result = self.taxids[a]

        # Store in cache
        self._lca_cache[cache_key] = result

        return result

    def get_full_lineage(
        self,
        taxid: int,
        *,
        major: Optional[bool] = None,
        with_meta: bool = False,
    ):
        """
        Get full lineage from root to taxid.

        Args:
            taxid: Taxon ID
            major: Filter to major ranks only (None uses default from initialization)
            with_meta: Return list of dicts with taxid, rank, name

        Returns:
            List of taxids or list of dicts if with_meta=True
        """
        i = self.idx_of_taxid.get(taxid)
        if i is None:
            return []

        parent_idx = self._parent_idx
        taxids = self.taxids
        ranks = self._ranks
        names = self._names
        root_idx = getattr(self, "_root_idx", None)

        # Collect path from taxid to root
        path_idx: List[int] = []
        seen = 0
        N = len(parent_idx)
        while True:
            path_idx.append(i)
            if root_idx is not None and i == root_idx:
                break
            p = parent_idx[i]
            if p == i:  # self-parent = root
                break
            i = p
            seen += 1
            if seen > N:  # Safety check for cycles
                break

        path_idx.reverse()

        # Apply major rank filter if requested
        use_major = self._major_default if major is None else major
        if use_major:
            path_idx = [j for j in path_idx if self._is_major_rank_idx(j)]

        # Return results
        if with_meta:
            return [{"taxid": taxids[j], "rank": ranks[j], "name": names[j]} for j in path_idx]
        else:
            return [taxids[j] for j in path_idx]

    def _is_major_rank_idx(self, idx: int) -> bool:
        """Check if index corresponds to a major rank."""
        return self._ranks[idx] in MAJOR_RANKS

    def is_major_rank_taxid(self, taxid: int) -> bool:
        """Check if given taxid has a major rank"""
        idx = self.idx_of_taxid.get(taxid)
        if idx is None:
            return False
        return self._is_major_rank_idx(idx)

    def get_lca_cache_stats(self) -> Dict[str, int]:
        """
        Get LCA cache statistics.

        Returns:
            Dictionary with cache hits, misses, and hit rate
        """
        total = self._lca_cache_hits + self._lca_cache_misses
        hit_rate = (self._lca_cache_hits / total * 100) if total > 0 else 0.0

        return {
            "cache_size": len(self._lca_cache),
            "hits": self._lca_cache_hits,
            "misses": self._lca_cache_misses,
            "total_queries": total,
            "hit_rate_percent": hit_rate,
        }

    def clear_lca_cache(self) -> None:
        """Clear the LCA cache and reset statistics."""
        self._lca_cache.clear()
        self._lca_cache_hits = 0
        self._lca_cache_misses = 0
