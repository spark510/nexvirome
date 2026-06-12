"""
Coverage-Based Virome Classifier v4
Depth Entropy + Masking-Aware Coverage Integrated Version
"""

from typing import Dict, Set, Tuple, Optional
import pandas as pd
import numpy as np

# Data models + pure helpers live in coverage_models.py (extracted for clarity).
# Imported into this namespace so existing references and external imports
# (e.g. `from ...coverage_based_classifier3 import SegmentCoverage`) keep working.
try:
    from virome_classifier.coverage_models import (
        calculate_breadth_fast, compute_depth_entropy,
        SegmentCoverage, StrainCoverage, TaxonCoverage,
        CoverageThresholds, HitQualityFilter, HAS_NUMBA,
    )
except ImportError:  # support running as a top-level module
    from coverage_models import (
        calculate_breadth_fast, compute_depth_entropy,
        SegmentCoverage, StrainCoverage, TaxonCoverage,
        CoverageThresholds, HitQualityFilter, HAS_NUMBA,
    )


# =============================================================================
# Classifier
# =============================================================================

class CoverageBasedClassifier3:

    def __init__(
        self,
        taxonomy_db,
        segment_info=None,
        hit_filter=None,
        thresholds=None,
        mask_bed_file=None,
        multi_mapping_mode="all",
        classification_rank="strain",
        use_depth_entropy=False,           # NEW
        simple_mode=False,                 # NEW: Simple unique/multi ratio mode
        blacklist_taxa=None,               # NEW: Blacklist taxa (names or IDs)
        breadth_denominator="aligned",     # NEW: "aligned" (hit segments only,
                                           # legacy) or "expected" (taxon whole-
                                           # genome length from DB; drops partial-
                                           # /short-reference FPs like Shamonda S)
        # Local depth parameters
        local_depth_alpha=0.5,             # Weight for alignment score
        local_depth_beta=0.5,              # Weight for depth score
        local_depth_temperature=1.0,       # Softmax temperature
        local_depth_iterations=2,          # Number of EM iterations
        verbose=False
    ):
        self.tax = taxonomy_db
        self.segment_info = segment_info or {}
        self.hit_filter = hit_filter or HitQualityFilter()
        self.thresholds = thresholds or CoverageThresholds()
        self.verbose = verbose
        self.use_depth_entropy = use_depth_entropy     # NEW
        self.simple_mode = simple_mode                 # NEW

        # Blacklist setup
        self.blacklist_names = set()
        self.blacklist_taxids = set()
        if blacklist_taxa:
            self._setup_blacklist(blacklist_taxa)

        valid_ranks = {"strain", "species", "genus", "family"}
        if classification_rank not in valid_ranks:
            raise ValueError(classification_rank)
        self.classification_rank = classification_rank

        valid_modes = {"all", "best_hit", "local_depth"}
        if multi_mapping_mode not in valid_modes:
            raise ValueError(multi_mapping_mode)
        self.multi_mapping_mode = multi_mapping_mode

        # Local depth parameters
        self.local_depth_alpha = local_depth_alpha
        self.local_depth_beta = local_depth_beta
        self.local_depth_temperature = local_depth_temperature
        self.local_depth_iterations = local_depth_iterations

        self.segment_lookup_df = None
        self._expected_segment_counts = None   # cached O(N) strain->#segments
        if self.segment_info:
            self._build_segment_lookup()

        # Breadth denominator policy. "expected" divides breadth by the taxon's
        # whole-genome length (refseq_sequences.expected_genome_length), so a read
        # set that lands only on the short segment of a multipartite virus does not
        # inflate breadth. Loaded once as strain_taxid -> expected_genome_length.
        valid_den = {"aligned", "expected"}
        if breadth_denominator not in valid_den:
            raise ValueError(breadth_denominator)
        self.breadth_denominator = breadth_denominator

        # Masking
        self.mask_filter = None
        self.masked_targets = set()
        if mask_bed_file:
            from .alignment import MaskingFilter
            self.mask_filter = MaskingFilter.from_bed_file(mask_bed_file)
            self.masked_targets = set(self.mask_filter._mask_dict.keys())

        # Expected-genome-length map (loaded last so self.tax/verbose are set).
        self._expected_genome_length = {}      # strain_taxid -> bp
        if breadth_denominator == "expected":
            self._expected_genome_length = self._load_expected_genome_length()

    def _load_expected_genome_length(self):
        """Map taxid -> expected_genome_length from the DB. The coverage loop keys
        strains by the (often species-level) taxid the reads roll up to — e.g.
        reads on Shamonda (NCBI 159150) are reported under its species
        Orthobunyavirus schmallenbergense (NCBI 3052437). The expected lengths are
        stored per refseq taxid (159150), so we ALSO register the value under the
        species ancestor (3052437), summing distinct segments across the species so
        a multipartite virus gets its whole-genome length at species level too.
        Empty dict (=> aligned fallback) if the column is absent."""
        out = {}
        conn = getattr(self.tax, "conn", None) or getattr(self.tax, "_conn", None)
        if conn is None:
            db_path = getattr(self.tax, "db_path", None) or getattr(self.tax, "_db_path", None)
            if db_path is not None:
                import sqlite3
                conn = sqlite3.connect(str(db_path))
        if conn is None:
            return out
        try:
            cols = {r[1] for r in conn.execute("PRAGMA table_info(refseq_sequences)")}
            if "expected_genome_length" not in cols:
                if self.verbose:
                    print("  [breadth_denominator=expected] column missing; "
                          "falling back to aligned per-taxon")
                return out
            # per refseq taxid (the value compute_expected_genome_length stored)
            for taxid, egl in conn.execute(
                "SELECT taxid, MAX(expected_genome_length) FROM refseq_sequences "
                "WHERE expected_genome_length IS NOT NULL GROUP BY taxid"):
                if taxid is not None and egl:
                    out[int(taxid)] = int(egl)
        except Exception as e:
            if self.verbose:
                print(f"  [breadth_denominator=expected] load failed: {e}; using aligned")
            return out
        # also register under each taxon's species ancestor, taking the LARGEST
        # whole-genome length among the refseq taxa under that species (so the
        # denominator is never smaller than any single isolate's full genome).
        species_egl = {}
        for taxid, egl in list(out.items()):
            sp = None
            try:
                sp = self.tax.get_taxid_at_rank(taxid, "species")
            except Exception:
                sp = None
            if sp:
                species_egl[int(sp)] = max(species_egl.get(int(sp), 0), int(egl))
        for sp, egl in species_egl.items():
            # don't clobber a directly-stored value that is already >=
            out[sp] = max(out.get(sp, 0), egl)
        return out

    def _setup_blacklist(self, blacklist_taxa):
        """
        Setup blacklist from taxon names or IDs.

        Args:
            blacklist_taxa: List of taxon names (str) or IDs (int)
        """
        if isinstance(blacklist_taxa, (list, tuple, set)):
            blacklist_list = blacklist_taxa
        else:
            blacklist_list = [blacklist_taxa]

        for item in blacklist_list:
            if isinstance(item, int):
                # Direct taxon ID
                self.blacklist_taxids.add(item)
            elif isinstance(item, str):
                # Taxon name - try to resolve to ID
                self.blacklist_names.add(item.lower())

                # Try to find taxid by name
                try:
                    taxid = self.tax.get_taxid_by_name(item)
                    if taxid:
                        self.blacklist_taxids.add(taxid)
                        if self.verbose:
                            print(f"  [BLACKLIST] Resolved '{item}' → TaxID {taxid}")
                except:
                    # If name resolution fails, keep in blacklist_names for later matching
                    if self.verbose:
                        print(f"  [BLACKLIST] Added name pattern: '{item}'")
                    pass

        if self.verbose and (self.blacklist_taxids or self.blacklist_names):
            print(f"Blacklist initialized: {len(self.blacklist_taxids)} taxids, {len(self.blacklist_names)} name patterns")

    def _is_blacklisted(self, taxid, taxon_name=""):
        """
        Check if a taxon is blacklisted.

        Args:
            taxid: Taxon ID
            taxon_name: Taxon name (optional)

        Returns:
            bool: True if blacklisted
        """
        # Check taxid
        if taxid in self.blacklist_taxids:
            return True

        # Check name patterns
        if taxon_name and self.blacklist_names:
            taxon_name_lower = taxon_name.lower()
            for pattern in self.blacklist_names:
                if pattern in taxon_name_lower:
                    return True

        return False

    # Build segment lookup
    def _build_segment_lookup(self):
        rows = []
        for acc, info in self.segment_info.items():
            strain_taxid = info.strain_taxid if hasattr(info, "strain_taxid") else info.species_taxid

            rows.append({
                'target': acc,
                'ref_strain_taxid': strain_taxid,
                'ref_species_taxid': self.tax.get_taxid_at_rank(strain_taxid, "species"),
                'ref_genus_taxid': self.tax.get_taxid_at_rank(strain_taxid, "genus"),
                'ref_family_taxid': self.tax.get_taxid_at_rank(strain_taxid, "family"),
                'ref_virus_id': str(strain_taxid),
                'ref_segment_id': info.segment_id,
                'ref_segment_length': info.segment_length
            })

        self.segment_lookup_df = pd.DataFrame(rows).set_index("target")

    # =========================================================================
    # Step1 (분리된 필터링 단계)
    # =========================================================================
    def step1a_quality_filter(self, hits_df: pd.DataFrame) -> pd.DataFrame:
        """
        Step 1a: Hit quality filtering only.

        Args:
            hits_df: Input hits DataFrame

        Returns:
            Filtered hits DataFrame
        """
        if self.verbose:
            print("Step 1a: Quality filtering...")

        filtered = self.hit_filter.filter_df(hits_df)

        if self.verbose:
            pct = len(filtered) / len(hits_df) * 100 if len(hits_df) > 0 else 0
            print(f"  {len(filtered):,} / {len(hits_df):,} hits passed ({pct:.1f}%)")

        return filtered

    def step1b_contamination_filter(self, hits_df: pd.DataFrame) -> pd.DataFrame:
        """
        Step 1b: Contamination sequence filtering only.

        Removes hits to known contamination sequences (NC_032111.1: BeAn 58058 virus).

        Args:
            hits_df: Input hits DataFrame

        Returns:
            Filtered hits DataFrame
        """
        if self.verbose:
            print("Step 1b: Contamination filtering...")

        # Known contamination sequences - BeAn 58058 virus
        CONTAMINATION_SEQUENCES = {
            'NC_032111.1',  # BeAn 58058 virus complete genome
            'NC_022518.1',
        }

        if 'target' in hits_df.columns:
            before_contam = len(hits_df)
            filtered = hits_df[~hits_df['target'].isin(CONTAMINATION_SEQUENCES)].copy()
            removed_contam = before_contam - len(filtered)

            if self.verbose and removed_contam > 0:
                print(f"  Removed {removed_contam:,} contamination sequence hits "
                      f"(NC_032111.1: BeAn 58058 virus)")
        else:
            filtered = hits_df.copy()
            if self.verbose:
                print("  [WARNING] 'target' column not found, skipping contamination filter")

        return filtered

    def step1c_masking_filter(self, hits_df: pd.DataFrame) -> pd.DataFrame:
        """
        Step 1c: Masking-based filtering only.

        Removes hits that heavily overlap with masked regions (if use_depth_entropy=False).

        Args:
            hits_df: Input hits DataFrame

        Returns:
            Filtered hits DataFrame
        """
        if self.verbose:
            print("Step 1c: Masking filtering...")

        filtered = hits_df.copy()

        # 3) MaskingFilter 기반 hit 제거는
        #    - MaskingFilter가 존재하고
        #    - use_depth_entropy == False 일 때만 적용
        if self.mask_filter:
            before = len(filtered)

            try:
                result = self.mask_filter.filter_by_unmasked_coverage(
                    filtered,
                    min_coverage=self.thresholds.min_unmasked_coverage
                )

                if not self.use_depth_entropy:
                    # Depth entropy를 사용하지 않을 때만 hit-level 제거 수행
                    filtered = result.passed
                    removed = before - len(filtered)

                    if self.verbose and removed > 0:
                        print(f"  Removed {removed:,} masked-region-dominant hits "
                            f"({removed / before * 100:.1f}%) [hit-level masking]")
                else:
                    # use_depth_entropy=True 인 경우에는 hit를 제거하지 않음
                    # masking은 Step4에서 depth/entropy 계산 시 position-level로 반영
                    removed = before - len(result.passed)
                    print(f"  MaskingFilter evaluated hits (no removal, use_depth_entropy=True). "
                        f"[would remove {removed:,} hits if hit-level masking were enabled]")

            except Exception as e:
                if self.verbose:
                    print(f"  [WARNING] MaskingFilter failed: {e}")
                    print(f"  Falling back to target-based masking")

                # fallback: target 단위 제거도 use_depth_entropy=False 일 때만
                if self.masked_targets and not self.use_depth_entropy:
                    filtered2 = filtered[~filtered['target'].isin(self.masked_targets)].copy()
                    removed = len(filtered) - len(filtered2)
                    filtered = filtered2
                    if self.verbose and removed > 0:
                        print(f"  Removed {removed:,} masked target hits "
                            f"({removed / before * 100:.1f}%) [target-level fallback]")

        return filtered

    def step1_filter_hits(self, hits_df: pd.DataFrame) -> pd.DataFrame:
        """
        Step 1: Complete filtering pipeline (calls step1a, step1b, step1c).

        This is a convenience method that runs all filtering steps in sequence:
        1. Quality filtering (hit_filter)
        2. Contamination filtering (NC_032111.1)
        3. Masking filtering (if enabled)

        For fine-grained control, use individual methods:
        - step1a_quality_filter()
        - step1b_contamination_filter()
        - step1c_masking_filter()

        Args:
            hits_df: Input hits DataFrame

        Returns:
            Fully filtered hits DataFrame
        """
        if self.verbose:
            print("Step 1: Complete hit filtering pipeline...")
            print("-" * 60)

        # 1a) Quality filter
        filtered = self.step1a_quality_filter(hits_df)

        # 1b) Contamination filter
        filtered = self.step1b_contamination_filter(filtered)

        # 1c) Masking filter
        filtered = self.step1c_masking_filter(filtered)

        if self.verbose:
            print("-" * 60)
            total_removed = len(hits_df) - len(filtered)
            pct_removed = total_removed / len(hits_df) * 100 if len(hits_df) > 0 else 0
            print(f"Total: {len(filtered):,} / {len(hits_df):,} hits passed "
                  f"({total_removed:,} removed, {pct_removed:.1f}%)")

        return filtered


    # =========================================================================
    # Step2
    # =========================================================================
    def _rank_map_for(self, taxids, rank: str) -> dict:
        """{strain_taxid -> taxid_at_rank} computed once over the UNIQUE taxids
        (hit tables repeat the same strain taxid thousands of times, so a per-row
        apply re-walked the lineage redundantly — this caches per unique taxid)."""
        uniq = pd.unique(pd.Series(taxids))
        return {int(t): self.tax.get_taxid_at_rank(int(t), rank) for t in uniq}

    def step2_convert_to_high_lvl_taxon(self, hits_df: pd.DataFrame):
        if self.segment_lookup_df is None or self.segment_lookup_df.empty:
            hits_df['strain_taxid'] = hits_df['taxid'].fillna(0).astype(int)
            st = hits_df['strain_taxid']
            hits_df['species_taxid'] = st.map(self._rank_map_for(st, "species"))
            hits_df['genus_taxid'] = st.map(self._rank_map_for(st, "genus"))
            hits_df['family_taxid'] = st.map(self._rank_map_for(st, "family"))
        else:
            merged = hits_df.merge(
                self.segment_lookup_df,
                left_on='target',
                right_index=True,
                how='left'
            )
            # fallback
            raw = merged['taxid'].fillna(0).astype(int)
            merged['strain_taxid'] = merged['ref_strain_taxid'].fillna(raw).astype(int)
            merged['species_taxid'] = merged['ref_species_taxid'].fillna(0)
            merged['genus_taxid'] = merged['ref_genus_taxid']
            merged['family_taxid'] = merged['ref_family_taxid']

            # fill missing — cache rank lookups over the unique strain taxids that
            # still need filling, then map (avoids per-row lineage re-walk).
            mask = merged['species_taxid'] == 0
            if mask.any():
                merged.loc[mask, 'species_taxid'] = merged.loc[mask, 'strain_taxid'].map(
                    self._rank_map_for(merged.loc[mask, 'strain_taxid'], "species"))
            mask = merged['genus_taxid'].isna()
            if mask.any():
                merged.loc[mask, 'genus_taxid'] = merged.loc[mask, 'strain_taxid'].map(
                    self._rank_map_for(merged.loc[mask, 'strain_taxid'], "genus"))
            mask = merged['family_taxid'].isna()
            if mask.any():
                merged.loc[mask, 'family_taxid'] = merged.loc[mask, 'strain_taxid'].map(
                    self._rank_map_for(merged.loc[mask, 'strain_taxid'], "family"))

            hits_df = merged

        rank_map = {
            'strain': 'strain_taxid',
            'species': 'species_taxid',
            'genus': 'genus_taxid',
            'family': 'family_taxid'
        }
        hits_df['working_taxid'] = hits_df[rank_map[self.classification_rank]].fillna(0).astype(int)

        hits_df['virus_id'] = hits_df['strain_taxid'].astype(str)

        # Handle segment_id and segment_length (may not exist if no segment_info)
        if 'ref_segment_id' in hits_df.columns:
            hits_df['segment_id'] = hits_df['ref_segment_id'].fillna(hits_df['target'])
        else:
            hits_df['segment_id'] = hits_df['target']

        if 'ref_segment_length' in hits_df.columns:
            hits_df['segment_length'] = hits_df['ref_segment_length'].fillna(
                hits_df.get('tlen', 0)
            ).astype(int)
        else:
            hits_df['segment_length'] = hits_df.get('tlen', 0).astype(int)

        return hits_df[hits_df['working_taxid'] > 0].copy()

    # =========================================================================
    # Step3
    # =========================================================================
    def step3_create_candidates_vectorized(self, hits_df):
        if 'bits' not in hits_df:
            hits_df['bits'] = hits_df['fident']

        df = (
            hits_df.groupby(['query', 'working_taxid'], as_index=False)
            .agg({'bits': 'max', 'fident': 'max'})
            .rename(columns={'bits': 'best_score', 'fident': 'best_identity'})
        )
        df['candidate_count'] = df.groupby('query')['working_taxid'].transform('count')

        unique_set = set(df[df['candidate_count'] == 1]['query'])
        return df, unique_set

    # =========================================================================
    # Step3b
    # =========================================================================
    def step3b_resolve_multi_mapping(self, hits_df, candidates_df, unique_query_set):

        if self.multi_mapping_mode == "all":
            return hits_df

        # Multi-mapping read list
        all_q = set(candidates_df['query'].unique())
        multi_q = all_q - unique_query_set

        if self.multi_mapping_mode == "best_hit":
            # assign best
            multi_c = candidates_df[candidates_df['query'].isin(multi_q)]
            best = (
                multi_c
                .sort_values(['query', 'best_score', 'best_identity'], ascending=[True, False, False])
                .drop_duplicates('query')
            )
            q2best = dict(zip(best['query'], best['working_taxid']))

            unique_hits = hits_df[hits_df['query'].isin(unique_query_set)]
            multi_hits = hits_df[hits_df['query'].isin(multi_q)].copy()
            multi_hits['best_taxon'] = multi_hits['query'].map(q2best)
            multi_hits = multi_hits[multi_hits['working_taxid'] == multi_hits['best_taxon']]
            multi_hits = multi_hits.drop(columns=['best_taxon'])

            return pd.concat([unique_hits, multi_hits], ignore_index=True)

        elif self.multi_mapping_mode == "local_depth":
            # Use local depth-based assignment
            # This is handled in step4 with iteration
            return hits_df

        return hits_df

    # =========================================================================
    # Step3c: Local Depth Assignment
    # =========================================================================
    def step3c_local_depth_assignment(
        self,
        hits_df: pd.DataFrame,
        candidates_df: pd.DataFrame,
        unique_query_set: Set[str],
        current_depth_profiles: Dict[tuple, np.ndarray]
    ) -> pd.DataFrame:
        """
        Assign multi-mapping reads based on local depth profiles.

        Args:
            hits_df: All alignment hits
            candidates_df: Candidate summary (query, working_taxid, best_score, etc.)
            unique_query_set: Set of unique-mapping queries
            current_depth_profiles: Dict[(virus_id, segment_id)] -> depth array

        Returns:
            hits_df with assignment_weight column added
        """

        if self.verbose:
            print("  Local depth assignment...")

        # Multi-mapping reads
        all_q = set(candidates_df['query'].unique())
        multi_q = all_q - unique_query_set

        if len(multi_q) == 0:
            # No multi-mapping reads
            hits_df['assignment_weight'] = 1.0
            return hits_df

        # Initialize weights
        hits_df['assignment_weight'] = 1.0

        # For each multi-mapping read
        multi_hits = hits_df[hits_df['query'].isin(multi_q)].copy()

        for query, query_hits in multi_hits.groupby('query'):
            # Get candidates for this read
            candidates = query_hits[['working_taxid', 'virus_id', 'segment_id',
                                     'tstart', 'tend', 'bits', 'fident']].to_dict('records')

            if len(candidates) <= 1:
                continue

            # Calculate scores for each candidate
            scores = []

            for cand in candidates:
                virus_id = cand['virus_id']
                segment_id = cand['segment_id']
                tstart = int(cand['tstart'])
                tend = int(cand['tend'])
                bits = cand['bits']
                fident = cand['fident']

                # Alignment score (normalized)
                aln_score = bits / 100.0  # Simple normalization

                # Depth score
                depth_key = (virus_id, segment_id)
                if depth_key in current_depth_profiles:
                    depth_array = current_depth_profiles[depth_key]
                    seg_len = len(depth_array)

                    # Local depth in alignment region
                    tstart_clip = max(0, min(tstart, seg_len))
                    tend_clip = max(0, min(tend, seg_len))

                    if tend_clip > tstart_clip:
                        local_depth = depth_array[tstart_clip:tend_clip].mean()
                        global_depth = depth_array.mean()

                        # Score based on how close local depth is to global
                        # Higher score if local depth is close to global (smooth coverage)
                        if global_depth > 0:
                            depth_score = -abs(local_depth - global_depth) / (global_depth + 1e-6)
                        else:
                            depth_score = 0.0
                    else:
                        depth_score = 0.0
                else:
                    depth_score = 0.0

                # Combined score
                combined_score = (
                    self.local_depth_alpha * aln_score +
                    self.local_depth_beta * depth_score
                )

                scores.append({
                    'working_taxid': cand['working_taxid'],
                    'virus_id': virus_id,
                    'segment_id': segment_id,
                    'score': combined_score
                })

            # Softmax to get weights
            if len(scores) > 0:
                score_values = np.array([s['score'] for s in scores])
                # Apply temperature
                exp_scores = np.exp(score_values / self.local_depth_temperature)
                weights = exp_scores / exp_scores.sum()

                # Assign weights back to hits_df
                for i, score_info in enumerate(scores):
                    mask = (
                        (hits_df['query'] == query) &
                        (hits_df['working_taxid'] == score_info['working_taxid']) &
                        (hits_df['virus_id'] == score_info['virus_id']) &
                        (hits_df['segment_id'] == score_info['segment_id'])
                    )
                    hits_df.loc[mask, 'assignment_weight'] = weights[i]

        if self.verbose:
            weighted_reads = len(multi_q)
            print(f"    Assigned weights to {weighted_reads:,} multi-mapping reads")

        return hits_df

    def _get_expected_segment_counts(self) -> Dict[int, int]:
        """strain_taxid -> number of distinct segments (computed once, O(N), cached).

        Replaces the previous O(N^2) inline set-comprehension that re-scanned all
        segment_info entries for every entry (~53s/sample on the profiled run)."""
        if getattr(self, "_expected_segment_counts", None) is not None:
            return self._expected_segment_counts
        counts: Dict[int, set] = {}
        if self.segment_info:
            for info in self.segment_info.values():
                strain = getattr(info, "strain_taxid", info.species_taxid)
                counts.setdefault(strain, set()).add(info.segment_id)
        self._expected_segment_counts = {k: len(v) for k, v in counts.items()}
        return self._expected_segment_counts

    # =========================================================================
    # Step4
    # =========================================================================
    def _build_depth_arrays(self, seg_hits, seg_len, unique_query_set):
        """Per-position depth for one segment, split unique vs multi-mapping.
        Returns (depth_unique, depth_multi, depth_total) float32 arrays of length
        seg_len. Vectorised difference-array fill — identical to depth[s:e]+=w per
        hit. Extracted from step4 for readability (no behaviour change)."""
        depth_unique = np.zeros(seg_len, dtype=np.float32)
        depth_multi = np.zeros(seg_len, dtype=np.float32)
        depth_total = np.zeros(seg_len, dtype=np.float32)

        _ts = seg_hits['tstart'].to_numpy()
        _te = seg_hits['tend'].to_numpy()
        if 'assignment_weight' in seg_hits.columns:
            _ws = seg_hits['assignment_weight'].to_numpy()
        else:
            _ws = np.ones(len(seg_hits), dtype=np.float32)
        # Prefer the precomputed boolean column (set once in step4); fall back to
        # the set lookup if a caller passes seg_hits without it.
        if '_is_unique' in seg_hits.columns:
            _is_uniq = seg_hits['_is_unique'].to_numpy()
        else:
            _is_uniq = seg_hits['query'].isin(unique_query_set).to_numpy()

        _s = np.clip(_ts.astype(np.int64), 0, seg_len)
        _e = np.clip(_te.astype(np.int64), 0, seg_len)
        _valid = _e > _s
        if _valid.any():
            sv = _s[_valid]; ev = _e[_valid]
            wv = _ws[_valid].astype(np.float32)
            uq = _is_uniq[_valid]

            def _fill(track, sel):
                if not sel.any():
                    return
                d = np.zeros(seg_len + 1, dtype=np.float32)
                np.add.at(d, sv[sel], wv[sel])
                np.add.at(d, ev[sel], -wv[sel])
                track += np.cumsum(d[:-1])

            _fill(depth_unique, uq)
            _fill(depth_multi, ~uq)
            depth_total += depth_unique + depth_multi
        return depth_unique, depth_multi, depth_total

    def step4_calculate_coverage(
        self,
        hits_df: pd.DataFrame,
        unique_query_set: Set[str],
        candidates_df: pd.DataFrame = None
    ) -> Dict[int, TaxonCoverage]:
        """
        Step 4: coverage 계산 (strain 단위) 후 working taxon 수준으로 집계.

        - use_depth_entropy=False AND mask_filter=None 일 때:
            기존 interval 기반 breadth/depth 계산 (calculate_breadth_fast 사용)
        - use_depth_entropy=True 또는 mask_filter 존재할 때:
            segment 길이만큼 depth 배열을 만들고
            - 모든 hit interval을 depth에 반영
            - mask_filter가 있으면 mask 구간 depth를 0으로 설정
            - breadth = depth > 0 인 위치 수
            - avg_depth = depth.sum() / 유효 길이(마스킹 제외)
            - depth_entropy = normalized entropy
        """

        if self.verbose:
            print("Step 4: Calculating coverage...")

        # expected_segment_counts (strain -> #distinct segments).
        # Computed once and cached: the previous inline version was O(N^2) over
        # segment_info (~19k entries) and dominated runtime (~53s/sample). This is
        # input-independent, so cache it on the instance.
        expected_segment_counts = self._get_expected_segment_counts()

        strain_coverage: Dict[int, StrainCoverage] = {}

        # Precompute per-hit uniqueness ONCE (vectorised). Previously each segment
        # called seg_hits['query'].isin(unique_query_set) on a large Python set,
        # which pandas re-coerced to an object array every call — ~87% of step4
        # (16s of 19s). Doing it once here as a boolean column and slicing it in the
        # segment loop is behaviourally identical and removes that cost.
        if "_is_unique" not in hits_df.columns:
            hits_df = hits_df.copy()
            hits_df["_is_unique"] = hits_df["query"].isin(unique_query_set).to_numpy()

        # Strain 단위 coverage
        for virus_id, virus_hits in hits_df.groupby('virus_id'):
            strain_taxid = int(virus_hits['strain_taxid'].iloc[0])
            segment_objs: Dict[str, SegmentCoverage] = {}
            total_len = 0
            _seg_groups = {}  # seg_id -> seg_hits, reused for masked-count (avoids re-groupby)

            # segment별 coverage
            for seg_id, seg_hits in virus_hits.groupby('segment_id'):
                _seg_groups[seg_id] = seg_hits
                # segment 길이 결정 (기존 코드 로직 그대로 유지)
                unique_lengths = seg_hits['segment_length'].unique()
                if len(unique_lengths) == 1:
                    seg_len = int(unique_lengths[0])
                else:
                    length_counts = seg_hits.groupby('segment_length').size()
                    total_weighted = sum(length * count for length, count in length_counts.items())
                    seg_len = int(total_weighted / len(seg_hits))

                if seg_len <= 0:
                    continue

                # ---------------------------
                # 4A. breadth/depth 계산
                # ---------------------------
                use_depth_array = self.use_depth_entropy or (self.mask_filter is not None)

                if not use_depth_array:
                    # 빠른 path, interval 기반 breadth 계산
                    intervals = seg_hits[['tstart', 'tend']].sort_values('tstart').values
                    breadth_bp = calculate_breadth_fast(intervals)
                    total_aligned = np.sum(intervals[:, 1] - intervals[:, 0])
                    avg_depth = total_aligned / seg_len if seg_len > 0 else 0.0
                    depth_entropy = None
                    unique_breadth_bp = 0
                    multi_breadth_bp = 0
                    masked_breadth_bp = 0

                else:
                    # depth 배열 기반 path (masking-aware, entropy 계산 가능)
                    # 1) Unique와 multi-mapping을 구분하여 depth 계산 (헬퍼로 분리)
                    depth_unique, depth_multi, depth_total = self._build_depth_arrays(
                        seg_hits, seg_len, unique_query_set)

                    # 2) Masking 처리
                    mask_arr = np.zeros(seg_len, dtype=bool)
                    effective_len = seg_len

                    if self.mask_filter:
                        # target 또는 seg_id로 masked region 찾기
                        mask_region = self.mask_filter._mask_dict.get(seg_id, None)
                        if not mask_region and len(seg_hits) > 0:
                            # fallback: target 이름으로 시도
                            first_target = seg_hits['target'].iloc[0]
                            mask_region = self.mask_filter._mask_dict.get(first_target, None)

                        if mask_region:
                            # MaskedRegion 객체의 starts와 ends 배열을 사용
                            for ms, me in zip(mask_region.starts, mask_region.ends):
                                if me <= 0 or ms >= seg_len:
                                    continue
                                ms_c = max(0, int(ms))
                                me_c = min(seg_len, int(me))
                                if me_c <= ms_c:
                                    continue
                                mask_arr[ms_c:me_c] = True

                            # Case 1: use_depth_entropy=False → masked region을 depth에서 제거
                            # Case 2: use_depth_entropy=True → masked region은 표시만 하고 제거하지 않음
                            if not self.use_depth_entropy:
                                depth_unique[mask_arr] = 0
                                depth_multi[mask_arr] = 0
                                depth_total[mask_arr] = 0
                                effective_len = seg_len - int(mask_arr.sum())
                                if effective_len <= 0:
                                    effective_len = seg_len
                            else:
                                # entropy 모드에서는 masked region도 depth에 포함
                                effective_len = seg_len

                    # 3) breadth, depth, entropy 계산
                    covered_unique = depth_unique > 0
                    covered_multi = depth_multi > 0
                    covered_total = depth_total > 0
                    covered_masked = covered_total & mask_arr

                    unique_breadth_bp = int(covered_unique.sum())
                    multi_breadth_bp = int(covered_multi.sum())
                    masked_breadth_bp = int(covered_masked.sum())
                    breadth_bp = int(covered_total.sum())

                    total_aligned = float(depth_total.sum())
                    avg_depth = total_aligned / effective_len if effective_len > 0 else 0.0

                    if self.use_depth_entropy:
                        depth_entropy = compute_depth_entropy(depth_total)
                    else:
                        depth_entropy = None

                # SegmentCoverage 객체 생성
                segment_objs[seg_id] = SegmentCoverage(
                    segment_id=seg_id,
                    length=seg_len,
                    breadth_bp=breadth_bp,
                    breadth_ratio=breadth_bp / seg_len if seg_len > 0 else 0.0,
                    avg_depth=avg_depth,
                    hit_count=len(seg_hits),
                    depth_entropy=depth_entropy,
                    unique_breadth_bp=unique_breadth_bp,
                    multi_breadth_bp=multi_breadth_bp,
                    masked_breadth_bp=masked_breadth_bp,
                )

                total_len += seg_len

            if total_len == 0:
                continue

            # Breadth denominator. Legacy ("aligned"): total length of the segments
            # this strain's reads landed on. "expected": the taxon's whole-genome
            # length from the DB, so covering only the short segment of a
            # multipartite virus does not inflate breadth (e.g. Shamonda S 92bp /
            # 12104bp = 0.8% instead of /927bp = 9.8%). The numerator (covered bp)
            # is unchanged; only the denominator changes, and only in "expected"
            # mode — so the legacy path is byte-identical.
            breadth_denom = total_len
            if self.breadth_denominator == "expected":
                egl = self._expected_genome_length.get(strain_taxid)
                if egl and egl >= total_len:   # guard: never inflate breadth above 1
                    breadth_denom = egl

            # weighted breadth (strain 수준)
            weighted_breadth = sum(
                seg.breadth_ratio * seg.length for seg in segment_objs.values()
            ) / breadth_denom

            # unmasked weighted breadth (excluding masked regions)
            unmasked_weighted_breadth = sum(
                ((seg.breadth_bp - seg.masked_breadth_bp) / seg.length * seg.length)
                for seg in segment_objs.values()
                if seg.length > 0
            ) / breadth_denom

            # unique read count
            read_ids = virus_hits['query'].unique()
            # count distinct unique-mapping reads via the precomputed column when
            # available (avoids a per-read set lookup); identical result.
            if '_is_unique' in virus_hits.columns:
                unique_count = int(virus_hits.loc[virus_hits['_is_unique'], 'query'].nunique())
            else:
                unique_count = sum(1 for q in read_ids if q in unique_query_set)

            # masked read count (reads that overlap with masked regions)
            masked_read_count = 0
            if self.mask_filter:
                masked_queries = set()
                for seg_id, seg_hits in _seg_groups.items():  # reuse first-pass grouping
                    # Get masked regions for this segment
                    mask_region = self.mask_filter._mask_dict.get(seg_id, None)
                    if not mask_region and len(seg_hits) > 0:
                        first_target = seg_hits['target'].iloc[0]
                        mask_region = self.mask_filter._mask_dict.get(first_target, None)

                    if mask_region:
                        # Vectorised overlap test — behaviourally identical to the
                        # previous per-row iterrows() loop (the dominant hot path,
                        # ~19s/sample) but ~100x faster. Uses raw tstart/tend (no
                        # reversed-interval normalisation) to match the old result
                        # exactly: overlap = min(e,me) - max(s,ms) > 0.
                        s_arr = seg_hits['tstart'].to_numpy(dtype=np.int64)
                        e_arr = seg_hits['tend'].to_numpy(dtype=np.int64)
                        queries = seg_hits['query'].to_numpy()
                        overlapped = np.zeros(len(seg_hits), dtype=bool)
                        for ms, me in zip(mask_region.starts, mask_region.ends):
                            overlapped |= (np.minimum(e_arr, me) - np.maximum(s_arr, ms)) > 0
                        masked_queries.update(queries[overlapped].tolist())

                masked_read_count = len(masked_queries)

            # segments_detected
            segments_detected = sum(
                1
                for seg in segment_objs.values()
                if (
                    seg.breadth_ratio >= self.thresholds.min_segment_breadth
                    or seg.breadth_bp >= self.thresholds.min_segment_covered_bp
                )
            )

            expected_count = expected_segment_counts.get(strain_taxid, len(segment_objs))

            strain_coverage[strain_taxid] = StrainCoverage(
                strain_taxid=strain_taxid,
                virus_id=virus_id,
                total_genome_length=total_len,
                segment_count=len(segment_objs),
                expected_segment_count=expected_count,
                segment_coverages=segment_objs,
                weighted_breadth=weighted_breadth,
                unmasked_weighted_breadth=unmasked_weighted_breadth,
                segments_detected=segments_detected,
                unique_read_count=unique_count,
                multi_mapping_read_count=len(read_ids) - unique_count,
                total_read_count=len(read_ids),
                masked_read_count=masked_read_count,
            )

        # ---------------------------
        # 4B. working_taxid 수준으로 집계
        # ---------------------------

        taxon_coverage: Dict[int, TaxonCoverage] = {}

        # strain_taxid -> working_taxid 매핑
        strain_to_working: Dict[int, int] = {}
        for _, row in hits_df[['strain_taxid', 'working_taxid']].drop_duplicates().iterrows():
            strain_to_working[int(row['strain_taxid'])] = int(row['working_taxid'])

        for strain_taxid, strain_cov in strain_coverage.items():
            working_taxid = strain_to_working.get(strain_taxid, strain_taxid)

            if working_taxid not in taxon_coverage:
                # Get taxon name from taxonomy database
                taxon_name = ""
                try:
                    taxon_name = self.tax.get_name(working_taxid)
                except:
                    taxon_name = f"Unknown ({working_taxid})"

                taxon_coverage[working_taxid] = TaxonCoverage(
                    taxon_taxid=working_taxid,
                    taxon_rank=self.classification_rank,
                    strains={},
                    avg_genome_length=0.0,
                    weighted_breadth=0.0,
                    is_real=False,
                    taxon_name=taxon_name,
                )

            taxon_coverage[working_taxid].strains[strain_taxid] = strain_cov

        # taxon 수준 평균 genome length와 breadth 계산
        for taxon_taxid, taxon_cov in taxon_coverage.items():
            if not taxon_cov.strains:
                continue

            total_length = sum(s.total_genome_length for s in taxon_cov.strains.values())
            taxon_cov.avg_genome_length = total_length / len(taxon_cov.strains)

            taxon_cov.weighted_breadth = (
                sum(s.weighted_breadth for s in taxon_cov.strains.values())
                / len(taxon_cov.strains)
            )

            taxon_cov.unmasked_weighted_breadth = (
                sum(s.unmasked_weighted_breadth for s in taxon_cov.strains.values())
                / len(taxon_cov.strains)
            )

        if self.verbose:
            print(f"  Calculated coverage for {len(strain_coverage):,} strains")
            print(f"  Aggregated to {len(taxon_coverage):,} {self.classification_rank}-level taxa")

        return taxon_coverage


    # =========================================================================
    # Step5
    # =========================================================================
    def step5_judge_real_vs_fake(self, taxon_cov):

        real_taxa = set()
        fake_taxa = set()

        for tid, t in taxon_cov.items():
            total_uniq = sum(s.unique_read_count for s in t.strains.values())
            total_multi = sum(s.multi_mapping_read_count for s in t.strains.values())

            # ============================================================
            # SIMPLE MODE: All criteria must be met (AND logic)
            # ============================================================
            if self.simple_mode:
                simple_criteria = total_uniq > total_multi

                # Alternative criteria: High read count + unmasked breadth
                total_reads = total_uniq + total_multi
                high_count_ok = (
                    total_reads >= self.thresholds.min_read_count and
                    t.unmasked_weighted_breadth >= self.thresholds.min_unmasked_weighted_breadth
                )

                # All conditions must be satisfied (AND logic)
                is_real = simple_criteria and high_count_ok

                # Verbose logging
                if self.verbose and is_real:
                    print(f"  [SIMPLE] {tid} ({t.taxon_name}): "
                          f"unique={total_uniq}, multi={total_multi}, "
                          f"reads={total_reads}, unmasked_breadth={t.unmasked_weighted_breadth:.4f} → REAL")

                t.is_real = is_real
                for s in t.strains.values():
                    s.is_real = is_real

                if is_real:
                    real_taxa.add(tid)
                else:
                    fake_taxa.add(tid)

                continue
            # ============================================================

            # segment check (species 예외 포함)
            segment_ok = False
            for s in t.strains.values():
                req_seg = max(1, int(s.expected_segment_count * self.thresholds.min_segment_fraction))
                if s.segments_detected >= req_seg:
                    segment_ok = True
                    break

            breadth_ok = t.weighted_breadth >= self.thresholds.min_weighted_breadth
            unique_ok = (not self.thresholds.require_unique_reads) or total_uniq > 0

            # -------- Depth entropy check --------
            entropy_ok = True
            if self.use_depth_entropy:
                entropies = []
                for s in t.strains.values():
                    for seg in s.segment_coverages.values():
                        if seg.depth_entropy is not None:
                            entropies.append(seg.depth_entropy)
                best_entropy = max(entropies) if entropies else 0.0
                entropy_ok = best_entropy >= self.thresholds.min_depth_entropy
            # -------------------------------------

            # -------- Alternative criteria: High read count + unmasked breadth --------
            total_reads = total_uniq + total_multi
            high_count_ok = (
                total_reads >= self.thresholds.min_read_count and
                t.unmasked_weighted_breadth >= self.thresholds.min_unmasked_weighted_breadth
            )
            # --------------------------------------------------------------------------

            # --- FP-leak fix: unique-fraction gate (off when threshold <= 0) ---
            unique_fraction_ok = True
            if self.thresholds.min_unique_fraction > 0:
                ufrac = (total_uniq / total_reads) if total_reads > 0 else 0.0
                unique_fraction_ok = ufrac >= self.thresholds.min_unique_fraction

            # All criteria must be satisfied (AND logic)
            standard_ok = breadth_ok and segment_ok and unique_ok and entropy_ok
            is_real = standard_ok and high_count_ok and unique_fraction_ok

            # Verbose logging
            if self.verbose and is_real:
                print(f"  [NORMAL] {tid} ({t.taxon_name}): "
                      f"breadth={t.weighted_breadth:.4f}, segments={segment_ok}, "
                      f"reads={total_reads}, unmasked_breadth={t.unmasked_weighted_breadth:.4f} → REAL")

            t.is_real = is_real
            for s in t.strains.values():
                s.is_real = is_real

            if is_real:
                real_taxa.add(tid)
            else:
                fake_taxa.add(tid)

        # genus-competition removed — superseded by the relative-abundance cut
        # (apply_fp_postfilter, applied uniformly to lca/em/coverage in the CLI).
        return real_taxa, fake_taxa

    def _genus_of(self, taxid):
        """Genus taxid for a taxon (None if unavailable)."""
        try:
            return self.tax.get_taxid_at_rank(int(taxid), "genus")
        except Exception:
            return None

    def _apply_genus_competition(self, taxon_cov, real_taxa, fake_taxa):
        """Among same-genus REAL taxa, demote relatives whose unique support is a
        tiny fraction of the genus winner's — these are typically cross-mapping
        spillover from the dominant species (e.g. Cytomegalovirus relatives)."""
        # unique-read support per real taxon
        uniq = {}
        for tid in real_taxa:
            t = taxon_cov.get(tid)
            uniq[tid] = sum(s.unique_read_count for s in t.strains.values()) if t else 0
        # group by genus
        by_genus = {}
        for tid in real_taxa:
            g = self._genus_of(tid)
            if g is None:
                continue
            by_genus.setdefault(g, []).append(tid)
        ratio = self.thresholds.genus_competition_ratio
        demoted = set()
        for g, members in by_genus.items():
            if len(members) < 2:
                continue
            winner_u = max(uniq[m] for m in members)
            if winner_u <= 0:
                continue
            for m in members:
                if uniq[m] < ratio * winner_u:
                    demoted.add(m)
                    if self.verbose:
                        print(f"  [GENUS-COMP] demote {m} ({taxon_cov[m].taxon_name}): "
                              f"unique {uniq[m]} < {ratio}×winner({winner_u})")
        if demoted:
            for tid in demoted:
                taxon_cov[tid].is_real = False
                for s in taxon_cov[tid].strains.values():
                    s.is_real = False
            real_taxa = real_taxa - demoted
            fake_taxa = fake_taxa | demoted
        return real_taxa, fake_taxa

    # =========================================================================
    # Step6
    # =========================================================================
    def step6_assign_reads_vectorized(self, candidates_df, real_taxa):
        """
        Step 6: Assign reads to real taxa.
        Returns a DataFrame with query -> taxon assignments.
        """
        if self.verbose:
            print("Step 6: Assigning reads to real taxa...")

        # Filter candidates to only real taxa
        real_candidates = candidates_df[candidates_df['working_taxid'].isin(real_taxa)].copy()

        # For unique reads, direct assignment
        unique_reads = real_candidates[real_candidates['candidate_count'] == 1].copy()
        unique_reads['assignment_type'] = 'unique'

        # For multi-mapping reads, assign to best hit among real taxa
        multi_reads = real_candidates[real_candidates['candidate_count'] > 1].copy()
        if len(multi_reads) > 0:
            best_multi = (
                multi_reads
                .sort_values(['query', 'best_score', 'best_identity'], ascending=[True, False, False])
                .drop_duplicates('query', keep='first')
            )
            best_multi['assignment_type'] = 'multi'
        else:
            best_multi = pd.DataFrame(columns=unique_reads.columns)

        # Combine
        assignment = pd.concat([unique_reads, best_multi], ignore_index=True)

        if self.verbose:
            print(f"  Assigned {len(assignment):,} reads to {len(real_taxa)} real taxa")
            print(f"    Unique: {len(unique_reads):,}, Multi: {len(best_multi):,}")

        return assignment

    # =========================================================================
    # Step7
    # =========================================================================
    def step7_calculate_abundance(self, assignment_df, taxon_cov):
        """
        Step 7: Calculate abundance (read counts per taxon).

        Returns a DataFrame with:
            - taxon_taxid, taxon_name
            - read_count: assigned reads (unique + multi-mapping)
            - normalized_abundance: read_count / total_reads
            - length_norm_abundance: (read_count / genome_length) / sum(read_count / genome_length)
            - genome_length: average genome length
            - RPK: Reads Per Kilobase
            - TPM: Transcripts Per Million
            - weighted_breadth, is_real
        """
        if self.verbose:
            print("Step 7: Calculating abundance...")

        if len(assignment_df) == 0:
            return pd.DataFrame(columns=[
                'taxon_taxid', 'taxon_name', 'read_count', 'normalized_abundance',
                'length_norm_abundance', 'genome_length', 'RPK', 'TPM',
                'weighted_breadth', 'is_real'
            ])

        # Count reads per taxon
        abundance = (
            assignment_df
            .groupby('working_taxid')
            .size()
            .reset_index(name='read_count')
            .rename(columns={'working_taxid': 'taxon_taxid'})
        )

        # Total assigned reads for normalization
        total_reads = abundance['read_count'].sum()

        # Add taxon info from coverage
        abundance['taxon_name'] = abundance['taxon_taxid'].map(
            lambda tid: taxon_cov[tid].taxon_name if tid in taxon_cov else f"Unknown ({tid})"
        )
        abundance['weighted_breadth'] = abundance['taxon_taxid'].map(
            lambda tid: taxon_cov[tid].weighted_breadth if tid in taxon_cov else 0.0
        )
        abundance['is_real'] = abundance['taxon_taxid'].map(
            lambda tid: taxon_cov[tid].is_real if tid in taxon_cov else False
        )

        # Get genome length (average across strains)
        abundance['genome_length'] = abundance['taxon_taxid'].map(
            lambda tid: taxon_cov[tid].avg_genome_length if tid in taxon_cov else 0.0
        )

        # Calculate normalized abundance (relative abundance)
        abundance['normalized_abundance'] = abundance['read_count'] / total_reads if total_reads > 0 else 0.0

        # Calculate length-normalized abundance (genome length corrected relative abundance)
        # This corrects for genome size bias: larger genomes don't get artificially inflated counts
        abundance['length_normalized_count'] = abundance.apply(
            lambda row: row['read_count'] / row['genome_length'] if row['genome_length'] > 0 else 0.0,
            axis=1
        )
        total_length_norm = abundance['length_normalized_count'].sum()
        abundance['length_norm_abundance'] = (
            abundance['length_normalized_count'] / total_length_norm
            if total_length_norm > 0 else 0.0
        )

        # Calculate TPM (Transcripts Per Million)
        # TPM = (read_count / (genome_length / 1000)) / sum(RPK) * 1,000,000
        abundance['RPK'] = abundance.apply(
            lambda row: row['read_count'] / (row['genome_length'] / 1000) if row['genome_length'] > 0 else 0.0,
            axis=1
        )

        total_rpk = abundance['RPK'].sum()
        abundance['TPM'] = abundance['RPK'] / total_rpk * 1_000_000 if total_rpk > 0 else 0.0

        # Drop intermediate column
        abundance = abundance.drop(columns=['length_normalized_count'])

        # Reorder columns
        abundance = abundance[[
            'taxon_taxid', 'taxon_name', 'read_count', 'normalized_abundance',
            'length_norm_abundance', 'genome_length', 'RPK', 'TPM',
            'weighted_breadth', 'is_real'
        ]]
        abundance = abundance.sort_values('read_count', ascending=False)

        if self.verbose:
            print(f"  Calculated abundance for {len(abundance)} taxa")
            print(f"  Total assigned reads: {total_reads:,}")
            print(f"  Normalized abundance sum: {abundance['normalized_abundance'].sum():.4f}")
            print(f"  Length-normalized abundance sum: {abundance['length_norm_abundance'].sum():.4f}")
            print(f"  TPM sum: {abundance['TPM'].sum():.0f}")

        return abundance


    # =========================================================================
    # Helper: Summary
    # =========================================================================
    def get_summary_dataframe(self, taxon_cov):
        """
        Create a summary DataFrame with all detected taxa and their metrics.

        Returns:
            DataFrame with columns: taxon_taxid, taxon_name, is_real, weighted_breadth,
                                    unique_reads, multi_reads, total_reads, masked_reads,
                                    segments_detected, depth_entropy metrics
        """
        rows = []

        for taxid, tcov in taxon_cov.items():
            # Aggregate strain-level metrics
            total_unique = sum(s.unique_read_count for s in tcov.strains.values())
            total_multi = sum(s.multi_mapping_read_count for s in tcov.strains.values())
            total_reads = sum(s.total_read_count for s in tcov.strains.values())
            total_masked = sum(s.masked_read_count for s in tcov.strains.values())
            total_segments_detected = sum(s.segments_detected for s in tcov.strains.values())

            # Aggregate unique/multi breadth
            total_unique_breadth = sum(
                sum(seg.unique_breadth_bp for seg in s.segment_coverages.values())
                for s in tcov.strains.values()
            )
            total_multi_breadth = sum(
                sum(seg.multi_breadth_bp for seg in s.segment_coverages.values())
                for s in tcov.strains.values()
            )
            total_masked_breadth = sum(
                sum(seg.masked_breadth_bp for seg in s.segment_coverages.values())
                for s in tcov.strains.values()
            )

            # Collect entropy values
            entropies = []
            for s in tcov.strains.values():
                for seg in s.segment_coverages.values():
                    if seg.depth_entropy is not None:
                        entropies.append(seg.depth_entropy)

            # Calculate entropy statistics
            best_entropy = max(entropies) if entropies else None
            avg_entropy = sum(entropies) / len(entropies) if entropies else None
            min_entropy = min(entropies) if entropies else None

            rows.append({
                'taxon_taxid': taxid,
                'taxon_name': tcov.taxon_name,
                'taxon_rank': tcov.taxon_rank,
                'is_real': tcov.is_real,
                'weighted_breadth': tcov.weighted_breadth,
                'unmasked_weighted_breadth': tcov.unmasked_weighted_breadth,
                'avg_genome_length': tcov.avg_genome_length,
                'unique_reads': total_unique,
                'multi_reads': total_multi,
                'total_reads': total_reads,
                'masked_reads': total_masked,
                'unique_breadth_bp': total_unique_breadth,
                'multi_breadth_bp': total_multi_breadth,
                'masked_breadth_bp': total_masked_breadth,
                'segments_detected': total_segments_detected,
                'strain_count': len(tcov.strains),
                'best_entropy': best_entropy,
                'avg_entropy': avg_entropy,
                'min_entropy': min_entropy
            })

        df = pd.DataFrame(rows)
        if len(df) > 0:
            df = df.sort_values('total_reads', ascending=False)

        return df

    # =========================================================================
    # Main
    # =========================================================================
    def classify(self, hits_df):

        filtered = self.step1_filter_hits(hits_df)
        taxon_hits = self.step2_convert_to_high_lvl_taxon(filtered)
        candidates, uniq_set = self.step3_create_candidates_vectorized(taxon_hits)
        resolved = self.step3b_resolve_multi_mapping(taxon_hits, candidates, uniq_set)
        taxon_cov = self.step4_calculate_coverage(resolved, uniq_set)
        real_taxa, fake_taxa = self.step5_judge_real_vs_fake(taxon_cov)

        # step6,7 your original code
        assignment = self.step6_assign_reads_vectorized(candidates, real_taxa)
        abundance = self.step7_calculate_abundance(assignment, taxon_cov)

        return assignment, abundance, taxon_cov
