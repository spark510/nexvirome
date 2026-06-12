"""
EM (Expectation-Maximization) Classifier for multi-mapping read resolution.

Algorithm:
    1. Parse & filter hits (quality + masking)
    2. Build read-to-species mapping (which reads map to which species)
    3. Initialize: uniform species abundance
    4. E-step: For each multi-mapping read, assign fractional counts
              proportional to current species abundance × alignment score
    5. M-step: Re-estimate species abundance from fractional counts,
              normalized by genome length (RPKM-like)
    6. Repeat E/M until convergence or max iterations
    7. Final assignment: each read → species with highest posterior probability
    8. Output: LCA-compatible DataFrame + abundance table

This resolves the key limitation of simple LCA: instead of pushing
multi-mapping reads up to higher taxonomic ranks, EM distributes them
to the most likely species based on global abundance patterns.
"""
from __future__ import annotations

from typing import Dict, Optional, Set
from collections import defaultdict
import numpy as np
import pandas as pd

from ..core import FilterCriteria, log_info, log_verbose
from ..alignment import AlignmentParser, MaskingFilter


class EMClassifier:
    """
    EM-based virome classifier.

    Iteratively estimates species abundance and resolves multi-mapping reads
    by distributing them proportionally to estimated abundance × alignment quality.
    """

    def __init__(
        self,
        taxonomy_db,
        segment_info=None,
        min_identity: float = 0.8,
        min_length: int = 50,
        max_evalue: float = 1e-3,
        min_query_coverage: float = 0.5,
        mask_bed_file: Optional[str] = None,
        min_unmasked_coverage: float = 0.1,
        max_iterations: int = 20,
        convergence_threshold: float = 1e-6,
        verbose: bool = False,
    ):
        self.tax = taxonomy_db
        self.segment_info = segment_info or {}
        self.filter_criteria = FilterCriteria(
            min_identity=min_identity,
            min_alignment_length=min_length,
            max_evalue=max_evalue,
            min_query_coverage=min_query_coverage,
        )
        self.mask_bed_file = mask_bed_file
        self.min_unmasked_coverage = min_unmasked_coverage
        self.max_iterations = max_iterations
        self.convergence_threshold = convergence_threshold
        self.verbose = verbose

    def classify(self, hits_df: pd.DataFrame):
        """
        Run the full EM classification pipeline.

        Args:
            hits_df: Combined R1+R2 alignment hits DataFrame

        Returns:
            (lca_df, abundance_df): LCA-compatible results and abundance table
        """
        # Step 1: Quality filter
        log_info("\n🔍 EM Step 1: Quality filtering...")
        parser = AlignmentParser()
        filtered = parser.filter(hits_df, self.filter_criteria)
        log_info(f"  {len(filtered):,}/{len(hits_df):,} hits passed quality filter")

        # Step 2: Masking filter
        if self.mask_bed_file:
            log_info("🎭 EM Step 2: Masking filter...")
            mask_filter = MaskingFilter.from_bed_file(self.mask_bed_file)
            result = mask_filter.filter_by_unmasked_coverage(
                filtered, min_coverage=self.min_unmasked_coverage
            )
            filtered = result.passed
            log_info(f"  {len(filtered):,} hits after masking")

        # Step 3: Resolve species-level taxid using segment_info
        log_info("🧬 EM Step 3: Resolving species-level taxonomy...")
        filtered = self._resolve_species_taxid(filtered)

        # Step 4: Build read-species mapping
        log_info("📊 EM Step 4: Building read-species mapping...")
        read_species, species_genome_len = self._build_read_species_map(filtered)

        unique_count = sum(1 for spp in read_species.values() if len(spp) == 1)
        multi_count = len(read_species) - unique_count
        n_species = len(species_genome_len)
        log_info(f"  {len(read_species):,} reads → {n_species} species")
        log_info(f"  Unique: {unique_count:,}, Multi-mapping: {multi_count:,}")

        if n_species == 0:
            log_info("⚠️  No species found, returning empty results")
            return pd.DataFrame(), pd.DataFrame()

        # Step 5: EM iterations
        log_info(f"🔄 EM Step 5: Running EM ({self.max_iterations} max iterations)...")
        species_list = sorted(species_genome_len.keys())
        sp_idx = {sp: i for i, sp in enumerate(species_list)}
        n_sp = len(species_list)

        # Genome lengths for normalization
        genome_lengths = np.array([species_genome_len[sp] for sp in species_list], dtype=np.float64)
        # Avoid div by zero; minimum 1kb
        genome_lengths = np.maximum(genome_lengths, 1000.0)

        # Initialize: uniform abundance (genome-length-corrected)
        abundance = np.ones(n_sp, dtype=np.float64) / n_sp

        # Pre-build read data structures for vectorized EM
        # For each read: list of (species_index, alignment_score)
        read_data = []
        for query, species_scores in read_species.items():
            entries = []
            for sp_taxid, score in species_scores.items():
                if sp_taxid in sp_idx:
                    entries.append((sp_idx[sp_taxid], score))
            if entries:
                read_data.append(entries)

        for iteration in range(self.max_iterations):
            # E-step: compute fractional assignments
            new_counts = np.zeros(n_sp, dtype=np.float64)

            for entries in read_data:
                if len(entries) == 1:
                    # Unique mapping — full count
                    new_counts[entries[0][0]] += 1.0
                else:
                    # Multi-mapping — distribute by abundance × score
                    indices = np.array([e[0] for e in entries])
                    scores = np.array([e[1] for e in entries], dtype=np.float64)

                    # Posterior: P(species|read) ∝ abundance[species] × score
                    posterior = abundance[indices] * scores
                    total = posterior.sum()
                    if total > 0:
                        posterior /= total
                    else:
                        posterior = np.ones_like(posterior) / len(posterior)

                    new_counts[indices] += posterior

            # M-step: re-estimate abundance (genome-length-corrected)
            # RPKM-like: abundance ∝ count / genome_length
            rpk = new_counts / genome_lengths
            total_rpk = rpk.sum()
            new_abundance = rpk / total_rpk if total_rpk > 0 else np.ones(n_sp) / n_sp

            # Check convergence
            delta = np.abs(new_abundance - abundance).max()
            abundance = new_abundance

            if self.verbose:
                log_verbose(f"  Iteration {iteration + 1}: max_delta={delta:.2e}")

            if delta < self.convergence_threshold:
                log_info(f"  Converged at iteration {iteration + 1} (delta={delta:.2e})")
                break
        else:
            log_info(f"  Reached max iterations ({self.max_iterations}), delta={delta:.2e}")

        # Step 6: Final read assignment (MAP — maximum a posteriori)
        log_info("📝 EM Step 6: Final read assignment...")
        assignments = {}  # query -> (species_taxid, posterior_prob)

        for query, species_scores in read_species.items():
            entries = []
            for sp_taxid, score in species_scores.items():
                if sp_taxid in sp_idx:
                    entries.append((sp_taxid, sp_idx[sp_taxid], score))

            if not entries:
                continue

            if len(entries) == 1:
                assignments[query] = (entries[0][0], 1.0)
            else:
                indices = np.array([e[1] for e in entries])
                scores = np.array([e[2] for e in entries], dtype=np.float64)
                posterior = abundance[indices] * scores
                total = posterior.sum()
                if total > 0:
                    posterior /= total
                best_idx = np.argmax(posterior)
                assignments[query] = (entries[best_idx][0], float(posterior[best_idx]))

        # Step 7: Build output DataFrames
        log_info("📊 EM Step 7: Building output tables...")

        # LCA-compatible format
        lca_rows = []
        for query, (taxid, confidence) in assignments.items():
            lca_rows.append({
                'query': query,
                'lca_taxid': taxid,
                'lca_name': self.tax.get_name(taxid) or f"Unknown ({taxid})",
                'lca_rank': self.tax.get_rank(taxid) or 'no rank',
                'qlen': 100,
                'read_count': 1,
                'n_hits': len(read_species[query]),
                'n_unique_taxids': len(read_species[query]),
                'all_taxids': ','.join(str(t) for t in read_species[query].keys()),
            })

        lca_df = pd.DataFrame(lca_rows)

        # Abundance table
        abundance_rows = []
        # Final counts from assignments
        final_counts = defaultdict(int)
        for taxid, _ in assignments.values():
            final_counts[taxid] += 1

        total_assigned = sum(final_counts.values())

        for i, sp_taxid in enumerate(species_list):
            count = final_counts.get(sp_taxid, 0)
            if count == 0:
                continue

            abundance_rows.append({
                'taxon_taxid': sp_taxid,
                'taxon_name': self.tax.get_name(sp_taxid) or f"Unknown ({sp_taxid})",
                'read_count': count,
                'normalized_abundance': count / total_assigned if total_assigned > 0 else 0,
                'genome_length': int(genome_lengths[i]),
                'RPK': count / (genome_lengths[i] / 1000),
                'em_abundance': float(abundance[i]),
            })

        abundance_df = pd.DataFrame(abundance_rows)
        if len(abundance_df) > 0:
            total_rpk = abundance_df['RPK'].sum()
            abundance_df['TPM'] = abundance_df['RPK'] / total_rpk * 1e6 if total_rpk > 0 else 0
            abundance_df = abundance_df.sort_values('read_count', ascending=False)

        log_info(f"✅ EM classification complete: {len(assignments):,} reads → "
                 f"{len(abundance_df)} species")

        return lca_df, abundance_df

    def _resolve_species_taxid(self, hits_df: pd.DataFrame) -> pd.DataFrame:
        """Add species_taxid column using segment_info or taxonomy lookup."""
        if not self.segment_info:
            # No segment info — resolve from taxid column
            hits_df = hits_df.copy()
            hits_df['species_taxid'] = hits_df['taxid'].apply(
                lambda x: self.tax.get_taxid_at_rank(int(x), "species") if pd.notna(x) and int(x) > 0 else 0
            )
            return hits_df[hits_df['species_taxid'] > 0]

        # Use segment_info for lookup
        hits_df = hits_df.copy()

        def _get_species(row):
            target = row['target']
            seg = self.segment_info.get(target)
            if seg:
                return seg.species_taxid
            taxid = row.get('taxid', 0)
            if pd.notna(taxid) and int(taxid) > 0:
                return self.tax.get_taxid_at_rank(int(taxid), "species")
            return 0

        hits_df['species_taxid'] = hits_df.apply(_get_species, axis=1)
        return hits_df[hits_df['species_taxid'] > 0]

    def _build_read_species_map(self, hits_df: pd.DataFrame):
        """
        Build mapping: query -> {species_taxid: best_alignment_score}
        Also collect genome lengths per species.
        """
        read_species: Dict[str, Dict[int, float]] = defaultdict(dict)
        species_genome_len: Dict[int, int] = {}

        for _, row in hits_df.iterrows():
            query = row['query']
            sp_taxid = int(row['species_taxid'])
            score = float(row.get('bits', row.get('fident', 1.0)))

            # Keep best score per species per read
            if sp_taxid not in read_species[query] or score > read_species[query][sp_taxid]:
                read_species[query][sp_taxid] = score

            # Collect genome length (use tlen or segment_info)
            if sp_taxid not in species_genome_len:
                seg = self.segment_info.get(row.get('target', ''))
                if seg:
                    # Sum all segment lengths for this species
                    total_len = sum(
                        s.segment_length for s in self.segment_info.values()
                        if s.species_taxid == sp_taxid
                    )
                    species_genome_len[sp_taxid] = total_len if total_len > 0 else int(row.get('tlen', 1000))
                else:
                    species_genome_len[sp_taxid] = int(row.get('tlen', 1000))

        return dict(read_species), species_genome_len
