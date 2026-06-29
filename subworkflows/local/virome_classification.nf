/*
 * Virome Classification Subworkflow
 *
 * 1. Run MMseqs2 easy-search on R1 and R2 separately
 * 2. Run virome_classifier for LCA classification
 * 3. Merge all samples into OTU table
 */

include { MMSEQS_CREATEDB          } from '../../modules/nf-core/mmseqs/createdb/main'
include { MMSEQS_EASYSEARCH        } from '../../modules/nf-core/mmseqs/easysearch/main'
// include { VIROME_CLASSIFY       } from '../../modules/local/virome_classify/main'  // DISABLED: LCA mode retired
include { VIROME_CLASSIFY_COVERAGE } from '../../modules/local/virome_classify_coverage/main'
include { OTU_MERGE                } from '../../modules/local/otu_merge/main'

workflow VIROME_CLASSIFICATION {
    take:
    ch_reads         // channel: [meta, [R1.fq, R2.fq]] - cleaned reads from host removal
    ch_mmseqs_db     // channel: [meta_db, path(db)] - MMseqs2 database
    ch_taxonomy_db   // path: taxonomy SQLite database
    ch_mask_bed      // path: masked regions BED file
    ch_segment_info  // path: segment info CSV (optional) - for coverage-based method

    main:
    ch_versions = Channel.empty()
    def method = params.virome_classification_method ?: 'coverage'  // Default to coverage (method B: best-hit + unmasked breadth + TPM); LCA is disabled

    //
    // Step 1: Separate R1 and R2 for individual MMseqs2 searches
    //
    ch_r1_reads = ch_reads.map { meta, reads ->
        def r1_meta = meta + [read: 'R1']
        [r1_meta, reads[0]]
    }

    ch_r2_reads = ch_reads.map { meta, reads ->
        def r2_meta = meta + [read: 'R2']
        [r2_meta, reads[1]]
    }

    //
    // Step 2: Run MMseqs2 on R1
    //
    MMSEQS_EASYSEARCH(
        ch_r1_reads.mix(ch_r2_reads),
        ch_mmseqs_db
    )
    ch_versions = ch_versions.mix(MMSEQS_EASYSEARCH.out.versions.first())

    //
    // Step 3: Group R1 and R2 results back together
    //
    ch_paired_results = MMSEQS_EASYSEARCH.out.tsv
        .map { meta, tsv ->
            def sample_id = meta.id
            def read_type = meta.read
            [sample_id, read_type, tsv]
        }
        .groupTuple(by: 0)
        .map { sample_id, read_types, tsvs ->
            // Reorder to ensure R1 is first, R2 is second
            def meta = [id: sample_id]
            def r1_idx = read_types.findIndexOf { it == 'R1' }
            def r2_idx = read_types.findIndexOf { it == 'R2' }
            [meta, tsvs[r1_idx], tsvs[r2_idx]]
        }

    //
    // Step 4: Run virome classification (LCA or Coverage-based)
    //
    // Coverage (method B) is the only supported classification path. The LCA
    // branch is disabled — see the commented-out block below to re-enable.
    if (method != 'coverage') {
        error "virome_classification_method='${method}' is not supported; only 'coverage' (method B) is enabled. (LCA was retired.)"
    }
    VIROME_CLASSIFY_COVERAGE(
        ch_paired_results,
        ch_taxonomy_db,
        ch_mask_bed,
        ch_segment_info
    )
    ch_versions = ch_versions.mix(VIROME_CLASSIFY_COVERAGE.out.versions.first())
    ch_lca_output = VIROME_CLASSIFY_COVERAGE.out.classification
    ch_kraken_output = VIROME_CLASSIFY_COVERAGE.out.kraken
    ch_kreport_output = VIROME_CLASSIFY_COVERAGE.out.kreport
    // coverage_abundance.tsv carries read_count + TPM (per-taxon); the plain
    // _abundance.tsv is only a fraction summary, so use the coverage one for the
    // TPM/read_count OTU merge.
    ch_abundance_output = VIROME_CLASSIFY_COVERAGE.out.coverage_abundance

    // ---- DISABLED: LCA classification branch (kept for reference) ----
    // } else {
    //     VIROME_CLASSIFY(
    //         ch_paired_results,
    //         ch_taxonomy_db,
    //         ch_mask_bed
    //     )
    //     ch_versions = ch_versions.mix(VIROME_CLASSIFY.out.versions.first())
    //     ch_lca_output = VIROME_CLASSIFY.out.lca
    //     ch_kraken_output = VIROME_CLASSIFY.out.kraken
    //     ch_kreport_output = VIROME_CLASSIFY.out.kreport
    //     ch_abundance_output = VIROME_CLASSIFY.out.abundance
    // }

    //
    // Step 5: Collect all classification results for OTU table generation
    //
    if (params.generate_otu_table) {
        ch_all_lca = ch_lca_output
            .map { meta, lca -> lca }
            .collect()

        // per-sample coverage_abundance.tsv (read_count + TPM) → merged into
        // sample × taxon abundance matrices alongside the read-count OTU tables
        ch_all_abundance = ch_abundance_output
            .map { meta, ab -> ab }
            .collect()
            .ifEmpty([])

        OTU_MERGE(
            ch_all_lca,
            ch_all_abundance,
            ch_taxonomy_db
        )
        ch_versions = ch_versions.mix(OTU_MERGE.out.versions)

        ch_otu_raw = OTU_MERGE.out.raw
        ch_otu_genus = OTU_MERGE.out.genus
        ch_otu_species = OTU_MERGE.out.species
    } else {
        ch_otu_raw = Channel.empty()
        ch_otu_genus = Channel.empty()
        ch_otu_species = Channel.empty()
    }

    emit:
    lca       = ch_lca_output       // channel: [meta, classification.csv]
    kraken    = ch_kraken_output    // channel: [meta, .kraken]
    kreport   = ch_kreport_output   // channel: [meta, .kreport]
    abundance = ch_abundance_output // channel: [meta, _abundance.tsv]
    otu_raw   = ch_otu_raw          // path: otu_table_raw.csv
    otu_genus = ch_otu_genus        // path: otu_table_genus.csv
    otu_species = ch_otu_species    // path: otu_table_species.csv
    versions  = ch_versions         // channel: [versions.yml]
}
