process VIROME_CLASSIFY_COVERAGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11'}"

    input:
    tuple val(meta), path(r1_result), path(r2_result)
    path taxonomy_db
    path mask_bed
    path segment_info, stageAs: 'segment_info.csv'

    output:
    tuple val(meta), path("${prefix}_lca_classification.csv"), emit: classification
    tuple val(meta), path("${prefix}_coverage_abundance.tsv"), emit: coverage_abundance, optional: true
    tuple val(meta), path("${prefix}_coverage_summary.tsv"), emit: coverage_summary, optional: true
    tuple val(meta), path("${prefix}.kraken"), emit: kraken, optional: true
    tuple val(meta), path("${prefix}.kreport"), emit: kreport, optional: true
    tuple val(meta), path("${prefix}_abundance.tsv"), emit: abundance, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Uses the validated coverage path: virome_classifier.cli.classify --mode coverage
    // (CoverageBasedClassifier3). Segment info is read from the taxonomy DB, so the
    // staged segment_info file is currently informational; threshold knobs are tuned
    // via task.ext.args until exposed as first-class CLI flags (Phase 2).
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def min_identity = task.ext.min_identity ?: params.min_identity ?: 0.85
    def min_length = task.ext.min_length ?: params.min_length ?: 60
    def min_query_coverage = task.ext.min_query_coverage ?: params.min_query_coverage ?: 0.5
    def min_unmasked_coverage = task.ext.min_coverage ?: params.min_coverage ?: 0.05
    def classification_rank = task.ext.classification_rank ?: params.classification_rank ?: 'species'
    def multi_mapping_mode = task.ext.multi_mapping_mode ?: params.multi_mapping_mode ?: 'all'
    def depth_entropy = (task.ext.use_depth_entropy ?: params.use_depth_entropy) ? '--use-depth-entropy' : ''
    """
    export PYTHONPATH="${moduleDir}/../../../scripts:\${PYTHONPATH:-}"
    python -m virome_classifier.cli.classify \\
        --mode coverage \\
        --r1 ${r1_result} \\
        --r2 ${r2_result} \\
        --taxonomy ${taxonomy_db} \\
        --mask ${mask_bed} \\
        --output . \\
        --sample ${prefix} \\
        --min-identity ${min_identity} \\
        --min-length ${min_length} \\
        --min-query-coverage ${min_query_coverage} \\
        --min-unmasked-coverage ${min_unmasked_coverage} \\
        --classification-rank ${classification_rank} \\
        --multi-mapping-mode ${multi_mapping_mode} \\
        ${depth_entropy} \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        virome_classifier: \$(python -c "import virome_classifier; print(virome_classifier.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_lca_classification.csv
    touch ${prefix}_coverage_abundance.tsv
    touch ${prefix}_coverage_summary.tsv
    touch ${prefix}.kraken
    touch ${prefix}.kreport
    touch ${prefix}_abundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
