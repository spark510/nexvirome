process VIROME_CLASSIFY {
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

    output:
    tuple val(meta), path("${prefix}_lca_classification.csv"), emit: lca
    tuple val(meta), path("${prefix}.kraken"), emit: kraken, optional: true  // LCA-based Kraken format
    tuple val(meta), path("${prefix}.kreport"), emit: kreport, optional: true  // LCA-based hierarchical report
    tuple val(meta), path("${prefix}_abundance.tsv"), emit: abundance, optional: true  // LCA-based abundance
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // LOCKED Method-B params (GOLDEN_RULE). Each maps 1:1 onto a classify CLI flag.
    def mode             = task.ext.mode             ?: params.virome_classification_method ?: 'lca'
    def min_identity     = task.ext.min_identity     ?: params.min_identity          ?: 0.85
    def min_length       = task.ext.min_length       ?: params.min_length            ?: 60
    def min_query_cov    = task.ext.min_query_cov    ?: params.min_query_coverage    ?: 0.5
    def max_evalue       = task.ext.max_evalue       ?: params.max_evalue            ?: 1e-3
    def min_unmasked_cov = task.ext.min_unmasked_cov ?: params.min_unmasked_coverage ?: 0.01
    def min_read_count   = task.ext.min_read_count   ?: params.min_read_count        ?: 3
    def min_rel_abund    = task.ext.min_rel_abund    ?: params.min_rel_abundance     ?: 0.0
    def mm_mode          = task.ext.mm_mode          ?: params.multi_mapping_mode    ?: 'best_hit'
    def class_rank       = task.ext.class_rank       ?: params.classification_rank   ?: 'species'
    """
    export PYTHONPATH="${moduleDir}/../../../scripts:\${PYTHONPATH:-}"
    python -m virome_classifier.cli.classify \\
        --r1 ${r1_result} \\
        --r2 ${r2_result} \\
        --taxonomy ${taxonomy_db} \\
        --mask ${mask_bed} \\
        --output . \\
        --sample ${prefix} \\
        --mode ${mode} \\
        --min-identity ${min_identity} \\
        --min-length ${min_length} \\
        --min-query-coverage ${min_query_cov} \\
        --max-evalue ${max_evalue} \\
        --min-unmasked-coverage ${min_unmasked_cov} \\
        --min-read-count ${min_read_count} \\
        --min-rel-abundance ${min_rel_abund} \\
        --multi-mapping-mode ${mm_mode} \\
        --classification-rank ${class_rank} \\
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
    touch ${prefix}.kraken
    touch ${prefix}.kreport
    touch ${prefix}_abundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
