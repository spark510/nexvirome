process OTU_MERGE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'biocontainers/python:3.11'}"

    input:
    path lca_files
    path taxonomy_db

    output:
    path "otu_tables/otu_table_raw.csv", emit: raw
    path "otu_tables/otu_table_genus.csv", emit: genus
    path "otu_tables/otu_table_species.csv", emit: species, optional: true
    path "otu_tables/otu_table_family.csv", emit: family, optional: true
    path "otu_tables/otu_table_phage_host.csv", emit: phage_host, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ranks = task.ext.ranks ?: "genus species"
    def min_count = task.ext.min_count ?: 10
    def min_samples = task.ext.min_samples ?: 1
    def normalize = task.ext.normalize ? '--normalize' : ''
    // Also emit phage→host-genus rolled-up table by default (additive; the
    // per-rank tables are unaffected). Disable with ext.phage_host = false.
    def phage_host = (task.ext.phage_host == null || task.ext.phage_host) ? '--phage-host' : ''
    """
    # Create input directory for LCA files
    mkdir -p lca_input

    # Copy all LCA files to input directory
    for file in ${lca_files}; do
        cp "\$file" lca_input/
    done

    # Run OTU merge
    export PYTHONPATH="${moduleDir}/../../../scripts:\${PYTHONPATH:-}"
    python -m virome_classifier.cli.merge_otu \\
        --input-dir lca_input \\
        --pattern "*_lca_classification.csv" \\
        --taxonomy ${taxonomy_db} \\
        --output otu_tables \\
        --ranks ${ranks} \\
        --min-count ${min_count} \\
        --min-samples ${min_samples} \\
        ${normalize} \\
        ${phage_host} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        virome_classifier: \$(python -c "import virome_classifier; print(virome_classifier.__version__)" 2>/dev/null || echo "dev")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p otu_tables
    touch otu_tables/otu_table_raw.csv
    touch otu_tables/otu_table_genus.csv
    touch otu_tables/otu_table_species.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
    END_VERSIONS
    """
}
