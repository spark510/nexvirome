/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { UNTAR                       } from '../modules/nf-core/untar/main'

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build/main'

// include { MMSEQS_CREATEDB         } from '../modules/nf-core/mmseqs/createdb/main'
// include { MMSEQS_EASYSEARCH       } from '../modules/nf-core/mmseqs/easysearch/main'
// include { SEQKIT_FQ2FA           } from '../modules/nf-core/seqkit/fq2fa/main'
include { SHORTREAD_PREPROCESSING } from '../subworkflows/local/shortread_preprocessing'

include { paramsSummaryMap        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText  } from '../subworkflows/local/utils_nexvirome_pipeline'
include { PROFILING               } from '../subworkflows/local/profiling'
include { VIROME_CLASSIFICATION   } from '../subworkflows/local/virome_classification'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
    nextflow run <repo>/main.nf \
     --input sample_sheet2.csv \
     --databases database.csv \
    --host_index /home/share/bowtie2_db/GCA_009914755.4-bowtie2/ \
     --outdir results_test -profile conda,test -resume
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow NEXVIROME {


    take:
    ch_samplesheet // Channel of [meta, [R1, R2]] or [meta, [R1]]
    // ch_databases   // channel: databases from --databases

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // //  // Validate and decompress databases
    // ch_dbs_for_untar = ch_databases
    //     .branch { db_meta, db_path ->
    //         if ( !db_meta.db_type ) {
    //             db_meta = db_meta + [ db_type: "short;long" ]
    //         }
    //         untar: db_path.name.endsWith( ".tar.gz" )
    //         skip: true
    //     }
    // // Filter the channel to untar only those databases for tools that are selected to be run by the user.
    // // Also, to ensure only untar once per file, group together all databases of one file
    // ch_inputdb_untar = ch_dbs_for_untar.untar
    //     .filter { db_meta, db_path ->
    //         params[ "run_${db_meta.tool}" ]
    //     }
    //     .groupTuple(by: 1)
    //     .map {
    //         meta, dbfile ->
    //             def new_meta = [ 'id': dbfile.baseName ] + [ 'meta': meta ]
    //         [new_meta , dbfile ]
    //     }

    // // Untar the databases
    // UNTAR (ch_inputdb_untar )
    // ch_versions = ch_versions.mix( UNTAR.out.versions.first() )

    // // Spread out the untarred and shared databases
    // ch_outputdb_from_untar = UNTAR.out.untar
    //     .map {
    //         meta, db ->
    //         [meta.meta, db]
    //     }
    //     .transpose(by: 0)

    // ch_final_dbs = ch_dbs_for_untar.skip
    //                 .mix( ch_outputdb_from_untar  )
    //                 .map { db_meta, db ->
    //                     def corrected_db_params = db_meta.db_params ? [ db_params: db_meta.db_params ] : [ db_params: '' ]
    //                     [ db_meta + corrected_db_params, db ]
    //                 }

    // ---- 1. Host removal 관련 입력 체크 & channel 준비 ----
    if (!params.skip_host_removal) {
        if (!params.host_fasta && !params.host_index) {
            error """
            [ERROR] Host removal is enabled, but neither `host_fasta` nor `host_index` was provided.
            Please provide either --host_fasta or --host_index, or use --skip_host_removal.
            """
        }
    }

    // 빈 채널 기본값
    // host removal 관련 채널 준비
    ch_bowtie2_index = Channel.empty()
    ch_host_fasta = Channel.empty()

    if (!params.skip_host_removal) {
        if (params.host_index) {
            // 인덱스만 제공 (예: params.host_index = "index_dir")
            def host_index_dir = file(params.host_index, type: 'dir')
            def meta_index = [id: host_index_dir.name]
            ch_bowtie2_index = Channel.of(tuple(meta_index, host_index_dir)).first()
            ch_host_fasta = Channel.of(tuple(meta_index, file("/dev/null"))).first()
        } else if (params.host_fasta) {
            def meta_fasta = [id: file(params.host_fasta).simpleName]
            ch_host_fasta = Channel.of(tuple(meta_fasta, file(params.host_fasta))).first()
            BOWTIE2_BUILD(ch_host_fasta)
            ch_bowtie2_index = BOWTIE2_BUILD.out.index
        } else {
            error "Host removal enabled, but neither --host_index nor --host_fasta provided."
        }
    } else {
        ch_bowtie2_index = Channel.empty()
        ch_host_fasta = Channel.empty()
    }

    // ---- 2. Short read preprocessing ----
    SHORTREAD_PREPROCESSING(
        ch_samplesheet,
        !params.skip_cutadapt,
        !params.skip_host_removal,
        ch_bowtie2_index,
        ch_host_fasta
    )
    ch_cleaned_reads = SHORTREAD_PREPROCESSING.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(SHORTREAD_PREPROCESSING.out.logs)
    ch_versions      = ch_versions.mix(SHORTREAD_PREPROCESSING.out.versions)

    // ---- 3. Virome Classification ----
    if (params.run_virome_classify && params.mmseqs_database && params.taxonomy_db && params.mask_bed) {
        // Prepare MMseqs2 database channel
        // MMseqs2 DB is a prefix, gather all related files
        def mmseqs_db_base = file(params.mmseqs_database)
        def mmseqs_db_dir = mmseqs_db_base.parent
        def mmseqs_db_prefix = mmseqs_db_base.name
        def meta_mmseqs = [id: mmseqs_db_prefix]
        ch_mmseqs_db = Channel.of(tuple(meta_mmseqs, mmseqs_db_dir)).first()

        // Prepare taxonomy and mask channels
        ch_taxonomy_db = Channel.fromPath(params.taxonomy_db, checkIfExists: true).first()
        ch_mask_bed = Channel.fromPath(params.mask_bed, checkIfExists: true).first()

        // Optional segment-info CSV (coverage method). When not provided, pass a
        // NO_FILE placeholder so the coverage module's input arity is satisfied;
        // CoverageBasedClassifier3 otherwise reads segment info from the taxonomy DB.
        ch_segment_info = params.segment_info ?
            Channel.fromPath(params.segment_info, checkIfExists: true).first() :
            Channel.fromPath("${projectDir}/assets/NO_FILE", checkIfExists: true).first()

        // Run virome classification
        VIROME_CLASSIFICATION(
            ch_cleaned_reads,
            ch_mmseqs_db,
            ch_taxonomy_db,
            ch_mask_bed,
            ch_segment_info
        )
        ch_versions = ch_versions.mix(VIROME_CLASSIFICATION.out.versions)
    }

    // ---- 4. (Optional) downstream analysis ----
    // ch_cleaned_reads 를 downstream 모듈로 넘기세요.

    // PROFILING ( ch_cleaned_reads, ch_final_dbs )
    // ch_versions = ch_versions.mix( PROFILING.out.versions )



    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'virome_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}