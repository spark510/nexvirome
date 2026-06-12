/*
========================================================================================
    IMPORT MODULES
========================================================================================
    - 이 서브워크플로우에서 사용할 모든 모듈을 불러옵니다.
----------------------------------------------------------------------------------------
*/

include { FASTQC          } from '../../modules/nf-core/fastqc/main'

include { TRIMMOMATIC     } from '../../modules/nf-core/trimmomatic/main'
include { FASTP           } from '../../modules/nf-core/fastp/main'
include { CUTADAPT        } from '../../modules/nf-core/cutadapt/main'

include { BOWTIE2_BUILD   } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN   } from '../../modules/nf-core/bowtie2/align/main'

workflow SHORTREAD_PREPROCESSING {

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        INPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    take:
    ch_input_reads     // Channel: Tuples of [meta, [read1, read2]] or [meta, [read1]] from samplesheet
    run_cutadapt       // Boolean: Whether to perform adapter removal using Cutadapt
    run_host_removal   // Boolean: Whether to perform host read removal using Bowtie2
    ch_host_index      // Channel: Bowtie2 index files (if host removal is enabled)
    ch_host_fasta      // Channel: Reference FASTA file used for Bowtie2 alignment (header)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        MAIN WORKFLOW
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    main:
    ch_trim_log     = Channel.empty()
    ch_cutadapt_log = Channel.empty()
    ch_bowtie_log   = Channel.empty()
    ch_versions     = Channel.empty()

    //
    // Step 1: Trimmomatic - Perform quality trimming
    //
    // TRIMMOMATIC(
    //     ch_input_reads
    // )
    // ch_trimmed_reads = TRIMMOMATIC.out.trimmed_reads
    // ch_trim_log      = TRIMMOMATIC.out.trim_log.collect { it[1] }
    // ch_versions      = ch_versions.mix(TRIMMOMATIC.out.versions.first())

    FASTP(
        ch_input_reads,
        [], // No adapter FASTA provided, will use default settings
        false,
        true,
        false
    )
    ch_trimmed_reads = FASTP.out.reads
    ch_trim_log      = FASTP.out.log.collect { it[1] }
    ch_versions      = ch_versions.mix(FASTP.out.versions.first())

    //
    // Step 2: Optional - Run Cutadapt for adapter removal
    //
    if (run_cutadapt) {
        CUTADAPT(
            ch_trimmed_reads
        )
        ch_cutadapt_reads = CUTADAPT.out.reads
        ch_cutadapt_log   = CUTADAPT.out.log.collect { it[1] }
        ch_versions       = ch_versions.mix(CUTADAPT.out.versions.first())
    } else {
        ch_cutadapt_reads = ch_trimmed_reads
    }

    //
    // Step 3: Optional - Run Bowtie2 for host (or PhiX) read removal
    //
    if (run_host_removal) {
        BOWTIE2_ALIGN(
            ch_cutadapt_reads,
            ch_host_index,
            ch_host_fasta,
            true, // save_unaligned
            true   // sort_bam
        )
        ch_processed_reads = BOWTIE2_ALIGN.out.fastq
        ch_bowtie_log  = BOWTIE2_ALIGN.out.log.collect { it[1] }
        ch_versions    = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
        ch_processed_reads.view()
    } else {
        ch_processed_reads = ch_cutadapt_reads
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        OUTPUTS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    emit:
    reads    = ch_processed_reads                           // Cleaned reads after trimming, optional adapter removal, and optional host filtering
    logs     = ch_trim_log.mix(ch_cutadapt_log).mix(ch_bowtie_log) // Logs from Trimmomatic, Cutadapt, and Bowtie2 (if applicable)
    versions = ch_versions                                // Software version tracking
}
