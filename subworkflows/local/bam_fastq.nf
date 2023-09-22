//
// Turns bam files to fastq and runs qc
//

include { FASTQC                      } from '../../modules/nf-core/fastqc/main'
include { SAMTOOLS_BAM2FQ             } from '../../modules/nf-core/samtools/bam2fq/main'
include { CHOPPER                     } from '../../modules/nf-core/chopper/main'

workflow BAM_TO_FASTQ {

    take:
    my_samples //


    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Run samtools/bam2fq
    //
    SAMTOOLS_BAM2FQ (
        my_samples,
        params.split
    )
    ch_reads    = SAMTOOLS_BAM2FQ.out.reads
    ch_versions = ch_versions.mix(SAMTOOLS_BAM2FQ.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_reads
    )
    ch_fastqc_zip = FASTQC.out.zip.collect{it[1]}.ifEmpty([])
    ch_versions   = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run chopper
    //
    CHOPPER (
        ch_reads
    )
    ch_chopper_fastq = CHOPPER.out.fastq // channel: [ val(meta), ?? ]
    ch_versions      = ch_versions.mix(CHOPPER.out.versions.first())

    emit:
    reads         = ch_reads
    chopper_fastq = ch_chopper_fastq
    versions      = ch_versions
    fastqc_zip    = ch_fastqc_zip
}
