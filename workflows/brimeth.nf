/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowBrimeth.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULES: Consisting of local modules
//
include   { QUALIMAP              } from '../modules/local/qualimap.nf'
include   { MODKIT_PILEUP         } from '../modules/local/modkit_pileup.nf'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
//include { INPUT_CHECK } from '../subworkflows/local/input_check'
include   { BAM_TO_FASTQ          } from '../subworkflows/local/bam_fastq.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
//include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
//include { SAMTOOLS_BAM2FQ             } from '../modules/nf-core/samtools/bam2fq/main'
//include { CHOPPER                     } from '../modules/nf-core/chopper/main'
include { MINIMAP2_ALIGN              } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'

include { BAM_STATS_SAMTOOLS          } from '../subworkflows/nf-core/bam_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow BRIMETH {

    ch_versions = Channel.empty()

    //
    // Folder input
    //
    // Single pod5 file is easist

    //
    // Download model
    //
    // https://github.com/christopher-hakkaart/nanoseq/blob/local_to_nfcore/modules/local/dorado_model.nf

    //
    // Run Dorado
    //
    // https://github.com/christopher-hakkaart/nanoseq/blob/local_to_nfcore/modules/local/dorado_basecaller.nf

    //
    // Create input channels
    //
    Channel
    .fromPath( params.input )
    .splitCsv( header: true )
    .map{ row -> [[id:row.sample], file(row.bam) ] }
    .set{ my_samples }

    //
    // SUBWORKFLOW: Bam to Fastq + QC
    //
    BAM_TO_FASTQ (
        my_samples
    )
    ch_fastqc_zip    = BAM_TO_FASTQ.out.fastqc_zip
    ch_chopper_fastq = BAM_TO_FASTQ.out.chopper_fastq
    ch_versions      = ch_versions.mix(BAM_TO_FASTQ.out.versions)

    //
    // MODULE: Run minimap2/align
    //
    MINIMAP2_ALIGN (
        ch_chopper_fastq,
        params.fasta,
        params.bam_format,
        params.cigar_paf_format,
        params.cigar_bam
    )
    ch_bam = MINIMAP2_ALIGN.out.bam // channel: [ val(meta), file(bai) ]
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    //
    // MODULE: Run Samtools index
    //
    SAMTOOLS_INDEX {
        ch_bam
    }
    ch_bai = SAMTOOLS_INDEX.out.bai // channel: [ val(meta), file(bai) ]
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // Mix bam and bai
    ch_bam_bai = ch_bam.join(ch_bai) // channel: [ val(meta), path(bam), path(bai) ]

    //
    // MODULE: Run Qualimap
    //
    QUALIMAP (
        ch_bam,
        //[[:],params.gff]
    )

    //
    // SUBWORKFLOW: Run Samtools_stats
    //
    BAM_STATS_SAMTOOLS (
        ch_bam_bai,
        [[:], params.fasta]
    )

    //
    // MODULE: Run Modkit
    //
    MODKIT_PILEUP (
        ch_bam_bai,
        params.fasta
    )
    ch_bed = MODKIT_PILEUP.out.bed // channel: [ val(meta), file(bed) ]
    ch_versions = ch_versions.mix(MODKIT_PILEUP.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowBrimeth.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowBrimeth.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(ch_fastqc_zip)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
