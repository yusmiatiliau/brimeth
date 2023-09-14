// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process MODKIT {
    tag "$meta.id"
    label 'process_low'

    container "docker.io/ontresearch/wf-somatic-methyl:shac81dbea5f824cc43fd7aeb9ad99b4efe1503216e"

    input:
        path(fasta)

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    path("*.CG.bed")    , emit: CG.bed
    path("*.CHG.bed")   , emit: CHG.bed
    path("*.CHH.bed")   , emit: CHH.bed
    path("*.6mA.bed")   , emit: 6mA.bed
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    modkit motif.bed \\
    $fasta \\
    ${param.motif}
    ${meta.id}.bed \\
    $args \\
    $args2 \\
    --log-filepath ${meta.id}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(modkit --version 2>&1) | sed 's/^.*modkit //; s/Using.*\$//' ))
    END_VERSIONS
    """

    //stub:
    //def args = task.ext.args ?: ''

    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    //"""
    //touch ${prefix}.bam

    //cat <<-END_VERSIONS > versions.yml
    //"${task.process}":
    //    : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    //END_VERSIONS
    //"""
}
