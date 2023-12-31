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

process QUALIMAP {
    tag "$meta.id"
    label 'process_high_memory'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconda::qualimap=2.2.2d"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--hdfd78af_2':
        'biocontainers/qualimap:2.2.2d--hdfd78af_2' }"

    input:
    tuple val(meta), path(bam)
    //path gff

    output:
    tuple val(meta), path("*")          , emit: results
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args     ?: ''
    def prefix  = task.ext.prefix   ?: "${meta.id}"
    //def regions = gff ? "--gff $gff": ''
    def memory  = (task.memory.mega*0.8).intValue() + 'M'
    //mem_size  = task.mem_size      ?: ''

    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    //def collect_pairs = meta.single_end ? '' : '--collect-overlap-pairs'
    //def memory = (task.memory.mega*0.8).intValue() + 'M'
    //def regions = gff ? "--gff $gff" : ''

    //def strandedness = 'non-strand-specific'
    //if (meta.strandedness == 'forward') {
    //    strandedness = 'strand-specific-forward'
    //} else if (meta.strandedness == 'reverse') {
    //    strandedness = 'strand-specific-reverse'
    //}

    // unset DISPLAY
    // mkdir -p tmp
    // export_JAVA_OPTIONS=Djava.io.tmpdir=./tmp

    """
    qualimap \\
        bamqc \\
        -bam $bam \\
        -outdir $prefix \\
        -nt $task.cpus \\
        --java-mem-size=$memory

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap --version 2>&1) | sed 's/^.*QualiMap //; s/Using.*\$//')
    END_VERSIONS
    """
}
