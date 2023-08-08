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
    label 'process_medium'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "bioconda::qualimap=2.2.2d"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--hdfd78af_2':
        'biocontainers/qualimap:2.2.2d--hdfd78af_2' }"

    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.html"), optional: true, emit: html
    tuple val(meta), path("*.txt"), optional: true, emit: txt
    path "versions.yml"           , emit: versions

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    output:
    tuple val(meta), path("${prefix}"), emit: results
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    //def collect_pairs = meta.single_end ? '' : '--collect-overlap-pairs'
    //def memory = (task.memory.mega*0.8).intValue() + 'M'
    //def regions = gff ? "--gff $gff" : ''

    //def strandedness = 'non-strand-specific'
    //if (meta.strandedness == 'forward') {
    //    strandedness = 'strand-specific-forward'
    //} else if (meta.strandedness == 'reverse') {
    //    strandedness = 'strand-specific-reverse'
    //}


    """
    qualimap \\
        bamqc \\
        -bam $ch_bam \\
        -@ $task.cpus \\
        -outdir $prefix \\
        -args \\
        

        qualimap bamqc -bam /nesi/nobackup/brins03581/Cen/SB09/dorado/SB09_5mC_merged_filteredQ15_aligned_sorted.bam
        -outdir /nesi/nobackup/brins03581/Cen/SB09/dorado/qualimap_SB09
        --java-mem-size=48G

    unset DISPLAY
    mkdir -p tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=$memory \\
        bamqc \\
        $args \\
        -bam $bam \\
        $regions \\
        -p $strandedness \\
        $collect_pairs \\
        -outdir $prefix \\
        -nt $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.suffix ? "${meta.id}${task.ext.suffix}" : "${meta.id}"
    """
    mkdir -p $prefix/css
    mkdir $prefix/images_qualimapReport
    mkdir $prefix/raw_data_qualimapReport
    cd $prefix/css
    touch agogo.css
    touch basic.css
    touch bgtop.png
    touch comment-close.png
    touch doctools.js
    touch down-pressed.png
    touch jquery.js
    touch plus.png
    touch qualimap_logo_small.png
    touch searchtools.js
    touch up.png
    touch websupport.js
    touch ajax-loader.gif
    touch bgfooter.png
    touch comment-bright.png
    touch comment.png
    touch down.png
    touch file.png
    touch minus.png
    touch pygments.css
    touch report.css
    touch underscore.js
    touch up-pressed.png
    cd ../images_qualimapReport/
    touch genome_coverage_0to50_histogram.png
    touch genome_coverage_quotes.png
    touch genome_insert_size_across_reference.png
    touch genome_mapping_quality_histogram.png
    touch genome_uniq_read_starts_histogram.png
    touch genome_coverage_across_reference.png
    touch genome_gc_content_per_window.png
    touch genome_insert_size_histogram.png
    touch genome_reads_clipping_profile.png
    touch genome_coverage_histogram.png
    touch genome_homopolymer_indels.png
    touch genome_mapping_quality_across_reference.png
    touch genome_reads_content_per_read_position.png
    cd ../raw_data_qualimapReport
    touch coverage_across_reference.txt
    touch genome_fraction_coverage.txt
    touch insert_size_histogram.txt
    touch mapped_reads_nucleotide_content.txt
    touch coverage_histogram.txt
    touch homopolymer_indels.txt
    touch mapped_reads_clipping_profile.txt
    touch mapping_quality_across_reference.txt
    touch duplication_rate_histogram.txt
    touch insert_size_across_reference.txt
    touch mapped_reads_gc-content_distribution.txt
    touch mapping_quality_histogram.txt
    cd ../
    touch genome_results.txt
    touch qualimapReport.html
    cd ../

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//')
    END_VERSIONS
    """

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        : \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
