process MODKIT_PILEUP {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/ontresearch/wf-human-variation-methyl:shaa6e616571797d97ae2736c7ebdcb4613fe77f263"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("*.bed"), emit: bed
    tuple val(meta), path("*.log"), emit: log
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args  ?: ''
    def args2   = task.ext.args2 ?: ''
    def args3   = task.ext.args3 ?: ''
    def ref     = params.motif_CG || params.motif_CHG || params.motif_CHH || params.motif_A ? "--ref ${reference}" : ''

    """
    modkit pileup \\
    $bam \\
    $args \\
    $args2 \\
    $args3 \\
    $ref \\
    ${meta.id}.bed \\
    --log-filepath ${meta.id}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        modkit: \$(modkit --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
