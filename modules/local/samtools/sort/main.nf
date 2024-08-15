process SAMTOOLS_SORT {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/samtools:1.6' :
        'ghcr.io/msk-access/samtools:1.6' }"

    input:
    tuple val(meta), path(tumor_bam)

    output:
    tuple val(meta), path("*_sorted.bam"), emit: sorted_bam
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when
    
    script:
    """
    samtools sort $tumor_bam -o ${meta.case_id}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -e "s/samtools v//g")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.case_id}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -e "s/samtools v//g")
    END_VERSIONS
    """
}
