process GUNZIP_FILES {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'staphb/htslib':
        'staphb/htslib' }"

    input:
    tuple val(meta), path(vcf)
    path(index)

    output:
    tuple val(meta), path("*.vcf"), emit: decomp_vcf
    tuple val(meta), path("*.tbi"), emit: decomp_tbi

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    bgzip -d ${vcf} > ${meta.id}.vcf
    bgzip -d ${index} > ${meta.id}.tbi
    """
}
