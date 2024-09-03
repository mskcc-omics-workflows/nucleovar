process BCFTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/bcftools:1.15.1':
        'ghcr.io/msk-access/bcftools:1.15.1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: tbi
    path("versions.yml")           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bgzip ${vcf}
    bcftools \\
        index \\
        $args \\
        --threads $task.cpus \\
        ${vcf}.gz



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    touch ${vcf}.gz
    touch ${vcf}.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
