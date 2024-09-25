process BCFTOOLS_CONCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/bcftools:1.15.1':
        'ghcr.io/msk-access/bcftools:1.15.1' }"

    input:
    tuple val(meta), path(vcf1), path(vcf2)

    output:
    tuple val(meta), path("*.gz"), emit: concat_vcf
    path  "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    bcftools index -t ${vcf1}
    bcftools index -t ${vcf2}

    bcftools concat \\
        --output ${prefix}_vardict_concat.vcf.gz \\
        $args \\
        --threads $task.cpus \\
        ${vcf1} ${vcf2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_vardict_concat.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
