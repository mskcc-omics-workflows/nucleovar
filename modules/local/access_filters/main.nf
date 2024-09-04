process ACCESS_FILTERS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.4' :
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.4' }"

    input:
    tuple val(meta), path(traceback_maf), path(anno_maf)

    output:
    tuple val(meta), path("*_filtered.maf"), emit: filtered_maf
    path("*_filtered_condensed.maf"), emit: condensed_filtered_maf
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools sort $tumor_bam -o ${meta.case_id}_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_filtered.maf
    touch ${meta.id}_filtered_condensed.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """
}
