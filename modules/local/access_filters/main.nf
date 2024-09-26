process ACCESS_FILTERS {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.3.1' :
        'ghcr.io/msk-access/postprocessing_variant_calls:0.3.1' }"

    input:
    tuple val(meta), path(traceback_maf), path(anno_maf), path(blocklist)

    output:
    tuple val(meta), path("*_filtered.maf"), emit: filtered_maf
    path("*_filtered_condensed.maf"), emit: condensed_filtered_maf
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pv maf filter access_filters -f ${traceback_maf} -a ${anno_maf} -ts ${meta.case_id}  -ns ${meta.control_id} -bl ${blocklist} --output ${meta.id}

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
