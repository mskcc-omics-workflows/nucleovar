process MAF_PROCESSING {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.7' :
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.7' }"

    input:
    tuple val(meta), path(genotyped_maf)
    path(rules_file)
    path(hotspots)

    output:
    tuple val(meta), path("*.maf"), emit: maf
    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pv maf tag maf_processing --maf ${genotyped_maf} --rules_json ${rules_file} --hotspots ${hotspots}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """
}
