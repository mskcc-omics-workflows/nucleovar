process TAG_BY_VARIANT_ANNOTATION {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.4' :
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.4' }"

    input:
    tuple val(meta), path(access_filtered_maf)
    path(canonical_tx_ref)

    output:
    path("annotated_exonic_variants.maf"), emit: annotated_exonic
    path("annotated_silent_variants.txt"), emit: annotated_silent
    path("annotated_nonpanel_exonic_variants.txt"), emit: annotated_nonpanel_exonic
    path("annotated_nonpanel_silent_variants.txt"), emit: annotated_nonpanel_silent
    path("annotated_dropped_variants.txt"), emit: annotated_dropped

    path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    pv maf tag by_variant_classification --maf ${access_filtered_maf} --canonical_tx_ref ${canonical_tx_ref} --output_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """

    stub:
    """
    touch annotated_exonic_variants.maf
    touch annotated_silent_variants.txt
    touch annotated_nonpanel_exonic_variants.txt
    touch annotated_nonpanel_silent_variants.txt
    touch annotated_dropped_variants.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """
}
