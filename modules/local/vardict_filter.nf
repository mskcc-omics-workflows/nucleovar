process VARDICT_FILTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.3':
        'ghcr.io/msk-access/postprocessing_variant_calls:0.2.3' }"


    input:
    val(meta)
    path(vardict_vcf_file)
    path(bams)


    output:
    path("*_STDfilter.vcf"),                     emit: std_vcf
    path("*_STDfilter_complex.vcf"),             emit: complex_variants_vardict_vcf
    path("*.txt"),                     emit: std_vardict_filter_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def paired_mode = bams instanceof List && bams.size() == 2 ? true : false


    """

    pv vardict case-control filter \
    --inputVcf ${vardict_vcf_file} \
    --tsampleName ${meta.case_id} \
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def paired_mode = bams instanceof List && bams.size() == 2 ? true : false

    """
    touch ${meta.id}_STDfilter.vcf
    touch ${meta.id}_STDfilter_complex.vcf
    touch ${meta.id}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
}
