process MUTECT2_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/msk-access/postprocessing_variant_calls:access_filters_fix"

    input:
    tuple val(meta), path(mutect_vcf_file)


    output:
    path("*_filtered.vcf"),                     emit: mutect_filtered_vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''


    """
    pv mutect2 case-control filter --inputVcf ${mutect_vcf_file} --tsampleName ${meta.case_id} --outDir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.mutect2_filtered.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
