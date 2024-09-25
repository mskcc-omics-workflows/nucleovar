process MUTECT_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/msk-access/postprocessing_variant_calls:mutect_filter_0.0.1"

    input:
    tuple val(meta), path(mutect_vcf_file), path(mutect_txt_file)
    path(reference_fasta)


    output:
    path("*_filtered.vcf"),                     emit: mutect_filtered_vcf
    path("*_filtered.txt"),                     emit: std_mutect_filter_output
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''


    """
    chmod -R 777 .
    pv mutect1 case-control filter --inputVcf ${mutect_vcf_file} --inputTxt ${mutect_txt_file} --refFasta ${reference_fasta} --tsampleName ${meta.case_id} --outDir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/python v//g")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_filtered.vcf
    touch ${meta.id}_filtered.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
