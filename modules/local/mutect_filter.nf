process MUTECT_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/msk-access/postprocessing_variant_calls:postprocessing_test"

    input:
    tuple val(meta), path(mutect_vcf_file), path(mutect_txt_file)
    path(reference_fasta)
    

    output:
    path("*.mutect.vcf"),                     emit: mutect_filtered_vcf
    path("*.mutect.txt"),                     emit: std_mutect_filter_output
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    

    """
    
    pv mutect1 case-control filter --inputVcf ${mutect_vcf_file} --inputTxt ${mutect_txt_file} --refFasta ${reference_fasta} --tsampleName ${meta.case_id}
    
    """
}