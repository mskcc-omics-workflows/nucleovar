process VARDICT_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/msk-access/postprocessing_variant_calls:0.2.3"

    input:
    tuple val(meta), path(vardict_vcf_file), path(bams)
    val(bam_filenames)
    

    output:
    path("*.vcf"),                     emit: filtered_vardict_vcf
    path("*_complex.vcf"),             emit: complex_variants_vardict_vcf
    path("*.txt"),                     emit: std_vardict_filter_output
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def paired_mode = bams instanceof List && bams.size() == 2 ? true : false
    

    """
    
    if [ "${paired_mode}" = true ]; then
        pv vardict case-control filter \
        --inputVcf ${vardict_vcf_file} \
        --tsampleName "${bam_filenames}" \
        ${args}
    else
        pv vardict single filter \
        --inputVcf ${vardict_vcf_file} \
        --tsampleName "${bam_filenames}" \
        ${args}
    fi
    
    """









    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     python: \$(python --version | sed 's/Python //g')
    // END_VERSIONS

    // stub:
    // task.ext.when == null || task.ext.when
    // def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    // def vardict_vcf_file = vardict_vcf_file ? "--inputVcf ${vardict_vcf_file}" : ''

    // """
    // mkdir vardict_filtered_output/${prefix} 
    // touch vardict_filtered_output/${prefix}/${prefix}.vcf
    // touch vardict_filtered_output/${prefix}/${prefix}.complex.vcf
    // touch vardict_filtered_output/${prefix}/${prefix}.txt
    // """
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     python: \$(python --version | sed 's/Python //g')
    // END_VERSIONS
}