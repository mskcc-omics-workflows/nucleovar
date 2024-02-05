process MUTECT_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/msk-access/postprocessing_variant_calls:0.2.3"

    input:
    tuple val(meta), path(vardict_vcf_file), path(bams)
    //tuple val(meta), path(mutect_vcf_file), path(mutect_txt_file)
    //path(reference_fasta)
    

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
    pv mutect1 --help
    
    
    """


    // pv vardict single filter \
    // --inputVcf ${vardict_vcf_file} \
    // --tsampleName "${bam_filenames}" \
    // ${args}

    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     python: \$(python --version | sed 's/Python //g')
    // END_VERSIONS
    // stub:
    // task.ext.when == null || task.ext.when
    // def args = task.ext.args ?: ''
    // def args2 = task.ext.args2 ?: ''

    // """
    // touch 
    // """
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     python: \$(python --version | sed 's/Python //g')
    // END_VERSIONS
}