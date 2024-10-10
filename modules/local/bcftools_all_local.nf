process BCFTOOLS_ALL_LOCAL {
    container 'ghcr.io/msk-access/bcftools:1.15.1'

    input:
    tuple val(meta), path(vcf), path(index)
    tuple val(meta), path(ref_fasta), path(ref_fasta_index)

    output:
    path("test.vcf.gz"), emit: vcf
    path('completed.txt')

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    bcftools norm --fasta-ref ${ref_fasta} --output test.vcf.gz --threads 1 ${vcf}
    touch completed.txt
    """
    //-f ref_fasta ${ref_fasta} ${vcf} --output ${prefix}.vcf.gz
    //-m ${vcf} --output ${prefix}.vcf.gz

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
