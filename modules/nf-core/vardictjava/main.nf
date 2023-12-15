process VARDICTJAVA {
    //tag "$meta.id"
    label 'vardictjava'

    // NOTE: may not need this because we have a docker image 
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vardict-java:1.8.3--hdfd78af_0':
        'biocontainers/vardict-java:1.8.3--hdfd78af_0' }"

    input:
    tuple val(tumor_sample_name), val(normal_sample_name), path(bams), path(bed), path(fasta)

    output:
    tuple val(sample_name), path("*.vcf"), emit: vardict_vcf
    //path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: '-c 1 -S 2 -E 3'
    def args2 = task.ext.args2 ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"

    def somatic = bams instanceof ArrayList && bams.size() <= 2 ? true : false
    def input = somatic ? "-b \"${bams[0]}|${bams[1]}\"" : "-b ${bams}"
    def filter = somatic ? "testsomatic.R" : "teststrandbias.R"
    def convert_to_vcf = somatic ? "var2vcf_paired.pl" : "var2vcf_valid.pl"
    def min_num_variant_reads = "-r ${params.min_num_variant_reads}"
    def column_for_region_start = "-S ${params.start_reg_col}"
    def column_for_region_end = "-E ${params.end_reg_col}"
    def allele_frequency_threshold = "-f ${params.allele_freq_thres}"
    def column_for_chromosome = "-c ${params.chrom_col}"
    def column_for_gene_name = "g ${params.gene_name_col}"

    """
    export JAVA_OPTS='"-Xms${task.memory.toMega()/4}m" "-Xmx${task.memory.toGiga()}g" "-Dsamjdk.reference_fasta=${fasta}"'

    vardict-java \\
        ${args} \\
        ${input} \\
        ${min_num_variant_reads} \\
        ${column_for_region_start} \\
        ${column_for_region_end} \\
        ${allele_frequency_threshold} \\
        ${column_for_chromosome} \\
        ${colum_for_gene_name} \\
        -th ${task.cpus} \\
        -G ${fasta} \\
        ${bed} \\
    | ${filter} \\
    | ${convert_to_vcf} \\
        ${args2} \\
    > ${prefix}.vcf 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
        var2vcf_valid.pl: \$( var2vcf_valid.pl -h | sed '2!d;s/.* //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: '-c 1 -S 2 -E 3'
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vardict-java: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
        var2vcf_valid.pl: \$( var2vcf_valid.pl -h | sed '2!d;s/.* //' )
    END_VERSIONS
    """
}
