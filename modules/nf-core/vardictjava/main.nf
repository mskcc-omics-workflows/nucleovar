process VARDICTJAVA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/vardict-java.yml"

    input:
    tuple val(meta), path(bams), path(bais), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions


    when:
    task.ext.when == null || task.ext.when


    script:
    def args = task.ext.args ?: '-c 1 -S 2 -E 3'
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def somatic = bams instanceof List && bams.size() == 2 ? true : false
    def input = somatic ? "-b \"${bams[0]}|${bams[1]}\"" : "-b ${bams}"
    def filter = somatic ? "testsomatic.R" : "teststrandbias.R"
    def convert_to_vcf = somatic ? "var2vcf_paired.pl" : "var2vcf_valid.pl"
    
    """
    export JAVA_OPTS='"-Xms${task.memory.toMega()/4}m" "-Xmx${task.memory.toGiga()}g" "-Dsamjdk.reference_fasta=${fasta}"'
    vardict-java \\
        ${args} \\
        ${input} \\
        -th ${task.cpus} \\
        -G ${fasta} \\
        ${bed} \\
    | ${filter} \\
    | ${convert_to_vcf} \\
        ${args2} \\
        "C-2HXC96-P001-d01_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex|DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex" \\
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


// export JAVA_OPTS='"-Xms${task.memory.toMega()/4}m" "-Xmx${task.memory.toGiga()}g" "-Dsamjdk.reference_fasta=${fasta}"'
//     vardict-java \\
//         ${args} \\
//         ${input} \\
//         -th ${task.cpus} \\
//         -G ${fasta} \\
//         ${bed} \\
//     | ${filter} \\
//     | ${convert_to_vcf} \\
//         ${args2} \\
//         "C-2HXC96-P001-d01_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex|DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex" \\
//     > ${prefix}.vcf

//     cat <<-END_VERSIONS > versions.yml
//     "${task.process}":
//         vardict-java: \$( realpath \$( command -v vardict-java ) | sed 's/.*java-//;s/-.*//' )
//         var2vcf_valid.pl: \$( var2vcf_valid.pl -h | sed '2!d;s/.* //' )
//     END_VERSIONS