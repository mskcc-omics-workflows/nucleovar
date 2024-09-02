process BEDTOOLS_GENOMECOV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bedtools_coreutils:ba273c06a3909a15':
        'community.wave.seqera.io/library/bedtools_coreutils:a623c13f66d5262b' }"

    input:
    tuple val(meta), path(sorted_tumor_bam)
    path(fasta_index)

    output:
    path("*.bed"), emit: target_bed_file
    //path  "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bedtools genomecov -ibam ${sorted_tumor_bam} -bg > target.bedgraph
    bedtools merge -i target.bedgraph > ${meta.case_id}_target.bed

    """

    stub:
    """
    touch ${meta.case_id}_target.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
