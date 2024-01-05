process MUTECT {
    label 'process_high'
    container 'ghcr.io/msk-access/test:latest'

    input:
    tuple val(meta), path(bams), path(bais)
    tuple val(bed), val(fasta), val(fai)

    output:
    tuple val(meta), path("*.vcf"), emit: vcf

    when:
    task.ext.when == null || task.ext.when

    script:
    def mutectargs = task.ext.mutectargs ?: ''
    def resourceargs = task.ext.resourceargs ?: ''

    """
    java -Xmx28g -Xms256m -XX:-UseGCOverheadLimit -jar /usr/bin/mutect.jar -T MuTect \
    --reference ${fasta} \
    --normal_sample_name ${control_sample_name} \
    --tumor_sample_name ${case_sample_name} \
    --input_file_normal ${control_bam} \
    --input_file_tumor ${case_bam} \
    ${mutectargs} \
    ${resourceargs}

    """
}