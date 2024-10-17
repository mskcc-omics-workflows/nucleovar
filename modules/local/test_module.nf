process TEST_MODULE {

    container "ghcr.io/msk-access/bcftools:1.15.1"

    input:
    val(samplename)
    path(vardict_concat_vcf)
    path(mutect_vcf)


    output:
    path("completed.txt")

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''


    """
    bcftools concat --output ${samplename}.norm.vcf.gz ${vardict_concat_vcf} ${mutect_vcf}


    """
}
