

//import bcftools subworkflow for mutect filtered regular
//include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { BCFTOOLS_CONCAT2 as BCFTOOLS_CONCAT_W_MUTECT } from '../../modules/local/bcftools_concat2/main'
include { GUNZIP_FILES } from '../../modules/local/gunzip_files'



workflow BCFTOOLS_CONCAT_WITH_MUTECT {
    take:
    samplename
    vardict_concat_vcf
    mutect_vcf
    vardict_index
    mutect_index



    main:

    samplename
        .combine(vardict_concat_vcf)
        .combine(mutect_vcf)
        .combine(vardict_index)
        .combine(mutect_index)
        .map{ sample,vcf1,vcf2,index1,index2 -> [sample,vcf1,vcf2] }
        .set{ inputs_for_bcftools_concat }



    BCFTOOLS_CONCAT_W_MUTECT( inputs_for_bcftools_concat )

    sample_plus_final_concat_vcf = BCFTOOLS_CONCAT_W_MUTECT.out.concat_vcf


    emit:
    sample_plus_final_concat_vcf

}




