

//import bcftools subworkflow for vardict filtered regular
//include { GUNZIP } from '../../modules/nf-core/gunzip/main'
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main'
include { GUNZIP_FILES } from '../../modules/local/gunzip_files'



workflow BCFTOOLS_CONCAT_VARDICTS {
    take:
    vardict_norm_and_sorted_vcf_standard
    vardict_index
    vardict_norm_and_sorted_vcf_complexvar
    vardict_complexvar_index



    main:

    vardict_norm_and_sorted_vcf_standard.map{ sample,vcf -> sample}.set{ sampleid_for_bcftools }
    //vardict_norm_and_sorted_vcf_standard.map{ sample,vcf -> vcf}.set{ standard_vcf_for_bcftools }
    //vardict_norm_and_sorted_vcf_complexvar.map{ sample,vcf -> vcf}.set{ complexvar_vcf_for_bcftools }




    sampleid_for_bcftools
        .combine(standard_vcf_for_bcftools)
        .combine(complexvar_vcf_for_bcftools)
        .combine(vardict_index)
        .combine(vardict_complexvar_index)
        .map{ sample,vcf1,vcf2,index1,index2 -> [sample,[first:vcf1,second:vcf2],[first:index1,second:index2]] }
        .set{ inputs_for_bcftools_concat }



    BCFTOOLS_CONCAT( inputs_for_bcftools_concat )

    // sample_plus_vardict_concat_vcf = BCFTOOLS_CONCAT.out.vcf


    // emit:
    // sample_plus_vardict_concat_vcf

}




