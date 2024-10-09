//
// Check input samplesheet and get read channels
//

include { CALL_VARIANTS_CASECONTROL     } from '../../subworkflows/local/call_variants_casecontrol'
include { BCFTOOLS_VARDICT     } from '../../subworkflows/local/bcftools_vardict'

workflow VARDICT_PROCESSING {
    take:
    sample_id_names
    duplex_bams
    fasta_ref
    fasta_index
    fasta_dict
    target_bed_file

    main:

    CALL_VARIANTS_CASECONTROL( sample_id_names,duplex_bams,Channel.from(fasta_ref),Channel.from(fasta_index),Channel.from(fasta_dict),target_bed_file )

    vardict_filtered_vcf_standard = CALL_VARIANTS_CASECONTROL.out.std_vardict_vcf
    vardict_filtered_vcf_complexvar = CALL_VARIANTS_CASECONTROL.out.complex_variants_vardict_vcf

    BCFTOOLS_VARDICT( vardict_filtered_vcf_complexvar,vardict_filtered_vcf_standard,Channel.from(fasta_ref),Channel.from(fasta_index) )

    vardict_concat_vcf = BCFTOOLS_VARDICT.out.vardict_concat_vcf
    vardict_index = BCFTOOLS_VARDICT.out.vardict_index

    vardict_concat_vcf.map{ id,vcf -> vcf}.set{ vardict_concat_vcf_isolated }


    emit:
    vardict_concat_vcf_isolated
    vardict_index
}





