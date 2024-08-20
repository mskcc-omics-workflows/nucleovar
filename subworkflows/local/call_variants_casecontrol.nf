//
// Check input samplesheet and get read channels
//

include { VARDICTJAVA } from '../../modules/nf-core/vardictjava/main'
include { VARDICT_FILTER } from '../../modules/local/vardict_filter'

workflow CALL_VARIANTS_CASECONTROL {
    take:
    sample_id_names
    duplex_bams
    fasta
    fai
    dict
    bed

    main:

    duplex_bams
        .map{ meta,control_bam,control_bai,case_bam,case_bai -> tuple(meta,[file(control_bam),file(case_bam)],[file(control_bai),file(case_bai)])}
        .combine(bed)
        .set{ vardict_input_set1 }

    sample_id_names
        .combine(fasta)
        .set{ vardict_input_set2 }

    sample_id_names
        .combine(fai)
        .set{ vardict_input_set3 }




    VARDICTJAVA(vardict_input_set1,vardict_input_set2,vardict_input_set3)
    vardict_vcf = VARDICTJAVA.out.vcf
    vardict_vcf.map{ meta,vcf -> file(vcf)}.set{ vardict_vcf_isolated }

    duplex_bams
        .map{ meta,control_bam,control_bai,case_bam,case_bai -> tuple([file(control_bam),file(case_bam)])}
        .set{ bams_for_vardict_filter }

    VARDICT_FILTER( sample_id_names,vardict_vcf_isolated,bams_for_vardict_filter )

    vardict_filtered_vcfs = VARDICT_FILTER.out.filtered_vardict_vcf

    vardict_filter_output_txt = VARDICT_FILTER.out.std_vardict_filter_output


    emit:
    vardict_filtered_vcfs
    vardict_filter_output_txt
}





