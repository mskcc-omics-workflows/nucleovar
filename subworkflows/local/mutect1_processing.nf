//
// Check input samplesheet and get read channels
//

include { MUTECT1        } from '../../modules/msk/mutect1'
include { BCFTOOLS_REHEADER as MUTECT1_REHEADER } from '../../modules/local/bcftools/reheader/main'
include { MUTECT_FILTER     } from '../../modules/local/mutect_filter'
include { BCFTOOLS_MUTECT as BCFTOOLS_MUTECT1     } from '../../subworkflows/local/bcftools_mutect'


workflow MUTECT1_PROCESSING {
    take:
    target_bed_file
    fasta_ref
    fasta_index
    fasta_dict 
    duplex_bams 
    sample_id_names 
    sample_order_file

    main:

    target_bed_file
        .combine(Channel.from(fasta_ref))
        .combine(Channel.from(fasta_index))
        .combine(Channel.from(fasta_dict))
        .set{ input2_for_mutect }

    duplex_bams
        .map{ meta,control_bam,control_bai,case_bam,case_bai ->
            tuple(case_bam,control_bam,case_bai,control_bai)}
        .set{ bams_for_mutect }

    sample_id_names
        .combine(bams_for_mutect)
        .set{ input1_for_mutect }

    MUTECT1(input1_for_mutect,input2_for_mutect)

    mutect1_vcf = MUTECT1.out.mutect_vcf
    mutect1_txt = MUTECT1.out.standard_mutect_output
    
    mutect1_vcf.combine(sample_order_file).set{ input_for_mutect1_reheader }


    MUTECT1_REHEADER( input_for_mutect1_reheader )

    mutect1_ordered_vcf = MUTECT1_REHEADER.out.sample_reordered_vcf
    mutect1_txt.map{ meta,txt -> txt}.set{ mutect1_txt_only }

    mutect1_ordered_vcf.combine(mutect1_txt_only).set{ input_for_mutect_filter }

    
    MUTECT_FILTER(input_for_mutect_filter,Channel.from(fasta_ref))
    mutect_filtered_vcf = MUTECT_FILTER.out.mutect_filtered_vcf

    BCFTOOLS_MUTECT1( mutect_filtered_vcf,Channel.from(fasta_ref),Channel.from(fasta_index) )
    mutect_vcf = BCFTOOLS_MUTECT1.out.standard_norm_sorted_vcf
    mutect1_index = BCFTOOLS_MUTECT1.out.mutect_index

    mutect_vcf.map{ id,vcf -> vcf}.set{ mutect1_vcf_isolated }

    

    emit:
    mutect1_vcf_isolated 
    mutect1_index
}





