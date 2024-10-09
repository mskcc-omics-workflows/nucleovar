//
// Check input samplesheet and get read channels
//

include { MUTECT2     } from '../../modules/local/mutect2/main'
include { BCFTOOLS_REHEADER as MUTECT2_REHEADER } from '../../modules/local/bcftools/reheader/main'
include { MUTECT2_FILTER     } from '../../modules/local/mutect2_filter/main'
include { BCFTOOLS_MUTECT as BCFTOOLS_MUTECT2     } from '../../subworkflows/local/bcftools_mutect'


workflow MUTECT2_PROCESSING {
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
        .set{ input2_for_mutect2 }

    duplex_bams
        .map{ meta,control_bam,control_bai,case_bam,case_bai ->
            tuple(case_bam,control_bam,case_bai,control_bai)}
        .set{ bams_for_mutect2 }

    sample_id_names
        .combine(bams_for_mutect2)
        .set{ input1_for_mutect2 }

    MUTECT2(input1_for_mutect2,input2_for_mutect2)

    mutect2_vcf = MUTECT2.out.mutect2_vcf

    
    mutect2_vcf.combine(sample_order_file).set{ input_for_mutect2_reheader }


    MUTECT2_REHEADER( input_for_mutect2_reheader )

    mutect2_ordered_vcf = MUTECT2_REHEADER.out.sample_reordered_vcf

    mutect2_ordered_vcf.view()
    // MUTECT_FILTER(input_for_mutect_filter,Channel.from(fasta_ref))
    // mutect_filtered_vcf = MUTECT_FILTER.out.mutect_filtered_vcf

    // BCFTOOLS_MUTECT( mutect_filtered_vcf,Channel.from(fasta_ref),Channel.from(fasta_index) )
    // mutect_vcf = BCFTOOLS_MUTECT.out.standard_norm_sorted_vcf
    // mutect1_index = BCFTOOLS_MUTECT.out.mutect_index

    // mutect_vcf.map{ id,vcf -> vcf}.set{ mutect1_vcf_isolated }

    // mutect1_vcf_isolated.view()

    // emit:
    // mutect1_vcf_isolated 
    // mutect1_index
}





