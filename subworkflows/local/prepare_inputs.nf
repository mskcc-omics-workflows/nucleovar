//
// Check input samplesheet and get read channels
//

include { SAMTOOLS_SORT     } from '../../modules/local/samtools/sort/main.nf'
include { BEDTOOLS_GENOMECOV } from '../../modules/local/bedtools/genomecov/main'
include { BEDTOOLS_MERGE } from '../../modules/local/bedtools/merge/main'





workflow PREPARE_INPUTS {
    take:
    case_bams
    sample_id_names
    target_bed

    main:

    if (target_bed == null) {
        println "Target BED file not provided. Will generate BED file for TUMOR/CASE sample."

        // //invoke the samtools sort and bedtools modules to generate a BED file for the tumor sample
        case_bams.map{ bam,bai -> bam}.set{ case_bam_only }

        sample_id_names.combine(case_bam_only).set{ input_tumorbam_for_bedtools}
        SAMTOOLS_SORT( input_tumorbam_for_bedtools )
        sorted_tumor_bam = SAMTOOLS_SORT.out.sorted_bam
        target_bed_file = BEDTOOLS_GENOMECOV( sorted_tumor_bam, Channel.from(fasta_index) )
    }
    else {
        println "Target BED file is provided."
        target_bed_file = target_bed
    }

    emit:
    target_bed_file
}





