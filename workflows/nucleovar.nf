/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { BCFTOOLS_VARDICT     } from '../subworkflows/local/bcftools_vardict'
include { BCFTOOLS_MUTECT     } from '../subworkflows/local/bcftools_mutect'
include { BCFTOOLS_REHEADER as MUTECT1_REHEADER } from '../modules/local/bcftools/reheader/main'
include { BCFTOOLS_REHEADER as MUTECT2_REHEADER } from '../modules/local/bcftools/reheader/main'
include { BCFTOOLS_ANNOTATE } from '../modules/local/bcftools/annotate/main'
include { CALL_VARIANTS_CASECONTROL     } from '../subworkflows/local/call_variants_casecontrol'
include { BCFTOOLS_CONCAT_VARDICTS     } from '../subworkflows/local/bcftools_concat_vardicts'
include { MODULE4     } from '../subworkflows/local/module4'
include { GUNZIP_FILES     } from '../modules/local/gunzip_files'
include { MUTECT1        } from '../modules/msk/mutect1'
include { MUTECT2     } from '../modules/local/mutect2/main'
include { SAMTOOLS_SORT     } from '../modules/local/samtools/sort/main.nf'
include { BEDTOOLS_GENOMECOV } from '../modules/local/bedtools/genomecov/main'
include { BEDTOOLS_MERGE } from '../modules/local/bedtools/merge/main'
include { MUTECT_FILTER     } from '../modules/local/mutect_filter'
include { BCFTOOLS_CONCAT_WITH_MUTECT     } from '../subworkflows/local/bcftools_concat_with_mutect'
include { GENOME_NEXUS } from '../subworkflows/msk/genome_nexus/main'
include { TRACEBACK } from '../subworkflows/msk/traceback/main'
include { PVMAF_TAGTRACEBACK } from '../modules/msk/pvmaf/tagtraceback'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NUCLEOVAR {

    take:
    samplesheet // channel: samplesheet read in from --input
    sample_id_names
    standard_bams
    case_bams
    control_bams
    duplex_bams
    case_bams_for_traceback
    control_bams_for_traceback
    samplesheets_for_traceback
    aux_bams
    normal_bams


    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    // sanity check to see if tumor and normal samples are included-- kicks off case control method of running pipeline
    def target = (params.target_bed != "") ? "--target_bed $params.target_bed" : ""



    canonical_bed = Channel.from(params.canonical_bed)
    fasta_ref = params.fasta
    fasta_index = params.fai
    fasta_dict = params.dict

    // if (params.target_bed == null) {
    //     println "Target BED file not provided. Will generate BED file for TUMOR/CASE sample."
    //     sample_id_names
    //         .map { meta -> [meta.case_id,meta.control_id] }
    //         .map { items -> items.join('\n') }
    //         .view { data -> new File('sample_order.txt').text = data }
    //     sample_order_file = Channel.fromPath('sample_order.txt')
    //     //invoke the samtools sort and bedtools modules to generate a BED file for the tumor sample
    //     case_bams.map{ bam,bai -> bam}.set{ case_bam_only }

    //     sample_id_names.combine(case_bam_only).set{ input_tumorbam_for_bedtools}
    //     SAMTOOLS_SORT( input_tumorbam_for_bedtools )
    //     sorted_tumor_bam = SAMTOOLS_SORT.out.sorted_bam
    //     target_bed_file = BEDTOOLS_GENOMECOV( sorted_tumor_bam, Channel.from(fasta_index) )
    // }
    // else {
    //     println "Target BED file is provided."
    //     target_bed_file = params.target_bed
    // }

    // CALL_VARIANTS_CASECONTROL (sample_id_names,duplex_bams,Channel.from(fasta_ref),Channel.from(fasta_index),Channel.from(fasta_dict),target_bed_file)
    // vardict_filtered_vcfs = CALL_VARIANTS_CASECONTROL.out.vardict_filtered_vcfs

    // vardict_filtered_vcfs
    //     .map{ standard_vcf,complexvar_vcf -> standard_vcf}
    //     .set{ vardict_filtered_vcf_standard }

    // vardict_filtered_vcfs
    //     .map{ standard_vcf,complexvar_vcf -> complexvar_vcf}
    //     .set{ vardict_filtered_vcf_complexvar }



    // target_bed_file
    //     .combine(Channel.from(fasta_ref))
    //     .combine(Channel.from(fasta_index))
    //     .combine(Channel.from(fasta_dict))
    //     .set{ input2_for_mutect }


    // duplex_bams
    //     .map{ meta,control_bam,control_bai,case_bam,case_bai ->
    //         tuple(case_bam,control_bam,case_bai,control_bai)}
    //     .set{ bams_for_mutect }

    // sample_id_names
    //     .combine(bams_for_mutect)
    //     .set{ input1_for_mutect }



    // MUTECT1(input1_for_mutect,input2_for_mutect)
    // MUTECT2( input1_for_mutect,input2_for_mutect )

    // mutect1_vcf = MUTECT1.out.mutect_vcf
    // mutect1_txt = MUTECT1.out.standard_mutect_output
    // mutect2_vcf = MUTECT2.out.mutect2_vcf


    // mutect1_vcf.combine(sample_order_file).set{ input_for_mutect1_reheader }
    // mutect2_vcf.combine(sample_order_file).set{ input_for_mutect2_reheader }

    // MUTECT1_REHEADER( input_for_mutect1_reheader )
    // MUTECT2_REHEADER( input_for_mutect2_reheader )

    // mutect1_ordered_vcf = MUTECT1_REHEADER.out.sample_reordered_vcf

    // mutect1_ordered_vcf.combine(mutect1_txt).set{ input_for_mutect_filter }

    // MUTECT_FILTER(input_for_mutect_filter,Channel.from(fasta_ref))
    // mutect_filtered_vcf = MUTECT_FILTER.out.mutect_filtered_vcf



    // sample_id_names.combine(vardict_filtered_vcf_standard).set{ standard_vcf_for_bcftools }
    // sample_id_names.combine(vardict_filtered_vcf_complexvar).set{ complexvar_vcf_for_bcftools }
    // BCFTOOLS_VARDICT( vardict_filtered_vcf_complexvar,vardict_filtered_vcf_standard,Channel.from(fasta_ref),Channel.from(fasta_index) )

    // vardict_concat_vcf = BCFTOOLS_VARDICT.out.vardict_concat_vcf
    // vardict_index = BCFTOOLS_VARDICT.out.vardict_index




    // //     // // // // temp testing mutect filtered vcf (permission error in mutect filter)
    // //mutect_filtered_vcf = Channel.fromPath("/home/naidur/snvs_indels/C-PR83CF-L004-d04_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.mutect_filtered.vcf")



    // BCFTOOLS_MUTECT( mutect_filtered_vcf,Channel.from(fasta_ref),Channel.from(fasta_index) )
    // mutect_vcf = BCFTOOLS_MUTECT.out.standard_norm_sorted_vcf
    // mutect_index = BCFTOOLS_MUTECT.out.mutect_index


    // vardict_concat_vcf.map{ id,vcf -> vcf}.set{ vardict_concat_vcf_isolated }
    // mutect_vcf.map{ id,vcf -> vcf}.set{ mutect_vcf_isolated }

    // BCFTOOLS_CONCAT_WITH_MUTECT( sample_id_names,vardict_concat_vcf_isolated,mutect_vcf_isolated,vardict_index,mutect_index )
    // sample_plus_final_concat_vcf = BCFTOOLS_CONCAT_WITH_MUTECT.out.sample_plus_final_concat_vcf


    // sample_plus_final_concat_vcf.map{ meta,vcf -> vcf}.set{ mutect_vardict_concat_vcf }
    // mutect_vcf.map{ meta,vcf -> vcf}.set{ mutect_norm_sorted_vcf_isolated }
    // sample_id_names.combine(mutect_vardict_concat_vcf).combine(mutect_norm_sorted_vcf_isolated).set{ input_for_bcftools_annotate }

    // //NOTE: turn into argument on CLI
    // header_file = Channel.from("/home/naidur/snvs_indels/nucleovar/tests/resources/v1.0/mutect_annotate_concat_header.txt")


    // BCFTOOLS_ANNOTATE( input_for_bcftools_annotate,header_file )
    // annotated_vcf = BCFTOOLS_ANNOTATE.out.vcf


    // // // Genome nexus subworkflow
    // GENOME_NEXUS( annotated_vcf )

    // input_maf = GENOME_NEXUS.out.maf
    // input_maf.map{ meta,maf -> maf}.set{ test_maf_only }


    //     // traceback subworkflow
    test_maf_only = Channel.fromPath("/work/access/production/data/small_variants/C-PR83CF/C-PR83CF-L004-d04/current/C-PR83CF-L004-d04.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots.maf")
    mafs = Channel.from([patient:'test',id:"C-PR83CF-L004-d04.DONOR22-TP.combined-variants"]).merge(test_maf_only)

    case_bams_for_traceback.mix(control_bams_for_traceback).mix(aux_bams).mix(normal_bams).set{ bams }

    TRACEBACK( bams, mafs, fasta_ref, fasta_index )
    PVMAF_TAGTRACEBACK(TRACEBACK.out.genotyped_maf, [params.input, params.aux_bams])
    genotyped_maf = PVMAF_TAGTRACEBACK.out.maf
    // // maf_processing module (tag by rules)
    // MAF_PROCESSING( genotyped_maf, rules_file )
    // tagged_maf = MAF_PROCESSING.out.maf
    // // access filters
    // ACCESS_FILTERS( tagged_maf )
    // access_filtered_maf = ACCESS_FILTERS.out.maf
    // // mpath loading script module
    // TAG_BY_VARIANT_ANNOTATION( access_filtered_maf )



    //Collate and save software versions

    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:
    genotyped_maf
    versions       = ch_versions                 // channel: [ path(versions.yml) ]


}

def create_std_filename( it ) {

    def finalFileName = "vcf_standard.vcf.gz"

    file(it).moveTo(file(finalFileName))
}

def create_mutect_filename( it ) {

    def finalFileName = "vcf_mutect.vcf.gz"

    file(it).moveTo(file(finalFileName))
}

def create_complexvar_filename( it ) {

    def finalFileName = "vcf_complexvar.vcf.gz"

    file(it).moveTo(file(finalFileName))
}


def create_duplex_bams_channel(LinkedHashMap row) {
    // create meta map

    // set bams variable
    if (row.type.isEmpty()) {
        error "Sample ID is missing for one or more entries. Please double check samplesheet before running the workflow."
    } else {
        if (row.type == "case") {
            case_id = row.type
        if (row.type == "control") {
            control_id = row.type
        }
        sample_ids = [case_id,control_id]
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
