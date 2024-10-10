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

include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { PREPARE_INPUTS         } from '../subworkflows/local/prepare_inputs'
include { VARDICT_PROCESSING     } from '../subworkflows/local/vardict_processing'
include { MUTECT1_PROCESSING     } from '../subworkflows/local/mutect1_processing'
include { MUTECT2_PROCESSING     } from '../subworkflows/local/mutect2_processing'
include { BCFTOOLS_CONCAT_WITH_MUTECT     } from '../subworkflows/local/bcftools_concat_with_mutect'
include { GENOME_NEXUS } from '../subworkflows/msk/genome_nexus/main'
include { TRACEBACK } from '../subworkflows/msk/traceback/main'
include { BCFTOOLS_ANNOTATE as ANNOTATE_WITH_MUTECT1 } from '../modules/local/bcftools/annotate/main'
include { BCFTOOLS_ANNOTATE as ANNOTATE_WITH_VARDICT } from '../modules/local/bcftools/annotate/main'
include { VCF2MAF } from '../modules/local/vcf2maf'
include { PVMAF_TAGTRACEBACK } from '../modules/msk/pvmaf/tagtraceback'
include { MAF_PROCESSING     } from '../modules/local/maf_processing/main'
include { ACCESS_FILTERS } from '../modules/local/access_filters/main'
include { TAG_BY_VARIANT_ANNOTATION } from '../modules/local/tag_by_variant_annotation/main'




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
    input_ss
    aux_bams_ss
    sample_id_names
    sample_order_file
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

    canonical_bed = Channel.from(params.canonical_bed)
    fasta_ref = params.fasta
    fasta_index = params.fai
    fasta_dict = params.dict
    rules_file = Channel.fromPath(params.rules_json)
    mutect1_header_file = Channel.fromPath("${projectDir}/tests/resources/v1.0/mutect_annotate_concat_header.txt")
    mutect2_header_file = Channel.fromPath("${projectDir}/tests/resources/v1.0/mutect2_annotate_concat_header.txt")
    vardict_header_file = Channel.fromPath("${projectDir}/tests/resources/v1.0/vardict_annotate_concat_header.txt")
    blocklist = Channel.fromPath(params.blocklist)
    canonical_tx_ref = Channel.fromPath(params.canonical_tx_ref)
    hotspots = Channel.fromPath(params.hotspots)




    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    PREPARE_INPUTS( case_bams,sample_id_names,Channel.from(params.target_bed) )

    target_bed_file = PREPARE_INPUTS.out.target_bed_file

    // processing subworkflows for each variant caller
    VARDICT_PROCESSING( sample_id_names,duplex_bams,fasta_ref,fasta_index,fasta_dict,target_bed_file )
    vardict_concat_vcf_isolated = VARDICT_PROCESSING.out.vardict_concat_vcf_isolated
    vardict_index = VARDICT_PROCESSING.out.vardict_index

    MUTECT1_PROCESSING( target_bed_file,fasta_ref,fasta_index,fasta_dict,duplex_bams,sample_id_names,sample_order_file )
    mutect1_vcf_isolated = MUTECT1_PROCESSING.out.mutect1_vcf_isolated
    mutect1_index = MUTECT1_PROCESSING.out.mutect1_index

    MUTECT2_PROCESSING( target_bed_file,fasta_ref,fasta_index,fasta_dict,duplex_bams,sample_id_names,sample_order_file )


    // // concatenation of VarDict and MuTect VCFs
    BCFTOOLS_CONCAT_WITH_MUTECT( sample_id_names,vardict_concat_vcf_isolated,mutect1_vcf_isolated,vardict_index,mutect1_index )
    sample_plus_final_concat_vcf = BCFTOOLS_CONCAT_WITH_MUTECT.out.sample_plus_final_concat_vcf
    sample_plus_final_concat_vcf.map{ meta,vcf -> vcf}.set{ mutect_vardict_concat_vcf }

    mutect1_vcf_isolated.map{ meta,vcf -> vcf}.set{ mutect1_norm_sorted_vcf_isolated }
    sample_id_names.combine(mutect_vardict_concat_vcf).combine(mutect1_norm_sorted_vcf_isolated).set{ input_for_bcftools_annotate }



    // // // // annotate the concatenated VarDict/MuTect VCF against MuTect original VCF
    ANNOTATE_WITH_MUTECT1( input_for_bcftools_annotate,mutect1_header_file )
    annotated_with_mutect1_vcf = ANNOTATE_WITH_MUTECT1.out.vcf
    annotated_with_mutect1_vcf.combine(vardict_concat_vcf_isolated).set{ input_for_bcftools_annotate2 }
    ANNOTATE_WITH_VARDICT( input_for_bcftools_annotate2,vardict_header_file )
    annotated_vcf = ANNOTATE_WITH_VARDICT.out.vcf


    if (params.annotator == 'genomenexus') {
        println "User has specified genomenexus for annotation software flag. Proceeding with Genome Nexus Subworkflow..."
        GENOME_NEXUS( annotated_vcf )
        input_maf = GENOME_NEXUS.out.maf
    }
    else if (params.annotator == 'vcf2maf') {
        println "User has specified vcf2maf for annotation software flag. Proceeding with PERL vcf2maf module..."
        annotated_vcf.combine(Channel.from(ref_fasta)).set{ inputs_for_perlvcf2maf }
        VCF2MAF( inputs_for_perlvcf2maf )
        input_maf = VCF2MAF.out.maf
    }
    input_maf.map{ meta,maf -> maf}.set{ test_maf_only }



    // // // // // // traceback subworkflow
    input_maf.map{ meta,maf -> tuple([patient: 'test',id:"${meta.case_id}.${meta.control_id}.combined-variants"],maf)}.set{ mafs }
    case_bams_for_traceback.mix(control_bams_for_traceback).mix(aux_bams).mix(normal_bams).set{ bams }
    TRACEBACK( bams, mafs, fasta_ref, fasta_index )
    PVMAF_TAGTRACEBACK(TRACEBACK.out.genotyped_maf, [params.input, params.aux_bams])
    genotyped_maf = PVMAF_TAGTRACEBACK.out.maf


    // // // // // maf_processing module (tag by rules)
    MAF_PROCESSING( genotyped_maf, rules_file, hotspots)
    tagged_maf = MAF_PROCESSING.out.maf
    tagged_maf.map{ meta,maf -> maf}.set{tagged_maf_only}


    // // // // // access filters
    sample_id_names.combine(tagged_maf_only).combine(test_maf_only).combine(blocklist).set{ inputs_for_access_filters }
    ACCESS_FILTERS( inputs_for_access_filters )
    access_filtered_maf = ACCESS_FILTERS.out.filtered_maf
    access_filtered_condensed_maf = ACCESS_FILTERS.out.condensed_filtered_maf


    // // // // // // mpath loading script module
    TAG_BY_VARIANT_ANNOTATION( access_filtered_maf,canonical_tx_ref )

    annotated_exonic_maf_file = TAG_BY_VARIANT_ANNOTATION.out.annotated_exonic
    annotated_silent_file = TAG_BY_VARIANT_ANNOTATION.out.annotated_silent
    annotated_nonpanel_exonic_file = TAG_BY_VARIANT_ANNOTATION.out.annotated_nonpanel_exonic
    annotated_nonpanel_silent_file = TAG_BY_VARIANT_ANNOTATION.out.annotated_nonpanel_silent
    annotated_dropped_file = TAG_BY_VARIANT_ANNOTATION.out.annotated_dropped


    // // //Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    input_ss
    aux_bams_ss
    access_filtered_maf
    access_filtered_condensed_maf
    annotated_exonic_maf_file
    annotated_silent_file
    annotated_nonpanel_exonic_file
    annotated_nonpanel_silent_file
    annotated_dropped_file



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
