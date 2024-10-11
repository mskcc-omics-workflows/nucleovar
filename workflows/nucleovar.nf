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
include { GUNZIP_FILES     } from '../modules/local/gunzip_files'
include { MUTECT1        } from '../modules/msk/mutect1'
include { MUTECT2     } from '../modules/local/mutect2/main'
include { SAMTOOLS_SORT     } from '../modules/local/samtools/sort/main.nf'
include { BEDTOOLS_GENOMECOV } from '../modules/local/bedtools/genomecov/main'
include { BEDTOOLS_MERGE } from '../modules/local/bedtools/merge/main'
include { MUTECT_FILTER     } from '../modules/local/mutect_filter'
include { BCFTOOLS_CONCAT_WITH_MUTECT     } from '../subworkflows/local/bcftools_concat_with_mutect'
include { VCF2MAF } from '../modules/local/vcf2maf'
include { GENOME_NEXUS } from '../subworkflows/msk/genome_nexus/main'
include { TRACEBACK } from '../subworkflows/msk/traceback/main'
include { PVMAF_TAGTRACEBACK } from '../modules/msk/pvmaf/tagtraceback'
include { MAF_PROCESSING     } from '../modules/local/maf_processing/main'
include { ACCESS_FILTERS } from '../modules/local/access_filters/main'
include { TAG_BY_VARIANT_ANNOTATION } from '../modules/local/tag_by_variant_annotation/main'
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
    samplesheet
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

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    // sanity check to see if tumor and normal samples are included-- kicks off case control method of running pipeline
    def target = (params.target_bed != "") ? "--target_bed $params.target_bed" : ""

    canonical_bed = Channel.from(params.canonical_bed)
    fasta_ref = params.fasta
    fasta_index = params.fai
    fasta_dict = params.dict
    rules_file = Channel.fromPath(params.rules_json)
    header_file = Channel.fromPath(params.header_file)
    blocklist = Channel.fromPath(params.blocklist)
    canonical_tx_ref = Channel.fromPath(params.canonical_tx_ref)
    hotspots = Channel.fromPath(params.hotspots)


    if (params.target_bed == null) {
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
        target_bed_file = Channel.from(params.target_bed)
    }

    CALL_VARIANTS_CASECONTROL (sample_id_names,duplex_bams,Channel.from(fasta_ref),Channel.from(fasta_index),Channel.from(fasta_dict),target_bed_file)

    vardict_filtered_vcf_standard = CALL_VARIANTS_CASECONTROL.out.std_vardict_vcf
    vardict_filtered_vcf_complexvar = CALL_VARIANTS_CASECONTROL.out.complex_variants_vardict_vcf


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


    // run mutect1 and mutect2 variant callers
    MUTECT1(input1_for_mutect,input2_for_mutect)
    // MUTECT2( input1_for_mutect,input2_for_mutect )

    mutect1_vcf = MUTECT1.out.mutect_vcf
    mutect1_txt = MUTECT1.out.standard_mutect_output
    // mutect2_vcf = MUTECT2.out.mutect2_vcf


    //mutect1_vcf.combine(sample_order_file).set{ input_for_mutect1_reheader }
    // mutect2_vcf.combine(sample_order_file).set{ input_for_mutect2_reheader }
    
    // standardizes the order of samples printed to output VCF in all variant callers (matches what is there for VarDict)
    MUTECT1_REHEADER( mutect1_vcf )
    // MUTECT2_REHEADER( input_for_mutect2_reheader )

    mutect1_ordered_vcf = MUTECT1_REHEADER.out.sample_reordered_vcf
    mutect1_txt.map{ meta,txt -> txt}.set{ mutect1_txt_only }

    mutect1_ordered_vcf.combine(mutect1_txt_only).set{ input_for_mutect_filter }

    // filtering variant callers
    MUTECT_FILTER(input_for_mutect_filter,Channel.from(fasta_ref))
    mutect_filtered_vcf = MUTECT_FILTER.out.mutect_filtered_vcf



    sample_id_names.combine(vardict_filtered_vcf_standard).set{ standard_vcf_for_bcftools }
    sample_id_names.combine(vardict_filtered_vcf_complexvar).set{ complexvar_vcf_for_bcftools }

    // bcftools suite subworkflow for VarDict VCFs
    BCFTOOLS_VARDICT( vardict_filtered_vcf_complexvar,vardict_filtered_vcf_standard,Channel.from(fasta_ref),Channel.from(fasta_index) )

    vardict_concat_vcf = BCFTOOLS_VARDICT.out.vardict_concat_vcf
    vardict_index = BCFTOOLS_VARDICT.out.vardict_index

    // bcftools suite subworkflow for MuTect VCFs
    BCFTOOLS_MUTECT( mutect_filtered_vcf,Channel.from(fasta_ref),Channel.from(fasta_index) )
    mutect_vcf = BCFTOOLS_MUTECT.out.standard_norm_sorted_vcf
    mutect_index = BCFTOOLS_MUTECT.out.mutect_index


    vardict_concat_vcf.map{ id,vcf -> vcf}.set{ vardict_concat_vcf_isolated }
    mutect_vcf.map{ id,vcf -> vcf}.set{ mutect_vcf_isolated }

    // concatenation of VarDict and MuTect VCFs
    BCFTOOLS_CONCAT_WITH_MUTECT( sample_id_names,vardict_concat_vcf_isolated,mutect_vcf_isolated,vardict_index,mutect_index )
    sample_plus_final_concat_vcf = BCFTOOLS_CONCAT_WITH_MUTECT.out.sample_plus_final_concat_vcf


    sample_plus_final_concat_vcf.map{ meta,vcf -> vcf}.set{ mutect_vardict_concat_vcf }
    mutect_vcf.map{ meta,vcf -> vcf}.set{ mutect_norm_sorted_vcf_isolated }
    sample_id_names.combine(mutect_vardict_concat_vcf).combine(mutect_norm_sorted_vcf_isolated).set{ input_for_bcftools_annotate }

    // annotate the concatenated VarDict/MuTect VCF against MuTect original VCF
    BCFTOOLS_ANNOTATE( input_for_bcftools_annotate,header_file )
    annotated_vcf = BCFTOOLS_ANNOTATE.out.vcf


    if (params.annotator == 'genomenexus') {
        println "User has specified genomenexus for annotation software flag. Proceeding with Genome Nexus Subworkflow..."
        GENOME_NEXUS( annotated_vcf )
        input_maf = GENOME_NEXUS.out.maf
    }
    else if (params.annotator == 'vcf2maf') {
        println "User has specified vcf2maf for annotation software flag. Proceeding with PERL vcf2maf module..."

    }

    input_maf.map{ meta,maf -> maf}.set{ test_maf_only }


    // // // // traceback subworkflow
    input_maf.map{ meta,maf -> tuple([patient: 'test',id:"${meta.case_id}.${meta.control_id}.combined-variants"],maf)}.set{ mafs }



    case_bams_for_traceback.mix(control_bams_for_traceback).mix(aux_bams).mix(normal_bams).set{ bams }

    TRACEBACK( bams, mafs, fasta_ref, fasta_index )
    PVMAF_TAGTRACEBACK(TRACEBACK.out.genotyped_maf, [params.input, params.aux_bams])
    genotyped_maf = PVMAF_TAGTRACEBACK.out.maf


    // // // maf_processing module (tag by rules)
    MAF_PROCESSING( genotyped_maf, rules_file, hotspots)
    tagged_maf = MAF_PROCESSING.out.maf
    tagged_maf.map{ meta,maf -> maf}.set{tagged_maf_only}


    // // // access filters

    sample_id_names.combine(tagged_maf_only).combine(test_maf_only).combine(blocklist).set{ inputs_for_access_filters }

    ACCESS_FILTERS( inputs_for_access_filters )
    access_filtered_maf = ACCESS_FILTERS.out.filtered_maf
    access_filtered_condensed_maf = ACCESS_FILTERS.out.condensed_filtered_maf


    // // // // mpath loading script module
    TAG_BY_VARIANT_ANNOTATION( access_filtered_maf,canonical_tx_ref )

    // // //Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
