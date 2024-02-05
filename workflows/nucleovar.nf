/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowNucleovar.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { BCFTOOLS_ALL_VARDICT     } from '../subworkflows/local/bcftools_all_vardict'
include { BCFTOOLS_ALL_VARDICT_COMPLEXVAR     } from '../subworkflows/local/bcftools_all_vardict_complexvar'
include { BCFTOOLS_ALL_MUTECT     } from '../subworkflows/local/bcftools_all_mutect'

include { BCFTOOLS_ALL_LOCAL_VARDICT     } from '../subworkflows/local/bcftools_all_local_vardict'

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

// Info required for completion email and summary
def multiqc_report = []

workflow NUCLEOVAR {

    ch_versions = Channel.empty()
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    
    // read in input samplesheet and if single sample, run single sample workflow (just vardict). if double sample, run vardict and mutect
    
    // CALL_VARIANTS (
    //     params.input,params.bed,params.fasta,params.fai
    // )

    //vardict_filtered_vcf = CALL_VARIANTS.out.vardict_filtered_vcf
    //complex_variants_vardict_filtered_vcf = CALL_VARIANTS.out.complex_variants_vardict_filtered_vcf
    //mutect_filtered_vcf = CALL_VARIANTS.out.mutect_filtered_vcf

    // temporary testing purposes, reading in sample vardict and mutect filtered VCFs
    vardict_filtered_vcf = Channel.fromPath('/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-C-2HXC96-P001-d01_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.vardict.filtered.vcf')
    vardict_complexvar_filtered_vcf = Channel.fromPath('/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-C-2HXC96-P001-d01_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex_STDfilter_complex.vcf')
    mutect_filtered_vcf = Channel.fromPath('/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/C-2HXC96-P001-d01_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.mutect_filter.mutect.vcf')
    ref_fasta = Channel.fromPath('/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/Homo_sapiens_assembly19.fasta')
    ref_fasta_index = Channel.fromPath('/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/Homo_sapiens_assembly19.fasta.fai')

    
    
    BCFTOOLS_ALL_VARDICT( vardict_filtered_vcf,ref_fasta )
    meta_plus_vardict_norm_and_sorted_vcf = BCFTOOLS_ALL_VARDICT.out.meta_plus_vardict_norm_and_sorted_vcf
    vardict_index = BCFTOOLS_ALL_VARDICT.out.vardict_filtered_vcf_and_index

    BCFTOOLS_ALL_VARDICT_COMPLEXVAR( vardict_complexvar_filtered_vcf,ref_fasta )
    meta_plus_vardict_complexvar_norm_and_sorted_vcf = BCFTOOLS_ALL_VARDICT_COMPLEXVAR.out.meta_plus_vardict_complexvar_norm_and_sorted_vcf
    vardict_complexvar_index = BCFTOOLS_ALL_VARDICT_COMPLEXVAR.out.vardict_complexvar_filtered_vcf_and_index
    
    BCFTOOLS_ALL_MUTECT( mutect_filtered_vcf,ref_fasta )
    meta_plus_vardict_complexvar_norm_and_sorted_vcf = BCFTOOLS_ALL_MUTECT.out.meta_plus_vardict_complexvar_norm_and_sorted_vcf 
    mutect_index = BCFTOOLS_ALL_MUTECT.out.mutcet_filtered_vcf_and_index

    

    // BCF CONCATENATE WORKFLOW FOR VARDICTS
    meta_plus_vardict_norm_and_sorted_vcf
        .map{ sample,vcf -> vcf}
        .set{ vardict_regular_vcf }

    meta_plus_vardict_complexvar_norm_and_sorted_vcf
        .map{ sample,vcf -> vcf}
        .set{ vardict_complexvar_vcf }
    

    // BCF CONCATENATE WORKFLOW FOR COMBINED VARDICTS + MUTECT 

    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
    }


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
