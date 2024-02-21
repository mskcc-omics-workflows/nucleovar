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

//include { CALL_VARIANTS     } from '../subworkflows/local/call_variants'
include { BCFTOOLS_ALL_VARDICT     } from '../subworkflows/local/bcftools_all_vardict'
include { BCFTOOLS_ALL_VARDICT_COMPLEXVAR     } from '../subworkflows/local/bcftools_all_vardict_complexvar'
include { BCFTOOLS_ALL_MUTECT     } from '../subworkflows/local/bcftools_all_mutect'
include { CALL_VARIANTS_CASECONTROL     } from '../subworkflows/local/call_variants_casecontrol'
include { BCFTOOLS_CONCAT_VARDICTS     } from '../subworkflows/local/bcftools_concat_vardicts'
include { GUNZIP_FILES     } from '../modules/local/gunzip_files'
include { BCFTOOLS_CONCAT_WITH_MUTECT     } from '../subworkflows/local/bcftools_concat_with_mutect'




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
    def sampleSheet = file(params.input).readLines().collect { it.split(",") }
    def allColumnsHaveValue = sampleSheet.every { row ->
    row.every { cell -> cell.trim() }
}

    if (allColumnsHaveValue) {
        println "Running the SNPs/indels workflow in case-control mode."
        CALL_VARIANTS_CASECONTROL (params.input,params.bed,params.fasta,params.fai)
        vardict_filtered_vcfs = CALL_VARIANTS_CASECONTROL.out.vardict_filtered_vcf

        vardict_filtered_vcfs
            .map{ standard_vcf,complexvar_vcf -> standard_vcf}
            .set{ vardict_filtered_vcf_standard }

        vardict_filtered_vcfs
            .map{ standard_vcf,complexvar_vcf -> complexvar_vcf}
            .set{ vardict_filtered_vcf_complexvar }

        mutect_filtered_vcf = CALL_VARIANTS_CASECONTROL.out.mutect_filtered_vcf
        ref_fasta = CALL_VARIANTS_CASECONTROL.out.genome_fasta_file
        ref_fasta_index = CALL_VARIANTS_CASECONTROL.out.genome_fasta_index_file

        BCFTOOLS_ALL_VARDICT( vardict_filtered_vcf_standard,ref_fasta,ref_fasta_index )

        meta_plus_vardict_norm_and_sorted_vcf_standard = BCFTOOLS_ALL_VARDICT.out.meta_plus_vardict_norm_and_sorted_vcf
        vardict_filtered_index_isolated_standard = BCFTOOLS_ALL_VARDICT.out.vardict_filtered_index_final

        

        BCFTOOLS_ALL_MUTECT( mutect_filtered_vcf,ref_fasta,ref_fasta_index )
        meta_plus_mutect_norm_and_sorted_vcf = BCFTOOLS_ALL_MUTECT.out.meta_plus_mutect_norm_and_sorted_vcf
        mutect_filtered_index_isolated = BCFTOOLS_ALL_MUTECT.out.mutect_filtered_index_final


        meta_plus_vardict_norm_and_sorted_vcf_standard
            .map{ sample,vcf -> sample}
            .set{ vardict_samplename_ch }

        meta_plus_vardict_norm_and_sorted_vcf_standard
            .map{ sample,vcf -> vcf}
            .set{ vardict_norm_and_sorted_vcf_standard }

        meta_plus_mutect_norm_and_sorted_vcf
            .map{ sample,vcf -> vcf}
            .set{ mutect_norm_and_sorted_vcf }


        vardict_norm_and_sorted_vcf_standard
            .map{ create_std_filename(it) }
            .set{ vcf_std_for_bcftools_concat }


        mutect_norm_and_sorted_vcf
            .map{ create_mutect_filename(it) }
            .set{ vcf_mutect_for_bcftools_concat }


        BCFTOOLS_CONCAT_WITH_MUTECT( vardict_samplename_ch,vcf_std_for_bcftools_concat,vcf_mutect_for_bcftools_concat,vardict_filtered_index_isolated_standard,mutect_filtered_index_isolated )
        





    } else {
        println "Running the SNPs/indels workflow in single sample mode."
        CALL_VARIANTS_SINGLE (params.input,params.bed,params.fasta,params.fai)

        vardict_filtered_vcfs = CALL_VARIANTS_SINGLE.out.vardict_filtered_vcf

        vardict_filtered_vcfs
            .map{ standard_vcf,complexvar_vcf -> standard_vcf}
            .set{ vardict_filtered_vcf_standard }

        vardict_filtered_vcfs
            .map{ standard_vcf,complexvar_vcf -> complexvar_vcf}
            .set{ vardict_filtered_vcf_complexvar }

        ref_fasta = CALL_VARIANTS_SINGLE.out.genome_fasta_file
        ref_fasta_index = CALL_VARIANTS_SINGLE.out.genome_fasta_index_file

        BCFTOOLS_ALL_VARDICT( vardict_filtered_vcf_standard,ref_fasta,ref_fasta_index )

        meta_plus_vardict_norm_and_sorted_vcf_standard = BCFTOOLS_ALL_VARDICT.out.meta_plus_vardict_norm_and_sorted_vcf
        vardict_filtered_index_isolated_standard = BCFTOOLS_ALL_VARDICT.out.vardict_filtered_index_final


        BCFTOOLS_ALL_VARDICT_COMPLEXVAR( vardict_filtered_vcf_complexvar,ref_fasta,ref_fasta_index )

        meta_plus_vardict_norm_and_sorted_vcf_complexvar = BCFTOOLS_ALL_VARDICT_COMPLEXVAR.out.meta_plus_vardict_norm_and_sorted_vcf
        vardict_filtered_index_isolated_complexvar = BCFTOOLS_ALL_VARDICT_COMPLEXVAR.out.vardict_filtered_index_final


        meta_plus_vardict_norm_and_sorted_vcf_standard
            .map{ sample,vcf -> sample}
            .set{ vardict_samplename_ch }

        meta_plus_vardict_norm_and_sorted_vcf_standard
            .map{ sample,vcf -> vcf}
            .set{ vardict_norm_and_sorted_vcf_standard }

        meta_plus_vardict_norm_and_sorted_vcf_complexvar
            .map{ sample,vcf -> vcf}
            .set{ vardict_norm_and_sorted_vcf_complexvar }

        vardict_norm_and_sorted_vcf_standard
            .map{ create_std_filename(it) }
            .set{ vcf_std_for_bcftools_concat }

        vardict_norm_and_sorted_vcf_complexvar
            .map{ create_complexvar_filename(it) }
            .set{ vcf_complexvar_for_bcftools_concat }

        
        BCFTOOLS_CONCAT_VARDICTS( vardict_samplename_ch,vardict_regular_norm_and_sorted_vcf,vardict_complexvar_norm_and_sorted_vcf,vardict_index,vardict_complexvar_index )
        
        vardict_concatenated_vcf = BCFTOOLS_CONCAT_VARDICTS.out.sample_plus_vardict_concat_vcf

    }



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
