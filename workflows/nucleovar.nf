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
include { CALL_VARIANTS_CASECONTROL     } from '../subworkflows/local/call_variants_casecontrol'
include { BCFTOOLS_CONCAT_VARDICTS     } from '../subworkflows/local/bcftools_concat_vardicts'
include { MODULE4     } from '../subworkflows/local/module4'
include { GUNZIP_FILES     } from '../modules/local/gunzip_files'
include { MUTECT1        } from '../modules/msk/mutect1'
include { MUTECT_FILTER     } from '../modules/local/mutect_filter'
include { BCFTOOLS_CONCAT_WITH_MUTECT     } from '../subworkflows/local/bcftools_concat_with_mutect'
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
    ch_samplesheet // channel: samplesheet read in from --input
    // TODO: remove hard coded paths
    // TODO: samplesheet needs to use new template format. 
    // Sample sheet should be process in nucleo var 
    // aka it' the controller for the samplesheet 
    // TODO: we should be using regular MSK-ACCESS images, not internals
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    // sanity check to see if tumor and normal samples are included-- kicks off case control method of running pipeline
    println params.input
    def sampleSheet = file(params.input).readLines().collect { it.split(",") }
    def allColumnsHaveValue = sampleSheet.every { row ->
    row.every { cell -> cell.trim() }
}

    if (allColumnsHaveValue) {
        println "Running the SNPs/indels workflow in case-control mode."


        bed = Channel.from(params.bed)
        fasta_ref = Channel.from(params.fasta)
        fasta_index = Channel.from(params.fai)
        fasta_dict = Channel.from(params.dict)



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

        BCFTOOLS_VARDICT( vardict_filtered_vcf_complexvar,vardict_filtered_vcf_standard,ref_fasta,ref_fasta_index )

        CALL_VARIANTS_CASECONTROL (params.input,fasta_ref,fasta_index,fasta_dict,bed)
        vardict_filtered_vcf_standard = CALL_VARIANTS_CASECONTROL.out.standard_vcf
        vardict_filtered_vcf_complexvar = CALL_VARIANTS_CASECONTROL.out.complexvar_vcf
        duplex_bams = CALL_VARIANTS_CASECONTROL.out.duplex_bams
        sample_ids = CALL_VARIANTS_CASECONTROL.out.sample_id_names_ch
        
        

        sample_ids.combine(vardict_filtered_vcf_standard).set{ standard_vcf_for_bcftools }
        sample_ids.combine(vardict_filtered_vcf_complexvar).set{ complexvar_vcf_for_bcftools }
        BCFTOOLS_VARDICT( vardict_filtered_vcf_complexvar,vardict_filtered_vcf_standard,fasta_ref,fasta_index )
    
        vardict_concat_vcf = BCFTOOLS_VARDICT.out.vardict_concat_vcf
        vardict_index = BCFTOOLS_VARDICT.out.vardict_index
        
        // // // // MUTECT1 MODULE
        // // // temporary code for putting together inputs for mutect1 module (will be deprececated when moving to new samplesheet)
        // Channel
        //     .fromPath(params.input)
        //     .splitCsv(header: true)
        //     .set{ input_samplesheet }

        // input_samplesheet
        //     .map{ row -> tuple(row.case_bam, row.control_bam, row.case_bai, row.control_bai) }
        //     .set{ bams_ch }

        // input_samplesheet
        //     .map{ row -> [case_id:row.case_sample_name,control_id:row.control_sample_name,id:"${row.case_sample_name}_${row.control_sample_name}"]}
        //     .set{ sample_id_names_ch }

        duplex_bams.map{ meta,control_bam,control_bai,case_bam,case_bai -> tuple(case_bam,control_bam,case_bai,control_bai)}.set{ bams_for_mutect }
        
        sample_ids
            .combine(bams_for_mutect)
            .set{ input1_for_mutect }

        bed
            .combine(fasta_ref)
            .combine(fasta_index)
            .combine(fasta_dict)
            .set{ input2_for_mutect }



        MUTECT1(input1_for_mutect,input2_for_mutect)
        mutect_vcf = MUTECT1.out.mutect_vcf
        mutect_txt = MUTECT1.out.standard_mutect_output

        mutect_txt.map{ meta,file -> file}.set{ mutect_txt_isolated }

        sample_ids.combine(mutect_vcf).combine(mutect_txt_isolated).map{ meta1,meta2,vcf,txt -> tuple(meta1,vcf,txt)}.set{ input1_for_mutect_filter }
        

        //MUTECT_FILTER(input1_for_mutect_filter,fasta_ref)

        
        // temp testing mutect filtered vcf (permission error in mutect filter)
        mutect_filtered_vcf = Channel.fromPath("/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-C-2HXC96-P001-d01_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.mutect_filter.mutect.vcf")
    
        BCFTOOLS_MUTECT( mutect_filtered_vcf,fasta_ref,fasta_index )
        mutect_vcf = BCFTOOLS_MUTECT.out.standard_norm_sorted_vcf
        mutect_index = BCFTOOLS_MUTECT.out.mutect_index


        vardict_concat_vcf.map{ id,vcf -> vcf}.set{ vardict_concat_vcf_isolated }
        BCFTOOLS_CONCAT_WITH_MUTECT( sample_ids,vardict_concat_vcf_isolated,mutect_vcf,vardict_index,mutect_index )
        sample_plus_final_concat_vcf = BCFTOOLS_CONCAT_WITH_MUTECT.out.sample_plus_final_concat_vcf

        //BCFTOOLS_ANNOTATE(vardict_concat_vcf,mutect_concat_vcf)

        //annotated_vcf = BCFTOOLS_ANNOTATE.out.vcf

        // testing inputs for traceback temporarily


        // // code to prepare simplex inputs as a channel

        // // code to prepare duplex inputs as a channel

        //rules_json = Channel.fromPath(params.rules_json)
        //MODULE4( annotated_vcf,bams_ch,fasta_ref,fasta_index,rules_json  )

    }
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    emit:
    //multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
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
