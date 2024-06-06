//
// Subworkflow with functionality specific to the msk/nucleovar pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    aux_bams          // string: Filepath to file containing curated simplex/duplex bams for traceback

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    Channel
        .fromPath(input)
        .splitCsv(header: true)
        .set { ch_samplesheet }

    ch_samplesheet
        .map {row -> row.sample_id }
        .collect()
        .map{ [id:"${it[1]}_${it[0]}",case_id:it[0],control_id:it[1]]}
        .set{ sample_id_names_ch }
    
    ch_samplesheet
        .branch{ row -> standard: row.type == "standard"
            return tuple([patient_id: row.patient_id,sample_id: row.sample_id],row.standard_bam,row.standard_bai) }
        .set{ standard_bams_ch }
    
    ch_samplesheet
        .branch{ row -> tumor: row.type == "case"
            return tuple(row.duplex_bam,row.duplex_bai) }
        .set{ case_bams_ch }

    ch_samplesheet
        .branch{ row -> tumor: row.type == "control"
            return tuple(row.duplex_bam,row.duplex_bai) }
        .set{ control_bams_ch }

    sample_id_names_ch
        .combine(control_bams_ch)
        .combine(case_bams_ch)
        .set{ duplex_bams_ch }

    ch_samplesheet
        .branch{ row -> tumor: row.type == "case"
            return tuple([patient: row.patient_id,id: row.sample_id],[],[],file(row.duplex_bam),file(row.duplex_bai),file(row.simplex_bam),file(row.simplex_bai)) }
        .set{ case_bams_for_traceback_ch }

    


    // donor38_duplex_bam = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_duplex_bams_dmp/versions/v2.0/DONOR38-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam")
    // donor38_duplex_bai = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_duplex_bams_dmp/versions/v2.0/DONOR38-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bai")

    // donor38_simplex_bam = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_simplex_bams_dmp/versions/v2.0/DONOR38-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam")
    // donor38_simplex_bai = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_simplex_bams_dmp/versions/v2.0/DONOR38-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bai")

    // aux_bams_ch = Channel.from(tuple([patient:'null',id:'DONOR38-TP'],
    //                 [],
    //                 [],
    //                 donor38_duplex_bam,
    //                 donor38_duplex_bai,
    //                 donor38_simplex_bam,
    //                 donor38_simplex_bai))

    Channel
        .fromPath(aux_bams)
        .splitCsv(header: true)
        .map{ row -> tuple(row.sample_id,row.simplex_path,row.duplex_path)}
        .set{ temp_ch }

    temp_ch
        .map{ sample_id,simplex_path,duplex_path -> tuple([patient:'null',id:sample_id],
                    [],
                    [],
                    duplex_path.toString(),
                    (file(duplex_path).parent/file(duplex_path).baseName +'.bai').toString(),
                    simplex_path.toString(),
                    (file(simplex_path).parent/file(simplex_path).baseName +'.bai').toString()) }
        .set{ aux_bams_ch }

    

    

    emit:
    samplesheet = ch_samplesheet
    sample_id_names = sample_id_names_ch
    standard_bams = standard_bams_ch
    case_bams = case_bams_ch
    control_bams = control_bams_ch
    duplex_bams = duplex_bams_ch
    versions    = ch_versions
    case_bams_for_traceback = case_bams_for_traceback_ch
    aux_bams = aux_bams_ch
    
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}