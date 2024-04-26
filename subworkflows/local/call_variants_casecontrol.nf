//
// Check input samplesheet and get read channels
//

include { VARDICTJAVA } from '../../modules/nf-core/vardictjava/main'
include { VARDICT_FILTER } from '../../modules/local/vardict_filter'

workflow CALL_VARIANTS_CASECONTROL {
    take:
    samplesheet
    fasta
    fai
    dict
    bed

    main:

    // taking in input samplesheet and splitting to create channels for vardictjava and mutect 
    

    Channel
            .fromPath(params.input)
            .splitCsv(header: true)
            .set{ input_samplesheet }

        // setting up sample names channel
        input_samplesheet
            .map {row -> row.sample_id }
            .collect()
            .map{ [id:"${it[1]}_${it[0]}",case_id:it[0],control_id:it[1]]}
            .set{ sample_id_names_ch }

        //setting up standard bams channel
        input_samplesheet
            .branch{ row -> standard: row.type == "standard"
                    return tuple([patient_id: row.patient_id,sample_id: row.sample_id],row.standard_bam,row.standard_bai) }
            .set{ standard_bams_ch }
        
        // setting up duplex bams channel
        input_samplesheet
            .branch{ row -> tumor: row.type == "case"
                            return tuple(row.duplex_bam,row.duplex_bai) }
            .set{ case_bams_ch }
        
        input_samplesheet
            .branch{ row -> tumor: row.type == "control"
                            return tuple(row.duplex_bam,row.duplex_bai) }
            .set{ control_bams_ch }
        
        // setting up inputs for vardict
        sample_id_names_ch
            .combine(control_bams_ch)
            .combine(case_bams_ch)
            .set{ duplex_bams }
        
        duplex_bams
            .map{ meta,control_bam,control_bai,case_bam,case_bai -> tuple(meta,[file(control_bam),file(case_bam)],[file(control_bai),file(case_bai)])}
            .combine(bed)
            .set{ vardict_input_set1 }

        sample_id_names_ch
            .combine(fasta)
            .set{ vardict_input_set2 }

        sample_id_names_ch
            .combine(fai)
            .set{ vardict_input_set3 }



    
    VARDICTJAVA(vardict_input_set1,vardict_input_set2,vardict_input_set3)
    vardict_vcf = VARDICTJAVA.out.vcf
    vardict_vcf.map{ meta,vcf -> file(vcf)}.set{ vardict_vcf_isolated }

    duplex_bams
        .map{ meta,control_bam,control_bai,case_bam,case_bai -> tuple([file(control_bam),file(case_bam)])}
        .set{ bams_for_vardict_filter }

    VARDICT_FILTER( sample_id_names_ch,vardict_vcf_isolated,bams_for_vardict_filter )

    vardict_filtered_vcfs = VARDICT_FILTER.out.filtered_vardict_vcf

    vardict_filtered_vcfs.map{ standard,complex_var -> standard}.set{ standard_vcf }
    vardict_filtered_vcfs.map{ standard,complex_var -> standard}.set{ complexvar_vcf }
    vardict_filter_output_txt = VARDICT_FILTER.out.std_vardict_filter_output

<<<<<<< HEAD
=======

inputs
    .map { create_msk_mutect1_inputs_channel(it)  }
    .set{ mutect1_input_set1_ch }

genome_fasta_dict_file = Channel.fromPath("/juno/cmo/access/production/resources/reference/current/Homo_sapiens_assembly19.dict")

name_ch
    .combine( standard_bed_file )
    .combine( genome_fasta_file )
    .combine( genome_fasta_index_file )
    .combine( genome_fasta_dict_file )
    .set{ mutect1_input_set2_ch }

//MUTECT1( mutect1_input_set1_ch,mutect1_input_set2_ch )
 

//MUTECT_FILTER(for_vardictfilter_ch)
//mutect_filtered_vcf = MUTECT_FILTER.out.filtered_mutect_vcf
//std_vardict_filter_output_txt = MUTECT_FILTER.out.std_mutect_filter_output


mutect_filtered_vcf = Channel.fromPath("/home/buehlere/snp_indels_conversion/test_nucleo_var/DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-C-2HXC96-P001-d01_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.mutect_filter.mutect.vcf")


emit:
vardict_filtered_vcf
complex_variants_vardict_filtered_vcf
mutect_filtered_vcf
genome_fasta_file
genome_fasta_index_file


}




def create_samplenames_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    if (row.control_sample_name && row.case_sample_name) {
        // Both sample name columns are non-empty
        meta.id = "${row.control_sample_name}-${row.case_sample_name}"
        inputs = [meta]
    } else if (row.control_sample_name) {
        // Only control_sample_name is non-empty
        meta.id = "${row.control_sample_name}"
        inputs = [meta]
    } else if (row.case_sample_name) {
        // Only case_sample_name is non-empty
        meta.id = "${row.case_sample_name}"
        inputs = [meta]
    } else {
        // Both are empty, handle this scenario as needed
       error "Sample Name columns are both empty. Please re-check your input samplesheet."
    }
}

def create_bams_channel(LinkedHashMap row) {
    // create meta map

    // set bams variable
    if (row.control_bam.isEmpty() || row.case_bam.isEmpty()) {
        
        if (row.control_bam.isEmpty()) {
            bams = row.case_bam
        } else if (row.case_bam.isEmpty()) {
            bams = row.control_bam
        }
    } else {
        if (!row.control_bam.isEmpty() && !row.case_bam.isEmpty()) {
            bams = [file(row.control_bam),file(row.case_bam)]
            baminput = [bams]

        }
    }
}

def create_bais_channel(LinkedHashMap row) {
>>>>>>> feature/call_variants_lsf
    

    
    emit:
    standard_vcf
    complexvar_vcf
    vardict_filter_output_txt
    standard_bams_ch
    duplex_bams
    sample_id_names_ch


}





