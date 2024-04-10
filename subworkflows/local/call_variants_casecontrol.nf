//
// Check input samplesheet and get read channels
//

include { VARDICTJAVA } from '../../modules/nf-core/vardictjava/main'
include { VARDICT_FILTER } from '../../modules/local/vardict_filter'
//include { MUTECT1 } from '../../modules/modules/modules/msk/mutect1'
include { MUTECT_FILTER } from '../../modules/local/mutect_filter'

workflow CALL_VARIANTS_CASECONTROL {
    take:
    samplesheet
    bed
    fasta
    fai

    main:

    // taking in input samplesheet and splitting to create channels for vardictjava and mutect 
    standard_bed_file = Channel.fromPath(bed)
    genome_fasta_file = Channel.fromPath(fasta)
    genome_fasta_index_file = Channel.fromPath(fai)

    input_sample_sheet = Channel.fromPath(samplesheet)
        .splitCsv ( header:true, sep:',' )
        .set{ inputs }


    inputs
        .map { create_samplenames_channel(it)  }
        .set{ name_ch }

    inputs
        .map { create_bams_channel(it)  }
        .set{ bam_ch }

    inputs
        .map { create_bais_channel(it)  }
        .set{ bai_ch }

    name_ch 
        .combine(bam_ch)
        .combine(bai_ch)
        .combine(standard_bed_file)
        .set { vardict_input_set1 }

    name_ch 
        .combine(genome_fasta_file)
        .set { vardict_input_set2 }

    name_ch 
        .combine(genome_fasta_index_file)
        .set { vardict_input_set3 }



    
VARDICTJAVA(vardict_input_set1,vardict_input_set2,vardict_input_set3)
vardict_vcf = VARDICTJAVA.out.vcf
vardict_vcf.combine(bam_ch).set{ for_vardictfilter_ch }


inputs
    .map { create_vardictfilter_names_channel(it)  }
    .set{ vardictfilter_samplenames_ch }


VARDICT_FILTER(for_vardictfilter_ch,vardictfilter_samplenames_ch)
vardict_filtered_vcf = VARDICT_FILTER.out.filtered_vardict_vcf
complex_variants_vardict_filtered_vcf = VARDICT_FILTER.out.complex_variants_vardict_vcf
std_vardict_filter_output_txt = VARDICT_FILTER.out.std_vardict_filter_output


inputs
    .map { create_msk_mutect1_inputs_channel(it)  }
    .set{ mutect1_input_set1_ch }

genome_fasta_dict_file = Channel.fromPath("/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/Homo_sapiens_assembly19.dict")

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


mutect_filtered_vcf = Channel.fromPath("/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/DONOR22-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex-C-2HXC96-P001-d01_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.mutect_filter.mutect.vcf")


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
    
    // set bais variable
    if (row.control_bai.isEmpty() || row.case_bai.isEmpty()) {
        
        if (row.control_bai.isEmpty()) {
            bais = [row.case_bai]
        } else if (row.case_bai.isEmpty()) {
            bais = [row.control_bai]
        }
    } else {
        if (!row.control_bai.isEmpty() && !row.case_bai.isEmpty()) {
            bais = [file(row.control_bai), file(row.case_bai)]
            baisinput = [bais]
        }
    }
}

def create_vardictfilter_names_channel(LinkedHashMap row) {
    // create meta map

    // set bams variable
    if (row.control_bam.isEmpty() || row.case_bam.isEmpty()) {
        
        if (row.control_bam.isEmpty()) {
            bams = file(row.case_bam).getBaseName()
        } else if (row.case_bam.isEmpty()) {
            bams = file(row.control_bam).getBaseName()
        }
    } else {
        if (!row.control_bam.isEmpty() && !row.case_bam.isEmpty()) {
            control_bam = file(row.control_bam).getName()
            case_bam = file(row.case_bam).getBaseName()
            samplename = "${control_bam}|${case_bam}"
            baminput = samplename

        }
    }
}


def create_msk_mutect1_inputs_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    if (row.control_sample_name && row.case_sample_name) {
        // Both sample name columns are non-empty
        meta.case_id = row.case_sample_name
        meta.control_id = row.control_sample_name
        inputs = [meta,row.case_bam,row.control_bam,row.case_bai,row.control_bai]
    } else if (row.control_sample_name) {
        // Only control_sample_name is non-empty
        error "Case sample ID is missing. MuTect1 requires both the case and control sample ID. Please confirm if you are running the case-control workflow."
    } else if (row.case_sample_name) {
        // Only case_sample_name is non-empty
        error "Control sample ID is missing. MuTect1 requires both the case and control sample ID. Please confirm if you are running the case-control workflow."
    } else {
        // Both are empty, handle this scenario as needed
       error "Sample Name columns are both empty. Please re-check your input samplesheet."
    }
}



