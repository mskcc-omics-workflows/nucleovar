//
// Check input samplesheet and get read channels
//

include { VARDICTJAVA } from '../../modules/nf-core/vardictjava/main'
include { VARDICT_FILTER } from '../../modules/local/vardict_filter'

workflow CALL_VARIANTS {
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

    inputs 
        .map{ row -> [row.case_sample_name,row.control_sample_name,row.case_bam,row.control_bam,row.case_bai,row.control_bai]
        }
        .set{ mutect1_input_set1 }
    
    standard_bed_file
        .combine(genome_fasta_file)
        .combine(genome_fasta_index_file)
        .set{ mutect1_input_set2 }



    
VARDICTJAVA(vardict_input_set1,vardict_input_set2,vardict_input_set3)
vardict_vcf = VARDICTJAVA.out.vcf
vardict_vcf.combine(bam_ch).set{ for_vardictfilter_ch }



// MUTECT1(mutect1_input_set1,mutect1_input_set2)
// mutect_vcf = MUTECT1.out.vcf



inputs
    .map { create_vardictfilter_names_channel(it)  }
    .set{ vardictfilter_samplenames_ch }


VARDICT_FILTER(for_vardictfilter_ch,vardictfilter_samplenames_ch)
//vardict_filtered_vcf = VARDICT_FILTER.out.filtered_vcf


//emit:
//vardict_sample_and_vcf
//mutect_vcf
//vardict_filtered_vcf


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



