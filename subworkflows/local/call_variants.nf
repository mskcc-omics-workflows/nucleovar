//
// Check input samplesheet and get read channels
//

include { VARDICTJAVA } from '../../modules/nf-core/vardictjava/main'

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


    VARDICTJAVA(vardict_input_set1,vardict_input_set2,vardict_input_set3)
    vardict_vcf = VARDICTJAVA.out.vcf

    // // // // MUTECT(mutectinputs,bed_fasta_fai)
    // // // // mutect_vcf = MUTECT.out.vcf
    // // // //
    emit:
    vardict_vcf


}


// def create_bams_and_bais_channel(LinkedHashMap row) {
//     // create meta map
//     def meta = [:]

//     // set bams variable
//     if (row.control_bam.isEmpty() || row.case_bam.isEmpty()) {
        
//         if (row.control_bam.isEmpty()) {
//             bams = row.case_bam
//         } else if (row.case_bam.isEmpty()) {
//             bams = row.control_bam
//         }
//     } else {
//         if (!row.control_bam.isEmpty() && !row.case_bam.isEmpty()) {
//             bams = [row.control_bam,row.case_bam]

//         }
//     }
//     // set bais variable
//     if (row.control_bai.isEmpty() || row.case_bai.isEmpty()) {
        
//         if (row.control_bai.isEmpty()) {
//             bais = [row.case_bai]
//         } else if (row.case_bai.isEmpty()) {
//             bais = [row.control_bai]
//         }
//     } else {
//         if (!row.control_bai.isEmpty() && !row.case_bai.isEmpty()) {
//             bais = [row.control_bai, row.case_bai]
//         }
//     }

//     if (row.control_sample_name && row.case_sample_name) {
//         // Both sample name columns are non-empty
//         meta.id = "${row.control_sample_name}|${row.case_sample_name}"
//         //inputs = [meta,[bams],bais]
//         bams = [bams]
//     } else if (row.control_sample_name) {
//         // Only control_sample_name is non-empty
//         meta.id = "${row.control_sample_name}"
//         inputs = [meta,[row.control_bam],[row.control_bai]]
//     } else if (row.case_sample_name) {
//         // Only case_sample_name is non-empty
//         meta.id = "${row.case_sample_name}"
//         inputs = [meta,[row.case_bam],[row.case_bai]]
//     } else {
//         // Both are empty, handle this scenario as needed
//        error "Sample Name columns are both empty. Please re-check your input samplesheet."
//     }
// }

def create_samplenames_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    if (row.control_sample_name && row.case_sample_name) {
        // Both sample name columns are non-empty
        meta.id = "${row.control_sample_name}|${row.case_sample_name}"
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
            bams = [row.control_bam,row.case_bam]
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
            bais = [row.control_bai, row.case_bai]
            baisinput = [bais]
        }
    }
}



