//
// Check input samplesheet and get read channels
//

//include { VARDICTJAVA } from '../../modules/nf-core/vardictjava/main'

workflow CALL_VARIANTS {
    take:
    input_dir
    samplesheet

    main:
    
    // taking in input directory where files are located and creating input channels for vardictjava
    def allFiles = new File(input_dir).listFiles()
    def bam_exists = allFiles.any { it.name.endsWith('.bam') }
    def bai_exists = allFiles.any { it.name.endsWith('.bai') }  
    def bed_exists = allFiles.any { it.name.endsWith('.bed') }  
    def fasta_exists = allFiles.any { it.name.endsWith('.fasta') }  
    def fai_exists = allFiles.any { it.name.endsWith('.fai') } 

    if (bam_exists && bai_exists && bed_exists && fasta_exists && fai_exists) {

        Channel
            .fromFilePairs("${input_dir}/*.{bam,bai}", size: 2)
            .set { bam_bai_files }
        
        Channel
            .fromPath("${input_dir}/*.bed")
            .set{ bed }
        Channel
            .fromPath("${input_dir}/*.fasta")
            .set{ fasta }
        Channel
            .fromPath("${input_dir}/*.fasta.fai")
            .set{ fai }

        sample_names = bam_bai_files.map { tuple -> tuple[0] }.collect()
        bams = bam_bai_files.map { tuple -> tuple[1][1] }.collect()
        bais = bam_bai_files.map { tuple -> tuple[1][0] }.collect()

    
    }
    else {
        error "Error: Mandatory Input files not found in ${input_dir}. Please make sure you are providing sample name(s), bam(s), bai(s), fasta, and BED file."
    }


    meta = sample_names.map { [it] }
    bamslist = bams.map { [it] }
    baislist = bais.map { [it] }

    meta1 = meta.combine(bamslist).combine(baislist).combine(bed)
    meta2 = meta.combine(fasta)
    meta3 = meta.combine(fai)

    // taking in input samplesheet and splitting to create channels for mutect 
    input_sample_sheet = Channel.fromPath(samplesheet)
        .splitCsv ( header:true, sep:',' )
        .set{ reads }

    reads
        .map { create_mutect_channel(it) }
        .view()
        .set { mutectinputs }


    // VARDICTJAVA(meta1,meta2,meta3)
    // vcf = VARDICTJAVA.out.vcf

    // MUTECT(mutectinputs,bed_fasta_fai)
    // mutect_vcf = MUTECT.out.vcf



}

def create_mutect_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    def bams = [:]
    def bais = [:]
    def bed_fasta_fai = [:]
    

    //samplenames
    if (row.control_sample_name.isEmpty() || row.case_sample_name.isEmpty()) {
        // Check if either column1 or column2 is empty
        if (row.control_sample_name.isEmpty()) {
            meta = [row.case_sample_name]
        } else if (row.case_sample_name.isEmpty()) {
            meta = [row.control_sample_name]
        }
    } else {
        if (!row.control_sample_name.isEmpty() && !row.case_sample_name.isEmpty()) {
            meta = [row.control_sample_name,row.case_sample_name]
        }
    }
    //bams
    if (row.control_bam.isEmpty() || row.case_bam.isEmpty()) {
        
        if (row.control_bam.isEmpty()) {
            bams = [row.case_bam]
        } else if (row.case_bam.isEmpty()) {
            bams = [row.control_bam]
        }
    } else {
        if (!row.control_bam.isEmpty() && !row.case_bam.isEmpty()) {
            bams = [row.control_bam,row.case_bam]
        }
    }
    // bais
    if (row.control_bai.isEmpty() || row.case_bai.isEmpty()) {
        
        if (row.control_bai.isEmpty()) {
            bais = [row.case_bai]
        } else if (row.case_bai.isEmpty()) {
            bais = [row.control_bai]
        }
    } else {
        if (!row.control_bai.isEmpty() && !row.case_bai.isEmpty()) {
            bais = [row.control_bai,row.case_bai]
        }
    }

    mutectinputs = [meta,bams,bais]
    bed_fasta_fai = [row.bed,row.fasta,row.fai]
    
}

    
    