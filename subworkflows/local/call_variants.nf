//
// Check input samplesheet and get read channels
//

include { VARDICTJAVA } from '../../modules/nf-core/vardictjava/main'

workflow CALL_VARIANTS {
    take:
    sample_names
    bams 
    bais
    bed
    fasta
    fai 

    main:
    meta = sample_names.map { [it] }
    bamslist = bams.map { [it] }
    baislist = bais.map { [it] }

    meta1 = meta.combine(bamslist).combine(baislist).combine(bed)
    meta2 = meta.combine(fasta)
    meta3 = meta.combine(fai)

    
    

    VARDICTJAVA(meta1,meta2,meta3)

    vcf = VARDICTJAVA.out.vcf

    emit:
    vcf


    
}


// // Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_meta_and_bams_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = sample_names
    result = [meta.id,bams]
}


/*
def create_meta2_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample_name
    fasta = [ meta, row.fasta ]
}

def create_meta3_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample_name
    fasta = [ meta, '/Users/naidur/ACCESS/access_pipeline/test_data/placeholder2.txt']
}
*/
