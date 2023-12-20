//
// Check input samplesheet and get read channels
//

include { VARDICTJAVA } from '../../modules/nf-core/vardictjava/main'

workflow CALL_VARIANTS {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    input_sample_sheet = Channel.fromPath(samplesheet)
        .splitCsv ( header:true, sep:',' )
        .set{ reads }

    reads
        .map { create_meta_channel(it) }
        .view()
        .set { meta }

    reads
        .map { create_meta2_channel(it) }
        .view()
        .set { meta2 }

    reads
        .map { create_meta3_channel(it) }
        .view()
        .set { meta3 }

    
    VARDICTJAVA(meta,meta2,meta3)

    vcf = VARDICTJAVA.out.vcf
    
    emit:
    vcf
    //fasta                                     // channel: [ val(meta), [ reads ] ]
     // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_meta_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample_name
    bams = [ meta,row.bam1,'/Users/naidur/ACCESS/access_pipeline/test_data/placeholder.txt', row.bed   ]
}

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
