

//import bcftools subworkflow for vardict filtered regular
//include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { VCF2MAF } from '../../modules/nf-core/vcf2maf/main'
include { TRACEBACK } from '../../subworkflows/msk/traceback/main'




workflow MODULE4 {
    take:
    input_maf
    case_bam_and_index
    fasta
    fasta_fai  
    

    main:
    //run vcf2maf perl module (placeholder at the moment while temporarily testing out traceback)
    //VCF2MAF(vcf,reference_fasta,reference_fasta_index)

    emptyLists = Channel.from([], [], [], [])
    def bamnames = [patient:'test',id:'sample']
    names = Channel.create(bamnames)
    
    //bams = Channel.from(bamnames).merge(case_bam).merge(case_bai).concat(emptyLists)
    case_bam_and_index
        .map{ it -> tuple([patient:'test',id:'sample'],file(it[0]),file(it[1]),[],[],[],[])}
        .set{ bams }
    
    
    mafs = Channel.from([patient:'test',id:"sample"]).merge(input_maf)
    //header = Channel.from([initial: file('/Users/naidur/Desktop/header.txt'), genotype: file('/Users/naidur/Desktop/genotype.txt')])
    headerfile = Channel.from(file('/Users/naidur/Desktop/header.txt'))
    genotypefile = Channel.from(file('/Users/naidur/Desktop/genotype.txt'))
    headerfile.combine(genotypefile).map{ [initial:headerfile,genotype:genotypefile]}.set{ header }
    
    //bams.map { tuple( 'patient': it[0]['patient'], *it ) }.view()
    

    TRACEBACK( bams,mafs,[initial:file('/Users/naidur/Desktop/maf_header.txt'),genotype:file('/Users/naidur/Desktop/maf_header_genotype.txt')],fasta,fasta_fai )


    
    
}


def mapMeta(metaMap) {
    // create meta map
    def meta = [:]
    meta.initial = '/Users/naidur/Desktop/header.txt'
    inputs = meta 
}

class MyChannel extends Channel {
    def getAt(String key) {
        delegate.header[key]
    }
}


