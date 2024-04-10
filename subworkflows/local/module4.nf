

//import bcftools subworkflow for vardict filtered regular
//include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { GUNZIP } from '../../modules/nf-core/gunzip/main' 
include { VCF2MAF } from '../../modules/nf-core/vcf2maf/main'
include { TRACEBACK } from '../../subworkflows/msk/traceback/main'




workflow MODULE4 {
    take:
    input_maf
    case_bam
    case_bai
    fasta
    fasta_fai  
    

    main:
    //def meta = [:]
    //meta.patient = 'null'
    
    
    //bams --> tumor: [[patient:null, id:'sample'], standard.bam, standard.bam.bai, [], [], [], []], normal: [[patient:null, id:'sample'], standard.bam, standard.bam.bai, [], [], [], []]
    //bams =  [[patient:null, id:'sample'], case_bam, case_bai, [], [], [], []]

    emptyLists = Channel.from([], [], [], [])
    def bamnames = [patient:'test',id:'sample']
    names = Channel.create(bamnames)
    
    bams = Channel.from(bamnames).merge(case_bam).merge(case_bai).merge(emptyLists)
    
    
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


