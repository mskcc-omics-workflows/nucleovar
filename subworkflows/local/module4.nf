

//import bcftools subworkflow for vardict filtered regular
//include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { VCF2MAF } from '../../modules/nf-core/vcf2maf/main'
include { TRACEBACK } from '../../subworkflows/msk/traceback/main'




workflow MODULE4 {
    take:
    vcf
    bams
    fasta
    fasta_fai 
    rules_json 
    

    main:
    // temporary input maf for testing purposes 
    input_maf = Channel.fromPath("/Users/naidur/ACCESS/access_pipeline/dev/postprocessing/test_data/tagged_by_hotspots.maf")

    //run vcf2maf perl module (placeholder at the moment while temporarily testing out traceback)
    vcf.combine(fasta).set{ input_for_vcf2maf }
    VCF2MAF( input_for_vcf2maf )
    maf = VCF2MAF.out.maf

    TAG_BY_ACCESS( maf,rules_json )





    emptyLists = Channel.from([], [], [], [])
    def bamnames = [patient:'test',id:'sample']
    names = Channel.create(bamnames)

    
    // standard channel input
    bams
        .map{ it -> tuple([patient:'test',id:'sample'],file(it[0]),file(it[1]),[],[],[],[])}
        .set{ bams }
    
    // simplex/duplex channel input 
    mafs = Channel.from([patient:'test',id:"sample"]).merge(tagged_maf)
    
    headerfile = Channel.from(file('/Users/naidur/Desktop/header.txt'))
    genotypefile = Channel.from(file('/Users/naidur/Desktop/genotype.txt'))
    headerfile.combine(genotypefile).map{ [initial:headerfile,genotype:genotypefile]}.set{ header }
    
    
    

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


