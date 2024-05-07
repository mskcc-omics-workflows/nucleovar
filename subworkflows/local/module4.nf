

//import bcftools subworkflow for vardict filtered regular
//include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { VCF2MAF } from '../../modules/nf-core/vcf2maf/main'
include { TRACEBACK } from '../../subworkflows/msk/traceback/main'




workflow MODULE4 {
    take:
    //vcf
    duplex_bams
    fasta
    fasta_fai 
    //rules_json 
    

    main:
    // temporary input maf for testing purposes 
    input_maf = Channel.fromPath("/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/C-2HXC96-P001-d01.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots.maf")
    simplex_bam1 = "/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/DONOR10-T_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam"
    simplex_bai1 = "/Users/naidur/ACCESS/access_pipeline/test_data/test_data/MSK_data/DONOR10-T_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bai"

    //simplex_bam1.combine(simplex_bai1).set{test_simplex_bams}

    //run vcf2maf perl module (placeholder at the moment while temporarily testing out traceback)
    // vcf.combine(fasta).set{ input_for_vcf2maf }
    // VCF2MAF( input_for_vcf2maf )
    // maf = VCF2MAF.out.maf

    // TAG_BY_ACCESS( maf,rules_json )


    emptyLists = Channel.from([], [], [], [])
    def bamnames = [patient:'test',id:'sample']
    names = Channel.create(bamnames)

    
    // standard channel input
    // duplex_bams
    //     .map{ meta,bam1,bai1,bam2,bai2 -> tuple([patient:'test',id:'C-2HXC96-P001-d01.DONOR22-TP.combined-variants'],bam1,bai1,[],[],[],[])}
    //     .view()
    //     .set{ duplex_bams_for_traceback }

    
    duplex_bams
        .map{ meta,bam1,bai1,bam2,bai2 -> tuple([patient:'test',id:'C-2HXC96-P001-d01.DONOR22-TP.combined-variants'],[],[],bam1,bai1,simplex_bam1,simplex_bai1)}
        .set{ duplex_simplex_bams_for_traceback }

    // code to extract the aux bams
    //map to match structure above and combine

    
    
    // simplex/duplex channel input 
    mafs = Channel.from([patient:'test',id:"C-2HXC96-P001-d01.DONOR22-TP.combined-variants"]).merge(input_maf)
    
    headerfile = Channel.from(file('/Users/naidur/Desktop/header.txt'))
    genotypefile = Channel.from(file('/Users/naidur/Desktop/genotype.txt'))
    headerfile.combine(genotypefile).map{ [initial:headerfile,genotype:genotypefile]}.set{ header }


    
    
    

    TRACEBACK( duplex_simplex_bams_for_traceback,mafs,[initial:file('/Users/naidur/Desktop/maf_header.txt'),genotype:file('/Users/naidur/Desktop/maf_header_genotype.txt')],fasta,fasta_fai )


    
    
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


