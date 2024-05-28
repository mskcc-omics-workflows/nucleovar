

//import bcftools subworkflow for vardict filtered regular
//include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { VCF2MAF } from '../../modules/nf-core/vcf2maf/main'
include { TRACEBACK } from '../../subworkflows/msk/traceback/main'




workflow MODULE4 {
    take:
    //vcf
    case_bams_for_traceback
    aux_bams
    fasta
    fasta_fai 
    maf_header
    maf_header_genotype
    //rules_json 
    

    main:
    // temporary input maf for testing purposes 
    input_maf = Channel.fromPath("/work/access/production/data/small_variants/C-PR83CF/C-PR83CF-L004-d04/current/C-PR83CF-L004-d04.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots.maf")
    mafs = Channel.from([patient:'test',id:"C-2HXC96-P001-d01.DONOR22-TP.combined-variants"]).merge(input_maf)


    //run vcf2maf perl module (placeholder at the moment while temporarily testing out traceback)
    // vcf.combine(fasta).set{ input_for_vcf2maf }
    // VCF2MAF( input_for_vcf2maf )
    // maf = VCF2MAF.out.maf

    // TAG_BY_ACCESS( maf,rules_json )


    // standard channel input (alternative version if using matched bams)
    // matched_bams
    //     .map{ meta,bam1,bai1,bam2,bai2 -> tuple([patient:'test',id:'C-2HXC96-P001-d01.DONOR22-TP.combined-variants'],bam1,bai1,[],[],[],[])}
    //     .view()
    //     .set{ standard_bams_for_traceback }


    case_bams_for_traceback.mix(aux_bams).set{ test }
    
    // simplex/duplex channel input 
    
    
    TRACEBACK( test,mafs,[initial:file(maf_header),genotype:file(maf_header_genotype)],fasta,fasta_fai )
    TRACEBACK.out.individual_genotyped_mafs.view()
    TRACEBACK.out.genotyped_maf.view()

    
    
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


