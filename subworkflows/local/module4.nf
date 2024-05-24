

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
    //rules_json 
    

    main:
    // temporary input maf for testing purposes 
    input_maf = Channel.fromPath("/work/access/production/data/small_variants/C-PR83CF/C-PR83CF-L004-d04/current/C-PR83CF-L004-d04.DONOR22-TP.combined-variants.vep_keptrmv_taggedHotspots.maf")


    

    donor10_simplex_bam = file("/home/naidur/access_pipeline/inputs/DONOR10-T_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam")
    donor10_simplex_bai = file("/home/naidur/access_pipeline/inputs/DONOR10-T_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bai")

    donor10_duplex_bam = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_duplex_bams_dmp/versions/v2.0/DONOR10-T_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam")
    donor10_duplex_bai = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_duplex_bams_dmp/versions/v2.0/DONOR10-T_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bai")

    donor38_duplex_bam = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_duplex_bams_dmp/versions/v2.0/DONOR38-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bam")
    donor38_duplex_bai = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_duplex_bams_dmp/versions/v2.0/DONOR38-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-duplex.bai")

    donor38_simplex_bam = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_simplex_bams_dmp/versions/v2.0/DONOR38-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bam")
    donor38_simplex_bai = file("/juno/cmo/access/production/resources/msk-access/v1.0/novaseq_curated_simplex_bams_dmp/versions/v2.0/DONOR38-TP_cl_aln_srt_MD_IR_FX_BR__aln_srt_IR_FX-simplex.bai")




    //simplex_bam1.combine(simplex_bai1).set{test_simplex_bams}

    //run vcf2maf perl module (placeholder at the moment while temporarily testing out traceback)
    // vcf.combine(fasta).set{ input_for_vcf2maf }
    // VCF2MAF( input_for_vcf2maf )
    // maf = VCF2MAF.out.maf

    // TAG_BY_ACCESS( maf,rules_json )


    emptyLists = Channel.from([], [], [], [])
    def bamnames = [patient:'test',id:'sample']
    names = Channel.create(bamnames)

    
    // standard channel input (alternative version if using matched bams)
    // matched_bams
    //     .map{ meta,bam1,bai1,bam2,bai2 -> tuple([patient:'test',id:'C-2HXC96-P001-d01.DONOR22-TP.combined-variants'],bam1,bai1,[],[],[],[])}
    //     .view()
    //     .set{ standard_bams_for_traceback }



    
    
    case_bams_for_traceback.mix(aux_bams).set{ test }
    // duplex_bams
    //     .map{ meta,bam1,bai1,bam2,bai2 -> tuple([patient:'test',id:'C-2HXC96-P001-d01.DONOR22-TP.combined-variants'],
    //         [],
    //         [],
    //         bam1,
    //         bai1,
    //         normal_sample_simplex_bam,
    //         normal_sample_simplex_bai)}
    //     .set{ duplex_simplex_bams_for_traceback }

    
    // simplex/duplex channel input 
    mafs = Channel.from([patient:'test',id:"C-2HXC96-P001-d01.DONOR22-TP.combined-variants"]).merge(input_maf)
    
    // fasta = file("/juno/cmo/access/production/resources/reference/current/Homo_sapiens_assembly19.fasta")
    // fasta_fai = file("/juno/cmo/access/production/resources/reference/current/Homo_sapiens_assembly19.fai")
    TRACEBACK( test,mafs,[initial:file('/home/naidur/access_pipeline/inputs/maf_header.txt'),genotype:file('/home/naidur/access_pipeline/inputs/maf_header_genotype.txt')],fasta,fasta_fai )
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


