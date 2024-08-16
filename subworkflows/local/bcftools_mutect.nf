

//import bcftools subworkflow for vardict filtered regular
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BGZIP } from '../../modules/nf-core/tabix/bgzip/main' 
include { BCFTOOLS_INDEX } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main' 
include { BCFTOOLS_SORT } from '../../modules/nf-core/bcftools/sort/main' 
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main' 
//include { BCFTOOLS_ALL_VARDICT     } from '../subworkflows/local/bcftools_all_vardict'

workflow BGZIP_INDEX {
    take:
    mutect_filtered_vcf
    ref_fasta
    ref_fasta_index
    

    main:
    mutect_filtered_vcf
        .map { create_samplenames_for_bcftools_mutect(it)  }
        .set{ sample_metamap }

    

    sample_metamap
        .combine( mutect_filtered_vcf )
        .set{ mutect_filtered_vcf_for_bcftools_ch }

    sample_metamap
        .combine( ref_fasta )
        .set{ meta_plus_fasta_ch }

    
    BCFTOOLS_INDEX( mutect_filtered_vcf_for_bcftools_ch )
    sample_vcf_and_index = BCFTOOLS_INDEX.out.tbi

    sample_vcf_and_index
        .map{ sample,vcf,index -> sample}
        .set{ sample_metamap }

    sample_metamap
        .combine(ref_fasta)
        .set{ meta_plus_fasta_ch }

    emit:
    sample_vcf_and_index
    meta_plus_fasta_ch

}



workflow BCFTOOLS_MUTECT {
    take:
    mutect_filtered_vcf
    ref_fasta
    ref_fasta_index

    

    main:

    BGZIP_INDEX( mutect_filtered_vcf,ref_fasta,ref_fasta_index )

    standard_sample_vcf_and_index = BGZIP_INDEX.out.sample_vcf_and_index
    standard_meta_plus_fasta = BGZIP_INDEX.out.meta_plus_fasta_ch
    
    
    BCFTOOLS_NORM( standard_sample_vcf_and_index,standard_meta_plus_fasta )

    standard_norm_sorted_vcf = BCFTOOLS_NORM.out.vcf
    

    standard_sample_vcf_and_index.map{sample,vcf,index -> index}.set{mutect_index}
    
    

    
    emit:
    standard_norm_sorted_vcf
    mutect_index
    

}

def create_samplenames_for_bcftools_mutect( it ) {
    // create meta map
        def meta = [:]
        meta.id = file(it).getBaseName()
        inputs = meta
    }

def create_index_filename_standard( it ) {

    def finalFileName = "vcf_standard.tbi.gz"

    file(it).moveTo(file(finalFileName))
}

def create_index_filename_complexvar( it ) {

    def finalFileName = "vcf_complexvar.tbi.gz"

    file(it).moveTo(file(finalFileName))
}