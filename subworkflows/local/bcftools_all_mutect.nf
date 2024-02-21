

//import bcftools subworkflow for mutect filtered regular
//include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIP } from '../../modules/nf-core/tabix/bgzip/main' 
include { BCFTOOLS_INDEX } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main' 
include { BCFTOOLS_SORT } from '../../modules/nf-core/bcftools/sort/main' 



workflow BCFTOOLS_ALL_MUTECT {
    take:
    mutect_filtered_vcf
    ref_fasta
    ref_fasta_index
    

    main:

    mutect_filtered_vcf
        .map { create_samplenames_for_bcftools_channel(it)  }
        .set{ sample_metamap }

    sample_metamap
        .combine( mutect_filtered_vcf )
        .set{ mutect_filtered_vcf_for_bcftools_ch }

    sample_metamap
        .combine( ref_fasta )
        .set{ meta_plus_fasta_ch }

    TABIX_BGZIP( mutect_filtered_vcf_for_bcftools_ch )

    compressed_mutect_filtered_vcf = TABIX_BGZIP.out.output
    

    BCFTOOLS_INDEX( compressed_mutect_filtered_vcf )
    mutect_filtered_index = BCFTOOLS_INDEX.out.tbi
    

    mutect_filtered_index
        .map{ sample,index -> index}
        .set{ mutect_filtered_index_isolated }
    
    compressed_mutect_filtered_vcf
        .combine( mutect_filtered_index_isolated )
        .set{ mutect_filtered_vcf_and_index_for_bcftools }


    BCFTOOLS_NORM( mutect_filtered_vcf_and_index_for_bcftools, meta_plus_fasta_ch )
    mutect_normalized_vcf = BCFTOOLS_NORM.out.vcf 
    

    sample_metamap
        .combine(mutect_normalized_vcf)
        .set{ meta_plus_mutect_normalized_vcf_ch }
    
    
    BCFTOOLS_SORT( meta_plus_mutect_normalized_vcf_ch )
    mutect_sorted_vcf = BCFTOOLS_SORT.out.vcf
    

    sample_metamap
        .combine(mutect_sorted_vcf)
        .set{ meta_plus_mutect_norm_and_sorted_vcf }


    mutect_filtered_index_isolated
        .map{ create_index_filenames(it) }
        .set{ mutect_filtered_index_final }

    

    
    emit:
    meta_plus_mutect_norm_and_sorted_vcf
    mutect_filtered_index_final
    



    

    
}

def create_samplenames_for_bcftools_channel( it ) {
    // create meta map
        def meta = [:]
        meta.id = file(it).getSimpleName()
        //meta.id = temp[0..-11]
        inputs = meta
    }

def create_index_filenames( it ) {

    def finalFileName = "vcf_mutect.vcf.gz.tbi"

    file(it).moveTo(file(finalFileName))

    
}


