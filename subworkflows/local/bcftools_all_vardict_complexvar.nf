

//import bcftools subworkflow for vardict filtered regular
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main' 
include { BCFTOOLS_SORT } from '../../modules/nf-core/bcftools/sort/main' 



workflow BCFTOOLS_ALL_VARDICT_COMPLEXVAR {
    take:
    vardict_complexvar_filtered_vcf
    ref_fasta
    

    main:

    vardict_complexvar_filtered_vcf
        .map { create_samplenames_for_bcftools_channel(it)  }
        .set{ sample_metamap }

    sample_metamap
        .combine( vardict_complexvar_filtered_vcf )
        .set{ vardict_complexvar_filtered_vcf_for_bcftools_ch }


    TABIX_BGZIPTABIX( vardict_complexvar_filtered_vcf_for_bcftools_ch  )

    vardict_complexvar_filtered_vcf_and_index = TABIX_BGZIPTABIX.out.gz_tbi

    sample_metamap
        .combine( ref_fasta )
        .set{ meta_plus_fasta_ch }


    BCFTOOLS_NORM( vardict_complexvar_filtered_vcf_and_index, meta_plus_fasta_ch )
    vardict_complexvar_normalized_vcf = BCFTOOLS_NORM.out.vcf 

    sample_metamap
        .combine(vardict_complexvar_normalized_vcf)
        .set{ meta_plus_vardict_normalized_vcf_ch }
    

    
    BCFTOOLS_SORT( meta_plus_vardict_normalized_vcf_ch )
    vardict_complexvar_sorted_vcf = BCFTOOLS_SORT.out.vcf

    sample_metamap
        .combine(vardict_complexvar_sorted_vcf)
        .set{ meta_plus_vardict_complexvar_norm_and_sorted_vcf }

    
    emit:
    meta_plus_vardict_complexvar_norm_and_sorted_vcf



    

    
}

def create_samplenames_for_bcftools_channel( it ) {
    // create meta map
        def meta = [:]
        meta.id = file(it).getBaseName()
        inputs = meta
    }


