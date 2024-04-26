

//import bcftools subworkflow for vardict filtered regular
include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { BGZIP } from '../../modules/nf-core/tabix/bgzip/main' 
include { BCFTOOLS_INDEX } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main' 
include { BCFTOOLS_SORT } from '../../modules/nf-core/bcftools/sort/main' 
include { BCFTOOLS_CONCAT } from '../../modules/nf-core/bcftools/concat/main' 
//include { BCFTOOLS_ALL_VARDICT     } from '../subworkflows/local/bcftools_all_vardict'

workflow BGZIP_INDEX_STANDARD_ {
    take:
    vardict_filtered_vcf
    ref_fasta
    ref_fasta_index
    

    main:
    vardict_filtered_vcf
        .map { create_samplenames_for_bcftools_standard(it)  }
        .set{ sample_metamap }

    sample_metamap
        .combine( vardict_filtered_vcf )
        .set{ vardict_filtered_vcf_for_bcftools_ch }

    sample_metamap
        .combine( ref_fasta )
        .set{ meta_plus_fasta_ch }

    
    BCFTOOLS_INDEX( vardict_filtered_vcf_for_bcftools_ch )
    sample_vcf_and_index = BCFTOOLS_INDEX.out.tbi

    sample_vcf_and_index
        .map{ sample,vcf,index -> sample}
        .set{ sample_metamap }
    
    sample_vcf_and_index
        .map{ sample,vcf,index -> index}
        .set{ isolated_index }

    sample_metamap
        .combine(ref_fasta)
        .set{ meta_plus_fasta_ch }

    emit:
    sample_vcf_and_index
    meta_plus_fasta_ch

}


workflow BGZIP_INDEX_COMPLEXVAR_ {
    take:
    vardict_filtered_vcf
    ref_fasta
    ref_fasta_index
    

    main:
    vardict_filtered_vcf
        .map { create_samplenames_for_bcftools_complexvar(it)  }
        .set{ sample_metamap }

    sample_metamap
        .combine( vardict_filtered_vcf )
        .set{ vardict_filtered_vcf_for_bcftools_ch }

    sample_metamap
        .combine( ref_fasta )
        .set{ meta_plus_fasta_ch }

    
    BCFTOOLS_INDEX( vardict_filtered_vcf_for_bcftools_ch )
    sample_vcf_and_index = BCFTOOLS_INDEX.out.tbi

    sample_vcf_and_index
        .map{ sample,vcf,index -> sample}
        .set{ sample_metamap }
    
    sample_vcf_and_index
        .map{ sample,vcf,index -> index}
        .set{ isolated_index }

    sample_metamap
        .combine(ref_fasta)
        .set{ meta_plus_fasta_ch }
    
    emit:
    sample_vcf_and_index
    meta_plus_fasta_ch

}

workflow NORM_STANDARD_ {
    take:
    sample_vcf_and_index
    meta_plus_fasta_ch

    main:
    
    BCFTOOLS_NORM( sample_vcf_and_index, meta_plus_fasta_ch )
    vardict_normalized_vcf = BCFTOOLS_NORM.out.vcf

    sample_vcf_and_index.map{ sample,vcf,index -> sample}.set{ sample_metamap }

    sample_metamap
        .combine(vardict_normalized_vcf)
        .set{ meta_plus_vardict_normalized_vcf_ch }


    emit:
    meta_plus_vardict_normalized_vcf_ch

}

workflow NORM_COMPLEXVAR_ {
    take:
    sample_vcf_and_index
    meta_plus_fasta_ch

    main:
    BCFTOOLS_NORM( sample_vcf_and_index, meta_plus_fasta_ch )
    vardict_normalized_vcf = BCFTOOLS_NORM.out.vcf

    sample_vcf_and_index.map{ sample,vcf,index -> sample}.set{ sample_metamap }


    sample_metamap
        .combine(vardict_normalized_vcf)
        .set{ meta_plus_vardict_normalized_vcf_ch }
    
    emit:
    meta_plus_vardict_normalized_vcf_ch
    
}



workflow BCFTOOLS_VARDICT {
    take:
    vardict_filtered_vcf_complexvar
    vardict_filtered_vcf_std
    ref_fasta
    ref_fasta_index

    

    main:
    
    BGZIP_INDEX_STANDARD_( vardict_filtered_vcf_std,ref_fasta,ref_fasta_index )
    BGZIP_INDEX_COMPLEXVAR_( vardict_filtered_vcf_complexvar,ref_fasta,ref_fasta_index )

    standard_sample_vcf_and_index = BGZIP_INDEX_STANDARD_.out.sample_vcf_and_index
    standard_meta_plus_fasta = BGZIP_INDEX_STANDARD_.out.meta_plus_fasta_ch

    complexvar_sample_vcf_and_index = BGZIP_INDEX_COMPLEXVAR_.out.sample_vcf_and_index
    complexvar_meta_plus_fasta = BGZIP_INDEX_COMPLEXVAR_.out.meta_plus_fasta_ch
    
    
    NORM_STANDARD_( standard_sample_vcf_and_index,standard_meta_plus_fasta )
    NORM_COMPLEXVAR_( complexvar_sample_vcf_and_index,standard_meta_plus_fasta )

    standard_norm_sorted_vcf = NORM_STANDARD_.out.meta_plus_vardict_normalized_vcf_ch
    
    complexvar_norm_sorted_vcf = NORM_COMPLEXVAR_.out.meta_plus_vardict_normalized_vcf_ch

    standard_sample_vcf_and_index.map{sample,vcf,index -> sample}.set{sample_metamap_standard}
    complexvar_sample_vcf_and_index.map{sample,vcf,index -> sample}.set{sample_metamap_complexvar}



    standard_sample_vcf_and_index.map{sample,vcf,index -> index}.set{index_standard}
    complexvar_sample_vcf_and_index.map{sample,vcf,index -> index}.set{index_complexvar}

    standard_norm_sorted_vcf.map{ sample,vcf -> vcf}.set{vcf_standard}
    complexvar_norm_sorted_vcf.map{ sample,vcf -> vcf}.set{vcf_complexvar}

    
    
    sample_metamap_standard
        .combine(vcf_standard)
        .combine(vcf_complexvar)
        .set{ inputs_for_bcftools_concat_ch }
    
    
    BCFTOOLS_CONCAT( inputs_for_bcftools_concat_ch )
    vardict_concat_vcf = BCFTOOLS_CONCAT.out.concat_vcf

    standard_sample_vcf_and_index.map{ id,vcf,index -> index}.set{ vardict_index }

    emit:
    vardict_concat_vcf
    vardict_index
    

}

def create_samplenames_for_bcftools_standard( it ) {
    // create meta map
        def meta = [:]
        meta.id = file(it).getSimpleName()
        inputs = meta
    }

def create_samplenames_for_bcftools_complexvar( it ) {
    // create meta map
        def meta = [:]
        temp = file(it).getSimpleName()
        meta.id = temp[0..-19] + "_STDfilter_complex"
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