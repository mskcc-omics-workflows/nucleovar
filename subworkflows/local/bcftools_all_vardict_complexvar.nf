

//import bcftools subworkflow for vardict filtered regular
//include { TABIX_BGZIPTABIX } from '../../modules/nf-core/tabix/bgziptabix/main'
include { TABIX_BGZIP } from '../../modules/nf-core/tabix/bgzip/main' 
include { BCFTOOLS_INDEX } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_NORM } from '../../modules/nf-core/bcftools/norm/main' 
include { BCFTOOLS_SORT } from '../../modules/nf-core/bcftools/sort/main' 



workflow BCFTOOLS_ALL_VARDICT_COMPLEXVAR {
    take:
    vardict_filtered_vcf
    ref_fasta
    ref_fasta_index
    

    main:

    vardict_filtered_vcf
        .map { create_samplenames_for_bcftools_channel(it)  }
        .set{ sample_metamap }

    sample_metamap
        .combine( vardict_filtered_vcf )
        .set{ vardict_filtered_vcf_for_bcftools_ch }

    sample_metamap
        .combine( ref_fasta )
        .set{ meta_plus_fasta_ch }

    TABIX_BGZIP( vardict_filtered_vcf_for_bcftools_ch )

    compressed_vardict_filtered_vcf = TABIX_BGZIP.out.output
    

    BCFTOOLS_INDEX( compressed_vardict_filtered_vcf )
    vardict_filtered_index = BCFTOOLS_INDEX.out.tbi
    

    vardict_filtered_index
        .map{ sample,index -> index}
        .set{ vardict_filtered_index_isolated }
    
    compressed_vardict_filtered_vcf
        .combine( vardict_filtered_index_isolated )
        .set{ vardict_filtered_vcf_and_index_for_bcftools }


    BCFTOOLS_NORM( vardict_filtered_vcf_and_index_for_bcftools, meta_plus_fasta_ch )
    vardict_normalized_vcf = BCFTOOLS_NORM.out.vcf 
    

    sample_metamap
        .combine(vardict_normalized_vcf)
        .set{ meta_plus_vardict_normalized_vcf_ch }
    
    
    BCFTOOLS_SORT( meta_plus_vardict_normalized_vcf_ch )
    vardict_sorted_vcf = BCFTOOLS_SORT.out.vcf
    

    sample_metamap
        .combine(vardict_sorted_vcf)
        .set{ meta_plus_vardict_norm_and_sorted_vcf }

    
    vardict_filtered_index_isolated
        .map{ create_index_filenames(it) }
        .set{ vardict_filtered_index_final }

    

    
    emit:
    meta_plus_vardict_norm_and_sorted_vcf
    vardict_filtered_index_final
    
    
}

def create_samplenames_for_bcftools_channel( it ) {
    // create meta map
        def meta = [:]
        temp = file(it).getSimpleName()
        meta.id = temp[0..-19]
        inputs = meta
    }

def create_index_filenames( it ) {
    // create meta map
    // def newFileName = it.toString().replaceAll(/\.vcf.gz.tbi$/, ".tbi.gz")

    // // Copy the contents of the original file to the new file
    // def inputFile = file(it.toString())
    // def outputFile = file(newFileName)

    def finalFileName = "vcf_complexvar.vcf.gz.tbi"

    file(it).moveTo(file(finalFileName))
    //outputFile.text = inputFile.text
}