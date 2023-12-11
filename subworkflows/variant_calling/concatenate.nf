workflow CONCATENATE {
    take:
    vardict_filtering_output
    mutect_filtering_output

    main:
    
    // join the vardict and mutect filtering outputs using native nextflow 

    BCFTOOLS_CONCAT()
    ANNOTATE_CONCAT()


    emit:
    final_concatenated_output  // channel: [ layout ]
}