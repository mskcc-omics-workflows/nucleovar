workflow VARDICT   {
    take:
    tumor_sample_name
    normal_sample_name
    bams
    bed
    fasta   

    main:
    VARDICTJAVA()
    VARDICT_FILTER()

    emit:
    vardict_filter_output // channel: [ layout ]
}