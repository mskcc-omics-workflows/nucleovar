workflow FILTERING {
    take:
    vardict_output
    mutect_output

    main:
    BASIC_FILTERING_VARDICT()
    BASIC_FILTERING_MUTECT()


    emit:
    vardict_filtered_output = // channel: [ layout ]
    mutect_filtered_output = // channel: [ layout ]
}