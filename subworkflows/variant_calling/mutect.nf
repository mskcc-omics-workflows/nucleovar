workflow MUTECT   {
    take:
    // initial inputs

    main:
    // differentiate if the vardict or mutect subworkflow will be called (??)
    MUTECT()
    MUTECT_FILTER()

    emit:
    mutect_filter_output // channel: [ layout ]
}