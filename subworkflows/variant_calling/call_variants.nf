workflow CALL_VARIANTS   {
    take:
    // initial inputs

    main:
    // differentiate if the vardict or mutect subworkflow will be called (??)
    VARDICTJAVA()
    MUTECT()

    emit:
    vardict_output // channel: [ layout ]
    mutect_output // channel: [ layout ]
}