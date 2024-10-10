process MUTECT_FILTER {
    tag "$meta.id"
    label 'process_single'

    container "ghcr.io/msk-access/postprocessing_variant_calls:postprocessing_test"

    input:
    tuple val(meta), path(maf), path(rules_json_file)


    output:
    path("*.maf"),                     emit: tagged_maf
    //path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''


    """

    pv maf tag by_rules --maf ${maf} --rules ${rules_json_file}

    """
}
