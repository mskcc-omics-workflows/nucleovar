process VARDICT_FILTER {
    label 'vardict_filter'

    // conda "conda-forge::python=3.8.3"
    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/python:3.8.3' :
    //     'biocontainers/python:3.8.3' }"

    input:
    path vardict_java_output

    output:
    path '*.vcf'       , emit: vardict_filtered_output
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in msk/nucleovar/bin/
    """
    check_samplesheet.py $samplesheet 



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
