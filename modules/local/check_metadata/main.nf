process CHECK_METADATA {
    tag "$id"
    label 'error_retry'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path metadata_json
    path metadata_schema

    output:
    path "valid.metadata.json"  , emit: json
    path "versions.yml"         , emit: versions

    script:
    """
    check_schema.py \\
        ${metadata_json} \\
        ${metadata_schema}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
