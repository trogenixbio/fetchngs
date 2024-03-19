process METADATA_MERGE {
    tag "$id"
    label 'error_retry'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path big_metadata
    path valid_metadata_json

    output:
    path "valid.metadata.json"  , emit: json
    path "versions.yml"         , emit: versions

    script:
    """
    merge_metadata.py \\
        ${big_metadata} \\
        ${valid_metadata_json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
