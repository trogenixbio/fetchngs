process SAMPLEJSON_TO_METADATA {
    tag "$sample_json"
    label 'error_retry'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path sample_json
    path metadata_schema
    path mappings_json

    output:
    path "metadata.json"        , emit: metadata_json
    path "*.tsv"                , emit: metadata_tsv
    path "versions.yml"         , emit: versions

    script:
    """
    samplejson_to_metadata.py \\
        ${sample_json} \\
        ${metadata_schema} \\
        metadata.json \\
        ${mappings_json}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
