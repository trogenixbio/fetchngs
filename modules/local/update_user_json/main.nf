process UPDATE_USER_JSON {
    label 'process_low'

    conda "conda-forge::python=3.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path metadata_json
    val cloud_prefix

    output:
    path "metadata_update.json" , emit: metadata_json_update
    path "versions.yml"         , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = "${params.outdir}"
    """
    update_user_json.py $metadata_json $prefix metadata_update.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
