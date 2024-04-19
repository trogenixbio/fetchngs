process UPDATE_USER_JSON {
    label 'process_low'

    conda "conda-forge::python=3.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path metadata_json
    path samplesheet_json
    val cloud_prefix
    val pub_internal
    val run_info_id

    output:
    path "metadata_update.json"    , emit: metadata_json_update
    path "samplesheet_update.json" , emit: samplesheet_json_update
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    update_user_json.py $metadata_json $samplesheet_json metadata_update.json samplesheet_update.json $cloud_prefix $pub_internal $run_info_id

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
