
process LOAD_USER_METADATA {
    tag "$metadata"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-693e24f156d01a5f55647120be99929b01b30949:609c862c3470382215fc1b2d9d8a4e9637b2e25f-0' :
        'biocontainers/mulled-v2-693e24f156d01a5f55647120be99929b01b30949:609c862c3470382215fc1b2d9d8a4e9637b2e25f-0' }"

    input:
    path input_metadata
    val is_excel

    output:
    path "metadata.json", emit: metadata

    script:
    def is_excel = is_excel ? "--is_excel": ""
    """
    load_user_metadata.py $is_excel $input_metadata metadata.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS

    """
}
