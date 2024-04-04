process DOWNLOAD_USER_DATA {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::python=3.9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    tuple val(meta), path(fastq, stageAs: 'fastq/*')

    output:
    tuple val(meta), path("*fastq.gz"), emit: fastq
    tuple val(meta), path("*md5")     , emit: md5
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''

    def echo_md5_single = meta.md5_1 ? "echo '${meta.md5_1}  ${meta.id}.fastq.gz'": "echo 'No md5sum available for ${meta.id}.fastq.gz'"

    def echo_md5_1 = meta.md5_1 ? "echo '${meta.md5_1}  ${meta.id}_1.fastq.gz'": "echo 'No md5sum available for ${meta.id}_1.fastq.gz'"
    def echo_md5_2 = meta.md5_2 ? "echo '${meta.md5_2}  ${meta.id}_2.fastq.gz'": "echo 'No md5sum available for ${meta.id}_2.fastq.gz'"

    def md5_1 = meta.md5_1 ? "md5sum -c": "touch "
    def md5_2 = meta.md5_2 ? "md5sum -c": "touch "

    if (meta.single_end) {
        """
        $echo_md5_single > ${meta.id}.fastq.gz.md5
        cp $fastq ${meta.id}.fastq.gz
        $md5_1 ${meta.id}.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    } else {
        """
        $echo_md5_1 > ${meta.id}_1.fastq.gz.md5
        cp ${fastq[0]} ${meta.id}_1.fastq.gz
        $md5_1 ${meta.id}_1.fastq.gz.md5

        $echo_md5_2 > ${meta.id}_2.fastq.gz.md5
        cp ${fastq[1]} ${meta.id}_2.fastq.gz
        $md5_2 ${meta.id}_2.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
    }
}
