
process SRA_FASTQ_FTP {
    tag "$meta.id"
    label 'process_low'
    label 'error_retry'

    conda "conda-forge::wget=1.20.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.20.1' :
        'biocontainers/wget:1.20.1' }"

    input:
    tuple val(meta), val(fastq)

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
        wget \\
            $args \\
            -O ${meta.id}.fastq.gz \\
            ${fastq[0]}

        $echo_md5_single > ${meta.id}.fastq.gz.md5
        $md5_1 ${meta.id}.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            wget: \$(echo \$(wget --version | head -n 1 | sed 's/^GNU Wget //; s/ .*\$//'))
        END_VERSIONS
        """
    } else {
        """
        wget \\
            $args \\
            -O ${meta.id}_1.fastq.gz \\
            ${fastq[0]}

        $echo_md5_1 > ${meta.id}_1.fastq.gz.md5
        $md5_1 ${meta.id}_1.fastq.gz.md5

        wget \\
            $args \\
            -O ${meta.id}_2.fastq.gz \\
            ${fastq[1]}

        $echo_md5_2 > ${meta.id}_2.fastq.gz.md5
        $md5_2 ${meta.id}_2.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            wget: \$(echo \$(wget --version | head -n 1 | sed 's/^GNU Wget //; s/ .*\$//'))
        END_VERSIONS
        """
    }
}
