process BASESPACE_CLI {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'something' :
        'bensouthgate/illumina_bs_cli:latest' }"

    input:
    tuple val(meta), val(fastq)

    output:
    tuple val(meta), path("${meta.id}*.fastq.gz"), emit: fastq
    tuple val(meta), path("*md5")     , emit: md5
    path "*.json"                     , emit: json
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''

    def api_server = params.bs_api_server
    def access_token = params.bs_access_token

    def echo_md5_single = meta.md5_1 ? "echo '${meta.md5_1}  ${meta.id}.fastq.gz'": "echo 'No md5sum available for ${meta.id}.fastq.gz'"

    def echo_md5_1 = meta.md5_1 ? "echo '${meta.md5_1}  ${meta.id}_1.fastq.gz'": "echo 'No md5sum available for ${meta.id}_1.fastq.gz'"
    def echo_md5_2 = meta.md5_2 ? "echo '${meta.md5_2}  ${meta.id}_2.fastq.gz'": "echo 'No md5sum available for ${meta.id}_2.fastq.gz'"

    def md5_1 = meta.md5_1 ? "md5sum -c": "touch "
    def md5_2 = meta.md5_2 ? "md5sum -c": "touch "

    def suffix_1 = params.suffix_1 ?: "R1_001"
    def suffix_2 = params.suffix_2 ?: "R2_001"

    if (meta.single_end) {
        """
        bs download dataset --name ${fastq[0]} --api-server $api_server --access-token $access_token --extension=fastq.gz

        cp *fastq.gz ${meta.id}.fastq.gz
        $echo_md5_single > ${meta.id}.fastq.gz.md5
        $md5_1 ${meta.id}.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            BaseSpaceCLI: \$(bs --version | sed 's/BaseSpaceCLI //g')
        END_VERSIONS
        """
    } else {
        """
        bs download dataset --name ${fastq[0]} --api-server $api_server --access-token $access_token --extension=fastq.gz

        cp *${suffix_1}.fastq.gz ${meta.id}_1.fastq.gz
        $echo_md5_1 > ${meta.id}_1.fastq.gz.md5
        $md5_1 ${meta.id}_1.fastq.gz.md5

        cp *${suffix_2}.fastq.gz ${meta.id}_2.fastq.gz
        $echo_md5_2 > ${meta.id}_2.fastq.gz.md5
        $md5_2 ${meta.id}_2.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            BaseSpaceCLI: \$(bs --version | sed 's/BaseSpaceCLI //g')
        END_VERSIONS
        """
    }
}
