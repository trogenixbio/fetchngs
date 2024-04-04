process JSON_TO_SAMPLESHEET {
    tag "$sample_json"
    label 'error_retry'

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path sample_json
    val nf_core_pipeline

    output:
    path "samplesheet.csv"      , emit: samplesheet
    path "versions.yml"         , emit: versions

    script:
    def pipeline = "${nf_core_pipeline}" ?: ''
    def include  = ''
    def exclude  = 'strandedness,replicate,fasta'

    if (pipeline) {
        if (pipeline == 'rnaseq') {
            include = "strandedness"
            exclude = "replicate,fasta"
        } else if (pipeline == 'atacseq') {
            include = "replicate"
            exclude = "strandedness,fasta"
        } else if (pipeline == 'taxprofiler') {
            include = "fasta"
            exclude = "replicate,strandedness"
        }
    }

    """
    json_to_samplesheet.py \\
        ${sample_json} \\
        samplesheet.csv \\
        --include ${include} \\
        --exclude ${exclude}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
