
process SRA_RUNINFO_TO_FTP {

    conda "conda-forge::python=3.9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'biocontainers/python:3.9--1' }"

    input:
    path runinfo

    output:
    //path "*.tsv"       , emit: tsv
    path "*.json"       , emit: json
    path "versions.yml", emit: versions
    // sra_runinfo_to_ftp.py \\
    //     ${runinfo.join(',')} \\
    //     ${runinfo.toString().tokenize(".")[0]}.runinfo_ftp.tsv
    script:
    """
    sra_runinfo_to_ftp_json.py \\
        ${runinfo.join(',')} \\
        ${runinfo.toString().tokenize(".")[0]}.runinfo_ftp.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
