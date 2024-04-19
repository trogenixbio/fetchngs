//
// Subworkflow with functionality specific to the nf-core/fetchngs pipeline
//

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version             // boolean: Display version and exit
    help                // boolean: Display help text
    validate_params     // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs     // boolean: Do not use coloured log outputs
    nextflow_cli_args   //   array: List of positional nextflow CLI args
    outdir              //  string: The output directory where the results will be saved
    input               //  string: File containing SRA/ENA/GEO/DDBJ identifiers one per line to download their associated metadata and FastQ files
    ena_metadata_fields //  string: Comma-separated list of ENA metadata fields to fetch before downloading data
    metadata_sheet

    main:

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input ids.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Auto-detect input id type
    //

    if (!metadata_sheet) {
        ch_input = file(input)
        if (isSraId(ch_input)) {
            sraCheckENAMetadataFields(ena_metadata_fields)
        } else {
            error('Ids provided via --input not recognised please make sure they are either SRA / ENA / GEO / DDBJ ids!')
        }
        // Read in ids from --input file
        Channel
            .from(ch_input)
            .splitCsv(header:false, sep:'', strip:true)
            .map { it[0] }
            .unique()
            .set { ch_ids }
    } else {
        ch_ids = file(input)
    }

    emit:
    ids = ch_ids
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs)
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }

        sraCurateSamplesheetWarn()

        writeWorkflowSummary(summary_params, outdir)

    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Output workflow summary as standard
//
def writeWorkflowSummary(summary_params, outdir) {

    // Gather params from summary_params
    def summary = [:]
    for (group in summary_params.keySet()) {
        summary << summary_params[group]
    }

    // Define pipeline id variables
    def pipeline_version = summary["pipeline_version"]
    def wf_timestamp = summary["wf_timestamp"]

    // Define resume details
    def sessionID = workflow.sessionId
    def resumed = "FALSE"
    if (workflow.resume) {
        resumed = "RESUMED (SESSION ID ${sessionID})"
    }

    // Define completion statement
    def completion = "FAILED"
    if (workflow.success) {
        completion = "SUCCESS"
    }

    def jsonOutput = """{
    "completion": "${completion}",
    "exit-status": "${workflow.exitStatus}",
    "resume": "${resumed}",
    "session-id": "${sessionID}",
    "run-command": "${workflow.commandLine}",
    "pipeline-run-id": "${pipeline_version}-${wf_timestamp}",
    "pipeline-timestamp": "${wf_timestamp}",
    "pipeline-runinfo-outdir": "${outdir}",
    "user": "${workflow.userName}",
    "date-started": "${workflow.start}",
    "date-completed": "${workflow.complete}",
    "pipeline-script-file-path": "${workflow.scriptFile}",
    "pipeline-script-hash-id": "${workflow.scriptId}",
    "pipeline-git-repo-url": "${workflow.manifest.homePage}",
    "pipeline-git-branch-tag": "${workflow.manifest.version}",
    "nextflow-version": "${workflow.nextflow.version}",
    "nextflow-build": "${workflow.nextflow.build}",
    "nextflow-compile-timestamp": "${workflow.nextflow.timestamp}"\n}"""

    def summaryFilePathjson = "${outdir}/workflow_summary_${pipeline_version}-${wf_timestamp}.json"
    File summaryFilejson = new File(summaryFilePathjson)

    // Write the JSON string to the file
    try {
        summaryFilejson.withWriter('UTF-8') { writer ->
            writer.write(jsonOutput)
        }
    } catch (IOException e) {
        println "Error writing to file: ${e.message}"
    }
}

//
// Check if input ids are from the SRA
//
def isSraId(input) {
    def is_sra = false
    def total_ids = 0
    def no_match_ids = []
    def pattern = /^(((SR|ER|DR)[APRSX])|(SAM(N|EA|D))|(PRJ(NA|EB|DB))|(GS[EM]))(\d+)$/
    input.eachLine { line ->
        total_ids += 1
        if (!(line =~ pattern)) {
            no_match_ids << line
        }
    }

    def num_match = total_ids - no_match_ids.size()
    if (num_match > 0) {
        if (num_match == total_ids) {
            is_sra = true
        } else {
            error("Mixture of ids provided via --input: ${no_match_ids.join(', ')}\nPlease provide either SRA / ENA / GEO / DDBJ ids!")
        }
    }
    return is_sra
}

//
// Check and validate parameters
//
def sraCheckENAMetadataFields(ena_metadata_fields) {
    // Check minimal ENA fields are provided to download FastQ files
    def valid_ena_metadata_fields = ['run_accession', 'experiment_accession', 'library_layout', 'fastq_ftp', 'fastq_md5']
    def actual_ena_metadata_fields = ena_metadata_fields ? ena_metadata_fields.split(',').collect{ it.trim().toLowerCase() } : valid_ena_metadata_fields
    if (!actual_ena_metadata_fields.containsAll(valid_ena_metadata_fields)) {
        error("Invalid option: '${ena_metadata_fields}'. Minimally required fields for '--ena_metadata_fields': '${valid_ena_metadata_fields.join(',')}'")
    }
}

//
// Print a warning after pipeline has completed
//
def sraCurateSamplesheetWarn() {
    log.warn "=============================================================================\n" +
        "  Please double-check the samplesheet and metadata that has been auto-created\n" +
        "  by the pipeline before any system wide metadata update is run.\n\n" +
        "  Public databases don't reliably hold information such as strandedness\n" +
        "  information, controls etc\n\n" +
        "  All of the sample metadata obtained from the ENA has been appended\n" +
        "  as additional columns to help you manually curate the samplesheet/metadata before\n" +
        "  running nf-core/other pipelines.\n" +
        "==================================================================================="
}
