/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MULTIQC_MAPPINGS_CONFIG } from '../../modules/local/multiqc_mappings_config'
include { SRA_FASTQ_FTP           } from '../../modules/local/sra_fastq_ftp'
include { SRA_IDS_TO_RUNINFO      } from '../../modules/local/sra_ids_to_runinfo'
include { SRA_RUNINFO_TO_FTP      } from '../../modules/local/sra_runinfo_to_ftp'
include { ASPERA_CLI              } from '../../modules/local/aspera_cli'
include { SRA_TO_SAMPLESHEET      } from '../../modules/local/sra_to_samplesheet'
include { SAMPLE_TO_JSON          } from '../../modules/local/sample_to_json'
include { JSON_TO_METADATA        } from '../../modules/local/json_to_metadata'
include { CHECK_METADATA          } from '../../modules/local/check_metadata'
include { LOAD_USER_METADATA      } from '../../modules/local/load_user_metadata'
include { DOWNLOAD_USER_DATA      } from '../../modules/local/download_user_data'
include { UPDATE_USER_JSON        } from '../../modules/local/update_user_json'
include { SRA_FASTQ_FTP as SRA_FASTQ_FTP_INTERNAL } from '../../modules/local/sra_fastq_ftp'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS } from '../../subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools'
include { softwareVersionsToYAML                       } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SRA {

    take:
    ids // channel: [ ids ]

    main:
    ch_versions = Channel.empty()

    // Check if user supplied metadata sheet
    // If not assume SRA IDs supplied
    if (!params.metadata_sheet) {

        //
        // MODULE: Get SRA run information for public database ids
        //
        SRA_IDS_TO_RUNINFO (
            ids,
            params.ena_metadata_fields ?: ''
        )
        ch_versions = ch_versions.mix(SRA_IDS_TO_RUNINFO.out.versions.first())

        //
        // MODULE: Parse SRA run information, create file containing FTP links and read into workflow as [ meta, [reads] ]
        //

        SRA_RUNINFO_TO_FTP (
            SRA_IDS_TO_RUNINFO.out.json
        )
        ch_versions = ch_versions.mix(SRA_RUNINFO_TO_FTP.out.versions.first())

        SRA_RUNINFO_TO_FTP
            .out
            .json
            .splitJson()
            .map {
                meta ->
                    def meta_clone = meta.clone()
                    meta_clone.single_end = meta_clone.single_end.toBoolean()
                    return meta_clone
            }
            .unique()
            .set { ch_sra_metadata }

        if (!params.skip_fastq_download) {

            ch_sra_metadata
                .branch {
                    meta ->
                        def download_method = 'ftp'
                        // meta.fastq_aspera is a metadata string with ENA fasp links supported by Aspera
                            // For single-end: 'fasp.sra.ebi.ac.uk:/vol1/fastq/ERR116/006/ERR1160846/ERR1160846.fastq.gz'
                            // For paired-end: 'fasp.sra.ebi.ac.uk:/vol1/fastq/SRR130/020/SRR13055520/SRR13055520_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR130/020/SRR13055520/SRR13055520_2.fastq.gz'
                        if (meta.fastq_aspera && params.download_method == 'aspera') {
                            download_method = 'aspera'
                        }
                        if ((!meta.fastq_aspera && !meta.fastq_1) || params.download_method == 'sratools') {
                            download_method = 'sratools'
                        }

                        aspera: download_method == 'aspera'
                            return [ meta, meta.fastq_aspera.tokenize(';').take(2) ]
                        ftp: download_method == 'ftp'
                            return [ meta, [ meta.fastq_1, meta.fastq_2 ] ]
                        sratools: download_method == 'sratools'
                            return [ meta, meta.run_accession ]
                }
                .set { ch_sra_reads }

            //
            // MODULE: If FTP link is provided in run information then download FastQ directly via FTP and validate with md5sums
            //
            SRA_FASTQ_FTP (
                ch_sra_reads.ftp
            )
            ch_versions = ch_versions.mix(SRA_FASTQ_FTP.out.versions.first())

            //
            // SUBWORKFLOW: Download sequencing reads without FTP links using sra-tools.
            //
            FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS (
                ch_sra_reads.sratools,
                params.dbgap_key ? file(params.dbgap_key, checkIfExists: true) : []
            )
            ch_versions = ch_versions.mix(FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.versions.first())

            //
            // MODULE: If Aspera link is provided in run information then download FastQ directly via Aspera CLI and validate with md5sums
            //
            ASPERA_CLI (
                ch_sra_reads.aspera,
                'era-fasp'
            )
            ch_versions = ch_versions.mix(ASPERA_CLI.out.versions.first())

            // Isolate FASTQ channel which will be added to emit block
            SRA_FASTQ_FTP
                .out
                .fastq
                .mix(FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.reads)
                .mix(ASPERA_CLI.out.fastq)
                .map {
                    meta, fastq ->
                        def reads = fastq instanceof List ? fastq.flatten() : [ fastq ]
                        def meta_clone = meta.clone()

                        meta_clone.fastq_1 = reads[0] ? "${params.cloud_prefix}${params.outdir}/fastq/${reads[0].getName()}" : ''
                        meta_clone.fastq_2 = reads[1] && !meta.single_end ? "${params.cloud_prefix}${params.outdir}/fastq/${reads[1].getName()}" : ''

                        return meta_clone
                }
                .set { ch_sra_metadata }
        }

        // Add details to meta ready for json and samplesheet creation
        ch_sra_metadata
            .map {
                meta ->
                    def meta_clone = meta.clone()

                    meta_clone.remove("id")
                    meta_clone.remove("fastq_1")
                    meta_clone.remove("fastq_2")
                    meta_clone.remove("md5_1")
                    meta_clone.remove("md5_2")
                    meta_clone.remove("single_end")

                    // Add relevant fields to the beginning of the map
                    pipeline_map = [
                        sample  : "${meta.id.split('_')[0..-2].join('_')}",
                        fastq_1 : meta.fastq_1,
                        fastq_2 : meta.fastq_2
                    ]

                    // Add nf-core pipeline specific entries
                    def pipeline = params.nf_core_pipeline ?: ''
                    def rna_strandedness = params.nf_core_rnaseq_strandedness ?: 'auto'
                    if (pipeline) {
                        if (pipeline == 'rnaseq') {
                            pipeline_map << [ strandedness: rna_strandedness ]
                        } else if (pipeline == 'atacseq') {
                            pipeline_map << [ replicate: 1 ]
                        } else if (pipeline == 'taxprofiler') {
                            pipeline_map << [ fasta: '' ]
                        }
                    }
                    pipeline_map << meta_clone

                    return [ meta, pipeline_map ]
            }
            .set { ch_sra_metadata_pipeline_map }

        // Collect json file
        ch_sra_metadata_pipeline_map
            .map { it[1] }
            .collectFile( name:'testmeta.json') {
                "\t{\n" +
                '\t\t' + '"sample": ' + '"' + it.sample.toString() + '"' + ',' + '\n' +
                '\t\t' + '"fastq_1": ' + '"' + it.fastq_1.toString() + '"' + ',' + '\n' +
                '\t\t' + '"fastq_2": ' + '"' + it.fastq_2.toString() + '"' + ',' + '\n' +
                '\t\t' + '"run_accession": ' + '"' + it.run_accession.toString() + '"' + ',' + '\n' +
                '\t\t' + '"experiment_accession": ' + '"' + it.experiment_accession.toString() + '"' + ',' + '\n' +
                '\t\t' + '"sample_accession": ' + '"' + it.sample_accession.toString() + '"' + ',' + '\n' +
                '\t\t' + '"secondary_sample_accession": ' + '"' + it.secondary_sample_accession.toString() + '"' + ',' + '\n' +
                '\t\t' + '"study_accession": ' + '"' + it.study_accession.toString() + '"' + ',' + '\n' +
                '\t\t' + '"secondary_study_accession": ' + '"' + it.secondary_study_accession.toString() + '"' + ',' + '\n' +
                '\t\t' + '"submission_accession": ' + '"' + it.submission_accession.toString() + '"' + ',' + '\n' +
                '\t\t' + '"run_alias": ' + '"' + it.run_alias.toString() + '"' + ',' + '\n' +
                '\t\t' + '"experiment_alias": ' + '"' + it.experiment_alias.toString() + '"' + ',' + '\n' +
                '\t\t' + '"sample_alias": ' + '"' + it.sample_alias.toString() + '"' + ',' + '\n' +
                '\t\t' + '"study_alias": ' + '"' + it.study_alias.toString() + '"' + ',' + '\n' +
                '\t\t' + '"library_layout": ' + '"' + it.library_layout.toString() + '"' + ',' + '\n' +
                '\t\t' + '"library_selection": ' + '"' + it.library_selection.toString() + '"' + ',' + '\n' +
                '\t\t' + '"library_source": ' + '"' + it.library_source.toString() + '"' + ',' + '\n' +
                '\t\t' + '"library_strategy": ' + '"' + it.library_strategy.toString() + '"' + ',' + '\n' +
                '\t\t' + '"library_name": ' + '"' + it.library_name.toString() + '"' + ',' + '\n' +
                '\t\t' + '"instrument_model": ' + '"' + it.instrument_model.toString() + '"' + ',' + '\n' +
                '\t\t' + '"instrument_platform": ' + '"' + it.instrument_platform.toString() + '"' + ',' + '\n' +
                '\t\t' + '"base_count": ' + '"' + it.base_count.toString() + '"' + ',' + '\n' +
                '\t\t' + '"read_count": ' + '"' + it.read_count.toString() + '"' + ',' + '\n' +
                '\t\t' + '"tax_id": ' + '"' + it.tax_id.toString() + '"' + ',' + '\n' +
                '\t\t' + '"scientific_name": ' + '"' + it.scientific_name.toString() + '"' + ',' + '\n' +
                '\t\t' + '"sample_title": ' + '"' + it.sample_title.toString() + '"' + ',' + '\n' +
                '\t\t' + '"experiment_title": ' + '"' + it.experiment_title.toString() + '"' + ',' + '\n' +
                '\t\t' + '"study_title": ' + '"' + it.study_title.toString() + '"' + ',' + '\n' +
                '\t\t' + '"sample_description": ' + '"' + it.sample_description.toString() + '"' + ',' + '\n' +
                '\t\t' + '"fastq_md5": ' + '"' + it.fastq_md5.toString() + '"' + ',' + '\n' +
                '\t\t' + '"fastq_bytes": ' + '"' + it.fastq_bytes.toString() + '"' + ',' + '\n' +
                '\t\t' + '"fastq_ftp": ' + '"' + it.fastq_ftp.toString() + '"' + ',' + '\n' +
                '\t\t' + '"fastq_galaxy": ' + '"' + it.fastq_galaxy.toString() + '"' + ',' + '\n' +
                '\t\t' + '"fastq_aspera": ' + '"' + it.fastq_aspera.toString() + '"' + ',' + '\n' +
                '\t\t' + '"cell_type": ' + '"' + it.cell_type.toString() + '"' + ',' + '\n' +
                '\t\t' + '"tissue_type": ' + '"' + it.tissue_type.toString() + '"' + ',' + '\n' +
                '\t\t' + '"cell_line": ' + '"' + it.cell_line.toString() + '"' + ',' + '\n' +
                '\t\t' + '"strandedness": ' + '"' + it.strandedness.toString() + '"' + ',' + '\n' +
                '\t\t' + '"replicate": ' + '"' + it.replicate.toString() + '"' + ',' + '\n' +
                '\t\t' + '"fasta": ' + '"' + it.fasta.toString() + '"' + '\n' +
                "\t},\n"
            }
            .map {
                "[\n" +  it.text + "]"
            }
            .map {
                it.replace("},\n]", "}\n]")
            }
            .collectFile( name:'samplesheet.json', storeDir: "${params.outdir}/public/run_info/${params.pipeline_version}-${params.wf_timestamp}/metadata/samplesheet")
            .set { ch_samplesheet_json }

        //
        // MODULE: Stage FastQ files downloaded by SRA together and auto-create a samplesheet
        //

        SRA_TO_SAMPLESHEET (
            ch_sra_metadata_pipeline_map,
            params.sample_mapping_fields
        )

        // Merge samplesheets and mapping files across all samples
        SRA_TO_SAMPLESHEET
            .out
            .samplesheet
            .map { it[1] }
            .collectFile(name:'tmp_samplesheet.csv', newLine: true, keepHeader: true, sort: { it.baseName })
            .map { it.text.tokenize('\n').join('\n') }
            .collectFile(name:'samplesheet.csv', storeDir: "${params.outdir}/public/run_info/${params.pipeline_version}-${params.wf_timestamp}/metadata/samplesheet")
            .set { ch_samplesheet }

        //
        // MODULE: Convert to metadata schema structure
        // outputs a set of tsvs per metadata schema for users to double check

        JSON_TO_METADATA (
            ch_samplesheet_json,
            params.metadata_schema,
            params.mappings_json
        )
        ch_versions = ch_versions.mix(JSON_TO_METADATA.out.versions.first())

        //
        // MODULE: Manually check metadata
        // Reads in json, checks for things like Types and mandatory, and regex defined in schema

        CHECK_METADATA (
            JSON_TO_METADATA.out.metadata_json, // This will not be available to user supplied.
            params.metadata_schema
        )
        ch_versions = ch_versions.mix(CHECK_METADATA.out.versions.first())

        // Collect ID mappings for reference
        SRA_TO_SAMPLESHEET
            .out
            .mappings
            .map { it[1] }
            .collectFile(name:'tmp_id_mappings.csv', newLine: true, keepHeader: true, sort: { it.baseName })
            .map { it.text.tokenize('\n').join('\n') }
            .collectFile(name:'id_mappings.csv', storeDir: "${params.outdir}/public/run_info/${params.pipeline_version}-${params.wf_timestamp}/metadata/samplesheet")
            .set { ch_mappings }

        //
        // MODULE: Create a MutiQC config file with sample name mappings
        //
        ch_sample_mappings_yml = Channel.empty()
        if (params.sample_mapping_fields) {
            MULTIQC_MAPPINGS_CONFIG (
                ch_mappings
            )
            ch_versions = ch_versions.mix(MULTIQC_MAPPINGS_CONFIG.out.versions)
            ch_sample_mappings_yml = MULTIQC_MAPPINGS_CONFIG.out.yml
        }

    }

    //
    // User supplied data
    //

    // If user supplied create json from tsv or xlsx input file
    if (params.metadata_sheet) {

        //
        // MODULE: Load in user defined metadata
        //

        LOAD_USER_METADATA (
            params.metadata_sheet
        )
        ch_versions = ch_versions.mix(LOAD_USER_METADATA.out.versions.first())
        //
        // MODULE: Check Metadata
        //

        CHECK_METADATA (
            LOAD_USER_METADATA.out.metadata_json, // This will not be available to user supplied.
            params.metadata_schema
        )
        ch_versions = ch_versions.mix(CHECK_METADATA.out.versions.first())

        // Users may upload files directly to raw blob
        if (!params.skip_fastq_download) {
            // Stage in Fastq files
            LOAD_USER_METADATA
                .out
                .samplesheet_json
                .splitJson()
                .map {
                    meta ->
                    // add path(s) of the fastq file(s) to the meta map
                    def fastq_meta = []
                    if(params.download_method == "ftp") {
                        if (meta.fastq_2 == "") {
                            meta["single_end"] = true
                            fastq_meta = [ meta, [ meta.fastq_1 ] ]
                        } else {
                            meta["single_end"] = false
                            fastq_meta = [ meta, [ meta.fastq_1, meta.fastq_2 ] ]
                        }
                    }
                    else {
                        if (!file(meta.fastq_1).exists()) {
                            exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${meta.fastq_1}"
                        }
                        if (meta.fastq_2 == "") {
                            meta["single_end"] = true
                            fastq_meta = [ meta, [ file(meta.fastq_1) ] ]
                        } else {
                            if (!file(meta.fastq_2).exists()) {
                                exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${meta.fastq_2}"
                            }
                            meta["single_end"] = false
                            fastq_meta = [ meta, [ file(meta.fastq_1), file(meta.fastq_2) ] ]
                        }
                    }
                    return fastq_meta
                }
                .set { ch_user_metadata }

            if (params.download_method == "ftp") {
                // Stage in Fastq files
                SRA_FASTQ_FTP_INTERNAL (
                    ch_user_metadata
                )
                ch_versions = ch_versions.mix(SRA_FASTQ_FTP_INTERNAL.out.versions.first())
            } else {
                // Stage in Fastq files
                DOWNLOAD_USER_DATA (
                    ch_user_metadata
                )
                ch_versions = ch_versions.mix(DOWNLOAD_USER_DATA.out.versions.first())
            }

            // META DATA MUST BE UPDATED WITH NEW LOCATION
            UPDATE_USER_JSON (
                LOAD_USER_METADATA.out.metadata_json,
                params.cloud_prefix
            )
            ch_versions = ch_versions.mix(UPDATE_USER_JSON.out.versions.first())

            //
            // MODULE: Check Metadata
            //

            CHECK_METADATA (
                UPDATE_USER_JSON.out.metadata_json_update, // This will not be available to user supplied.
                params.metadata_schema
            )
            ch_versions = ch_versions.mix(CHECK_METADATA.out.versions.first())
        }
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/public/run_info/${params.pipeline_version}-${params.wf_timestamp}/pipeline_info", name: 'nf_core_fetchngs_software_mqc_versions.yml', sort: true, newLine: true)

    emit:
    samplesheet     = ch_samplesheet
    mappings        = ch_mappings
    sample_mappings = ch_sample_mappings_yml
    sra_metadata    = ch_sra_metadata
    versions        = ch_versions.unique()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
