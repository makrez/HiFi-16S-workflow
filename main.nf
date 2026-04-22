/*
===============================================================================

Minimal native preprocessing workflow for PacBio HiFi 16S.
Current scope:
- input loading
- metadata inspection
- FASTQ QC and filtering
- optional primer trimming
- QC aggregation

===============================================================================

Author: Refactor draft
Updated: 2026-04-21
*/
import groovy.yaml.YamlSlurper

nextflow.enable.dsl = 2

version = "0.1.0"

include { 
    download_gtdb_db
    download_silva_db
    download_gg2_db } from './modules/utils'

include {
    inspect_metadata
    QC_fastq
    cutadapt
    summarize_cutadapt
    QC_fastq_post_trim
    collect_QC
    collect_QC_skip_cutadapt
} from './modules/qc'

include {
    dada2_filter_ccs
    learn_errors
    dada2_denoise_independent
    dada2_make_seqtab
    dada2_remove_chimeras
    dada2_filter_asvs
    dada2_stats
    dada2_final_stats
} from './modules/dada2'

// ----------------------------------------------------------------------------
// Help text
// ----------------------------------------------------------------------------
def helpMessage() {
    return """
    Minimal native preprocessing workflow for PacBio HiFi 16S

    Required parameters:
      --input                Samples TSV with columns: sample, fastq
      --metadata             Metadata TSV with at least sample_name column

    Optional parameters:
      --outdir               Output directory [default: results]
      --filterQ              Minimum read quality filter [default: 20]
      --downsample           Downsample reads per sample, 0 disables [default: 0]
      --skip_primer_trim     Skip cutadapt trimming [default: false]
      --front_p              Forward primer sequence
      --adapter_p            Reverse primer sequence
      --download_db          Run database download workflow [default: false]
      --download_targets     Comma-separated database names, e.g. gtdb or silva,gg2
      -profile               standard / conda / docker / singularity

    Example:
      nextflow run main.nf \\
        --input test/samples.tsv \\
        --metadata test/metadata.tsv \\
        --outdir results_test \\
        -profile conda
    """.stripIndent()
}

// ----------------------------------------------------------------------------
// Early exits
// ----------------------------------------------------------------------------
if (params.help) {
    exit 0, helpMessage()
}

if (params.version) {
    exit 0, version
}

def db_manifest_file = file("${projectDir}/conf/databases.yml")
if (!db_manifest_file.exists()) {
    exit 1, "Database manifest not found: ${db_manifest_file}"
}

def db_manifest = new YamlSlurper().parse(db_manifest_file)

// ----------------------------------------------------------------------------
// Parameter parsing / defaults
// ----------------------------------------------------------------------------
def requested_dbs = params.download_targets ? params.download_targets.tokenize(',')*.trim() : []

def n_sample = null
def dynamic_min_asv_totalfreq = null
def dynamic_min_asv_sample = null
def trim_cutadapt = null
def log_text = null

// ----------------------------------------------------------------------------
// Validation
// ----------------------------------------------------------------------------
if (params.download_db) {
    if (!params.download_targets) {
        exit 1, "When --download_db true, you must also provide --download_targets"
    }

    def valid_dbs = ['silva', 'gtdb', 'gg2']

    requested_dbs.each { db ->
        if (!valid_dbs.contains(db)) {
            exit 1, "Invalid database '${db}'. Allowed values: ${valid_dbs.join(', ')}"
        }
        if (!db_manifest[db]) {
            exit 1, "Database '${db}' is missing from conf/databases.yml"
        }
    }

    if (requested_dbs.contains('silva') && !params.silva_dir) {
        exit 1, "Missing --silva_dir for SILVA download"
    }
    if (requested_dbs.contains('gtdb') && !params.gtdb_dir) {
        exit 1, "Missing --gtdb_dir for GTDB download"
    }
    if (requested_dbs.contains('gg2') && !params.gg2_dir) {
        exit 1, "Missing --gg2_dir for GG2 download"
    }
}
else {
    if (!params.input) {
        exit 1, "Missing required parameter: --input"
    }

    if (!params.metadata) {
        exit 1, "Missing required parameter: --metadata"
    }

    if (!file(params.input).exists()) {
        exit 1, "Input sample sheet not found: ${params.input}"
    }

    if (!file(params.metadata).exists()) {
        exit 1, "Metadata file not found: ${params.metadata}"
    }

    n_sample = file(params.input).countLines() - 1
    if (n_sample < 1) {
        exit 1, "Input sample sheet appears empty: ${params.input}"
    }

    if (n_sample == 1) {
        dynamic_min_asv_totalfreq = 0
        dynamic_min_asv_sample = 0
        log.info("Only 1 sample detected → disabling ASV filtering")
    } else {
        dynamic_min_asv_totalfreq = params.min_asv_totalfreq ?: 0
        dynamic_min_asv_sample = params.min_asv_sample ?: 0
    }

    trim_cutadapt = params.skip_primer_trim ? "No" : "Yes"

    log_text = """
Parameters for native pb-16S preprocessing
=========================================
Number of samples in samples TSV: ${n_sample}
Metadata file: ${params.metadata}
Filter input reads above Q: ${params.filterQ}
Downsample reads per sample (0 = disabled): ${params.downsample}
Trim primers with cutadapt: ${trim_cutadapt}
Forward primer: ${params.front_p}
Reverse primer: ${params.adapter_p}
Output directory: ${params.outdir}
Execution profile uses conda: ${params.enable_conda}
Execution profile uses containers: ${params.enable_container}
Pipeline version: ${version}
""".stripIndent()

    log.info(log_text)
}

// ----------------------------------------------------------------------------
// Workflows
// ----------------------------------------------------------------------------
workflow pb16S_preprocess {

    sample_sheet_ch = Channel.fromPath(params.input)
    metadata_ch = Channel.fromPath(params.metadata)

    inspect_metadata(sample_sheet_ch, metadata_ch)

    sample_ch = Channel
        .fromPath(params.input)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            if (!row['sample-id'] || !row['absolute-filepath']) {
                throw new IllegalArgumentException(
                    "Input TSV must contain columns 'sample-id' and 'absolute-filepath'"
                )
            }
            tuple(row['sample-id'] as String, file(row['absolute-filepath'] as String))
        }

    QC_fastq(sample_ch)

    def reads_for_dada2

    if (params.skip_primer_trim) {
        collect_QC_skip_cutadapt(
            QC_fastq.out.all_seqkit_stats.collect(),
            QC_fastq.out.all_seqkit_summary.collect()
        )

        reads_for_dada2 = QC_fastq.out.filtered_fastq
    }
    else {
        cutadapt(
            QC_fastq.out.filtered_fastq,
            params.front_p,
            params.adapter_p
        )

        summarize_cutadapt(cutadapt.out.cutadapt_report)

        QC_fastq_post_trim(cutadapt.out.cutadapt_fastq)

        collect_QC(
            QC_fastq.out.all_seqkit_stats.collect(),
            QC_fastq.out.all_seqkit_summary.collect(),
            summarize_cutadapt.out.summary_tocollect.collect(),
            QC_fastq_post_trim.out.all_seqkit_stats.collect()
        )

        reads_for_dada2 = cutadapt.out.cutadapt_fastq
    }

    def filter_script = file("${projectDir}/scripts/filter_and_trim.R", checkIfExists: true)

    dada2_filter_ccs(
        reads_for_dada2,
        filter_script
    )

    def learn_errors_script = file("${projectDir}/scripts/learn_errors.R", checkIfExists: true)

    learn_errors(
        dada2_filter_ccs.out.filtered_fastq.map { sampleID, fastq -> fastq }.collect(),
        learn_errors_script
    )

    def denoise_script = file("${projectDir}/scripts/denoise_independent.R", checkIfExists: true)

    dada2_denoise_independent(
        dada2_filter_ccs.out.filtered_fastq,
        learn_errors.out.error_model,
        denoise_script
    )

    def make_seqtab_script = file("${projectDir}/scripts/make_seqtab.R", checkIfExists: true)

    dada2_make_seqtab(
        dada2_denoise_independent.out.dada_rds.map { sampleID, rds -> rds }.collect(),
        make_seqtab_script
    )

    def remove_chimeras_script = file("${projectDir}/scripts/remove_chimeras.R", checkIfExists: true)

    dada2_remove_chimeras(
        dada2_make_seqtab.out.seqtab_rds,
        remove_chimeras_script
    )

    def filter_asvs_script = file("${projectDir}/scripts/filter_asvs.R", checkIfExists: true)

    dada2_filter_asvs(
        dada2_remove_chimeras.out.seqtab_nochim_rds,
        filter_asvs_script,
        dynamic_min_asv_totalfreq,
        dynamic_min_asv_sample
    )

    def dada2_stats_script = file("${projectDir}/scripts/dada2_stats.R", checkIfExists: true)

    dada2_stats(
        dada2_filter_asvs.out.seqtab_filtered_rds,
        metadata_ch,
        dada2_stats_script
    )

    def final_stats_script = file("${projectDir}/scripts/final_stats.R", checkIfExists: true)

    dada2_final_stats(
        dada2_filter_ccs.out.filter_stats.collect(),
        dada2_denoise_independent.out.denoise_stats.collect(),
        dada2_remove_chimeras.out.seqtab_nochim_rds,
        dada2_filter_asvs.out.seqtab_filtered_rds,
        metadata_ch,
        final_stats_script
    )
}

workflow download_databases {
    if (requested_dbs.contains('gtdb')) {
        download_gtdb_db(
            db_manifest.gtdb.nb.url,
            db_manifest.gtdb.nb.filename,
            db_manifest.gtdb.vsearch.seq_url
        )
    }

    if (requested_dbs.contains('silva')) {
        download_silva_db(
            db_manifest.silva.nb.url,
            db_manifest.silva.nb.filename,
            db_manifest.silva.vsearch.seq_url,
            db_manifest.silva.vsearch.tax_url
        )
    }

    if (requested_dbs.contains('gg2')) {
        download_gg2_db(
            db_manifest.gg2.nb.url,
            db_manifest.gg2.nb.filename,
            db_manifest.gg2.vsearch.seq_url,
            db_manifest.gg2.vsearch.tax_url
        )
    }
}

// ----------------------------------------------------------------------------
// Workflow dispatch
// ----------------------------------------------------------------------------
workflow {
    if (params.download_db) {
        download_databases()
    } else {
        pb16S_preprocess()
    }
}
