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
    parseRequestedDbs
    validateDownloadParams
    validatePreprocessParams
    buildRunLog
} from './modules/validation'

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

include { 
  taxonomy_nb_assign
  taxonomy_best
  merge_taxonomy_with_table
  add_md5_to_taxonomy_table 
  } from './modules/taxonomy'

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
// Parameter parsing and Validation
// ----------------------------------------------------------------------------

def requested_dbs = parseRequestedDbs(params.download_targets)

if (params.download_db) {
    validateDownloadParams(params, requested_dbs, db_manifest)
} else {
    def n_sample = validatePreprocessParams(params)
    log.info(buildRunLog(params, version, n_sample))
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
        params.min_asv_total_freq,
        params.min_asv_sample
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
    
    def taxonomy_nb_script = file("${projectDir}/scripts/taxonomy_nb_assign.R", checkIfExists: true)
    
    def nb_db_defs = ['silva', 'gtdb', 'gg2'].collect { db ->
        def db_dir = params["${db}_dir"]
        def db_filename = db_manifest[db].nb.filename
        tuple(
            db,
            file("${db_dir}/nb/${db_filename}", checkIfExists: true)
        )
    }
    
    nb_db_ch = Channel.fromList(nb_db_defs)
    
    asv_for_nb_ch = dada2_filter_asvs.out.asv_fasta.map { fasta ->
        tuple(fasta)
    }
    
    nb_inputs_ch = asv_for_nb_ch
        .combine(nb_db_ch)
        .map { asv_fasta, db_name, db_fasta ->
            tuple(asv_fasta, db_name, db_fasta, taxonomy_nb_script)
        }
    
    taxonomy_nb_assign(nb_inputs_ch)
    
    best_script = file("${projectDir}/scripts/taxonomy_best.R", checkIfExists: true)

    taxonomy_best(
        taxonomy_nb_assign.out.nb_tax.map { it[1] }.collect(),
        params.db_to_prioritize,
        best_script
    )

    merge_script = file("${projectDir}/scripts/merge_taxonomy_with_table.R", checkIfExists: true)

    merge_taxonomy_with_table(
        taxonomy_best.out.best_tax,
        dada2_filter_asvs.out.asv_table_tsv,
        merge_script
    )
    add_md5_to_taxonomy_table(
      merge_taxonomy_with_table.out.merged_no_id
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
