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

nextflow.enable.dsl = 2

version = "0.1.0"

include {
    inspect_metadata
    QC_fastq
    cutadapt
    summarize_cutadapt
    QC_fastq_post_trim
    collect_QC
    collect_QC_skip_cutadapt
} from './modules/qc'

// ----------------------------------------------------------------------------
// Help text
// ----------------------------------------------------------------------------
def helpMessage() {
    return """
    Minimal preprocessing workflow for PacBio HiFi 16S

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

// ----------------------------------------------------------------------------
// Basic parameter validation
// ----------------------------------------------------------------------------
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

// ----------------------------------------------------------------------------
// Dynamic settings
// ----------------------------------------------------------------------------
n_sample = file(params.input).countLines() - 1
if (n_sample < 1) {
    exit 1, "Input sample sheet appears empty: ${params.input}"
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

// ----------------------------------------------------------------------------
// Workflow
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

    if (params.skip_primer_trim) {

        collect_QC_skip_cutadapt(
            QC_fastq.out.all_seqkit_stats.collect(),
            QC_fastq.out.all_seqkit_summary.collect()
        )

        filtered_reads = QC_fastq.out.filtered_fastq

    } else {

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

        filtered_reads = cutadapt.out.cutadapt_fastq
    }

}

workflow {
    pb16S_preprocess()
}

