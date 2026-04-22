process dada2_filter_ccs {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/filtered_fastq", pattern: '*.filtered.fastq.gz', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/filter_stats", pattern: '*.filter_stats.tsv', mode: params.publish_dir_mode

    cpus 1

    input:
    tuple val(sampleID), path(trimmed_fastq)
    path filter_script

    output:
    tuple val(sampleID), path("${sampleID}.filtered.fastq.gz"), emit: filtered_fastq
    path "${sampleID}.filter_stats.tsv", emit: filter_stats

    script:
    def minLenArg = params.min_len != null ? params.min_len : 'NA'
    def maxLenArg = params.max_len != null ? params.max_len : 'NA'
    def maxEeArg  = params.max_ee  != null ? params.max_ee  : 'NA'

    """
    Rscript --vanilla ${filter_script} \\
      ${trimmed_fastq} \\
      ${sampleID}.filtered.fastq.gz \\
      ${sampleID}.filter_stats.tsv \\
      ${maxEeArg} \\
      ${minLenArg} \\
      ${maxLenArg} \\
      0 \\
      FALSE
    """
}

process learn_errors {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/error_model", pattern: 'errorfun.rds', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/error_model", pattern: 'plot_error_model.pdf', mode: params.publish_dir_mode

    input:
    path filtered_fastqs
    path learn_errors_script

    output:
    path "errorfun.rds", emit: error_model
    path "plot_error_model.pdf", emit: error_plot

    script:
    def learnNbasesArg = params.learn_nbases != null ? params.learn_nbases : '1e6'

    """
    Rscript --vanilla ${learn_errors_script} \\
      ${learnNbasesArg} \\
      ${filtered_fastqs.join(' ')}
    """
}

process dada2_denoise_independent {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/denoised", pattern: '*.dada.rds', mode: params.publish_dir_mode


    input:
    tuple val(sampleID), path(filtered_fastq)
    path error_model
    path denoise_script

    output:
    tuple val(sampleID), path("${sampleID}.dada.rds"), emit: dada_rds
    path "${sampleID}.denoise_stats.tsv", emit: denoise_stats

    script:
    def omegaCArg = params.omegac != null ? params.omegac : '1e-40'
    def bandSizeArg = params.band_size != null ? params.band_size : '16'
    def hpGapArg = params.homopolymer_gap_penalty != null ? params.homopolymer_gap_penalty : 'NA'

    """
    Rscript --vanilla ${denoise_script} \\
      ${filtered_fastq} \\
      ${error_model} \\
      ${sampleID}.dada.rds \\
      ${omegaCArg} \\
      ${bandSizeArg} \\
      ${hpGapArg}
    """
}

process dada2_make_seqtab {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/seqtab", pattern: 'seqtab.rds', mode: params.publish_dir_mode

    input:
    path dada_rds_files
    path make_seqtab_script

    output:
    path "seqtab.rds", emit: seqtab_rds

    script:
    """
    Rscript --vanilla ${make_seqtab_script} \\
      ${dada_rds_files.join(' ')}
    """
}

process dada2_remove_chimeras {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/chimera_removal", pattern: 'seqtab_nochim.rds', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/chimera_removal", pattern: 'chimera_stats.tsv', mode: params.publish_dir_mode

    input:
    path seqtab_rds
    path remove_chimeras_script

    output:
    path "seqtab_nochim.rds", emit: seqtab_nochim_rds
    path "chimera_stats.tsv", emit: chimera_stats

    script:
    def chimeraMethodArg = params.chimera_method != null ? params.chimera_method : 'consensus'
    def minParentFoldArg = params.min_parent_fold != null ? params.min_parent_fold : '1.0'

    """
    Rscript --vanilla ${remove_chimeras_script} \\
      ${seqtab_rds} \\
      seqtab_nochim.rds \\
      chimera_stats.tsv \\
      ${chimeraMethodArg} \\
      ${minParentFoldArg}
    """
}

process dada2_filter_asvs {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/filtered_asvs", pattern: 'seqtab_nochim_filtered.rds', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/filtered_asvs", pattern: 'dada2_ASV.fasta', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/filtered_asvs", pattern: 'dada2_table_filtered.tsv', mode: params.publish_dir_mode
    publishDir "${params.outdir}/dada2/filtered_asvs", pattern: 'filter_asv_stats.tsv', mode: params.publish_dir_mode

    input:
    path seqtab_nochim_rds
    path filter_asvs_script
    val min_asv_totalfreq
    val min_asv_sample

    output:
    path "seqtab_nochim_filtered.rds", emit: seqtab_filtered_rds
    path "dada2_ASV.fasta", emit: asv_fasta
    path "dada2_table_filtered.tsv", emit: asv_table_tsv
    path "filter_asv_stats.tsv", emit: filter_asv_stats

    script:
    """
    Rscript --vanilla ${filter_asvs_script} \\
      ${seqtab_nochim_rds} \\
      seqtab_nochim_filtered.rds \\
      dada2_table_filtered.tsv \\
      dada2_ASV.fasta \\
      filter_asv_stats.tsv \\
      ${min_asv_totalfreq} \\
      ${min_asv_sample}
    """
}
process dada2_stats {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/dada2/", mode: params.publish_dir_mode

    input:
    path seqtab_filtered_rds
    path metadata
    path dada2_stats_script

    output:
    path "sample_frequency_detail.tsv", emit: sample_frequency_detail
    path "dada2_qc.tsv", emit: dada2_qc_tsv
    path "rarefaction_depth_suggested.txt", emit: rarefaction_depth_file
    path "alpha_depth_suggested.txt", emit: alpha_depth_file

    script:
    """
    Rscript --vanilla ${dada2_stats_script} \\
      ${seqtab_filtered_rds} \\
      ${metadata} \\
      sample_frequency_detail.tsv \\
      dada2_qc.tsv \\
      rarefaction_depth_suggested.txt \\
      alpha_depth_suggested.txt
    """
}

process dada2_final_stats {
    conda (params.enable_conda ? "$projectDir/env/qiime2-amplicon-2024.10-py310-ubuntu-conda.yml" : null)
    container "quay.io/biocontainers/bioconductor-dada2:1.38.0--r45ha27e39d_0"

    publishDir "${params.outdir}/results", mode: params.publish_dir_mode

    cpus 1

    input:
    path filter_stats_files
    path denoise_stats_files
    path seqtab_nochim_rds
    path seqtab_filtered_rds
    path metadata
    path final_stats_script

    output:
    path "dada2_tracking.tsv", emit: dada2_tracking_tsv
    path "sample_frequency_detail.tsv", emit: sample_frequency_detail_tsv
    path "dada2_qc.tsv", emit: dada2_qc_tsv
    path "rarefaction_depth_suggested.txt", emit: rarefaction_depth_file
    path "alpha_depth_suggested.txt", emit: alpha_depth_file

    script:
    """
    Rscript --vanilla ${final_stats_script} \\
      ${seqtab_nochim_rds} \\
      ${seqtab_filtered_rds} \\
      ${metadata} \\
      dada2_tracking.tsv \\
      sample_frequency_detail.tsv \\
      dada2_qc.tsv \\
      rarefaction_depth_suggested.txt \\
      alpha_depth_suggested.txt \\
      ${filter_stats_files.join(' ')} \\
      --DENOISE_STATS-- \\
      ${denoise_stats_files.join(' ')}
    """
}
