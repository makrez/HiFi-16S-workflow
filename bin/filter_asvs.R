#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Biostrings)
})

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

parse_required_numeric <- function(x, arg_name) {
  value <- suppressWarnings(as.numeric(x))

  if (length(value) != 1 || !is.finite(value) || value < 0) {
    err_quit(sprintf(
      "Argument '%s' must be a non-negative finite number. Got: %s",
      arg_name, x
    ))
  }

  value
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  err_quit(
    paste(
      "Expected 8 positional arguments:",
      "1) input_seqtab_nochim_rds",
      "2) output_seqtab_filtered_rds",
      "3) output_table_tsv",
      "4) output_asv_fasta",
      "5) output_filter_stats_tsv",
      "6) min_asv_totalfreq",
      "7) min_asv_sample",
      "8) min_asv_prevalence (0-1)",
      sep = "\n"
    )
  )
}

# ------------------------------------------------------------------------------
# Parse arguments
# ------------------------------------------------------------------------------

input_seqtab_rds <- args[1]
output_seqtab_rds <- args[2]
output_table_tsv <- args[3]
output_asv_fasta <- args[4]
output_filter_stats_tsv <- args[5]

min_asv_totalfreq <- parse_required_numeric(
  args[6], "min_asv_totalfreq"
)

min_asv_sample <- parse_required_numeric(
  args[7], "min_asv_sample"
)

min_asv_prevalence <- parse_required_numeric(
  args[8], "min_asv_prevalence"
)

if (min_asv_sample != floor(min_asv_sample)) {
  err_quit("min_asv_sample must be an integer.")
}

if (min_asv_prevalence > 1) {
  err_quit("min_asv_prevalence must be between 0 and 1.")
}

# ------------------------------------------------------------------------------
# Read sequence table
# ------------------------------------------------------------------------------

if (!file.exists(input_seqtab_rds)) {
  err_quit(sprintf(
    "Input seqtab RDS does not exist: %s",
    input_seqtab_rds
  ))
}

seqtab <- readRDS(input_seqtab_rds)

if (!is.matrix(seqtab)) {
  err_quit("Loaded seqtab is not a matrix.")
}

if (nrow(seqtab) == 0 || ncol(seqtab) == 0) {
  err_quit("Input seqtab has zero rows or zero columns.")
}

if (
  !is.numeric(seqtab) ||
    any(!is.finite(seqtab)) ||
    any(seqtab < 0)
) {
  err_quit("Input seqtab must contain finite, non-negative counts.")
}

if (
  is.null(rownames(seqtab)) ||
    anyNA(rownames(seqtab)) ||
    any(!nzchar(rownames(seqtab))) ||
    anyDuplicated(rownames(seqtab))
) {
  err_quit(
    "Input seqtab must contain unique, non-empty sample names as row names."
  )
}

if (
  is.null(colnames(seqtab)) ||
    anyNA(colnames(seqtab)) ||
    any(!nzchar(colnames(seqtab))) ||
    anyDuplicated(colnames(seqtab))
) {
  err_quit(
    paste(
      "Input seqtab must contain unique, non-empty ASV sequences",
      "as column names."
    )
  )
}

input_n_samples <- nrow(seqtab)
input_n_asvs <- ncol(seqtab)
input_total_reads <- sum(seqtab)

# ------------------------------------------------------------------------------
# Calculate ASV abundance and prevalence
# ------------------------------------------------------------------------------

# Total number of reads per ASV across all samples
asv_total_reads <- colSums(seqtab)

# Number of samples in which each ASV is present
asv_sample_count <- colSums(seqtab > 0)

# Fraction of samples containing each ASV.
# Calculated before filtering or removing empty samples.
asv_prevalence <- asv_sample_count / input_n_samples

# ------------------------------------------------------------------------------
# Apply filtering thresholds
# ------------------------------------------------------------------------------

# Minimum total read count
keep_totalfreq <- if (min_asv_totalfreq > 0) {
  asv_total_reads >= min_asv_totalfreq
} else {
  rep(TRUE, input_n_asvs)
}

# Minimum number of samples containing the ASV
keep_samplecount <- if (min_asv_sample > 0) {
  asv_sample_count >= min_asv_sample
} else {
  rep(TRUE, input_n_asvs)
}

# Minimum fraction of samples containing the ASV
keep_prevalence <- if (min_asv_prevalence > 0) {
  asv_prevalence >= min_asv_prevalence
} else {
  rep(TRUE, input_n_asvs)
}

# ASVs must pass all enabled filters
keep_asvs <- (
  keep_totalfreq &
    keep_samplecount &
    keep_prevalence
)

if (!any(keep_asvs)) {
  err_quit("All ASVs were removed by filtering.")
}

# Keep all samples initially so statistics can also be reported for samples
# that lose all reads during ASV filtering.
seqtab_filt_all_samples <- seqtab[, keep_asvs, drop = FALSE]

# ------------------------------------------------------------------------------
# Calculate per-sample filtering statistics
# ------------------------------------------------------------------------------

reads_before <- rowSums(seqtab)
reads_after <- rowSums(seqtab_filt_all_samples)

asvs_before <- rowSums(seqtab > 0)
asvs_after <- rowSums(seqtab_filt_all_samples > 0)

stats_df <- data.frame(
  sample = rownames(seqtab),
  reads_before_asv_filtering = reads_before,
  reads_after_asv_filtering = reads_after,
  number_of_asvs_before_filtering = asvs_before,
  number_of_asvs_after_filtering = asvs_after,
  condition = "sample",
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# Remove empty samples
# ------------------------------------------------------------------------------

# Drop samples without reads after ASV filtering,
# matching the previous QIIME behavior.
keep_samples <- reads_after > 0

seqtab_filt <- seqtab_filt_all_samples[
  keep_samples, ,
  drop = FALSE
]

if (nrow(seqtab_filt) == 0) {
  err_quit("All samples were removed by filtering.")
}

# ------------------------------------------------------------------------------
# Write filtered sequence table
# ------------------------------------------------------------------------------

saveRDS(
  seqtab_filt,
  output_seqtab_rds
)

write.table(
  seqtab_filt,
  file = output_table_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = TRUE,
  col.names = NA
)

# ------------------------------------------------------------------------------
# Write representative ASV sequences
# ------------------------------------------------------------------------------

asv_seqs <- colnames(seqtab_filt)

repseqs <- DNAStringSet(asv_seqs)

names(repseqs) <- asv_seqs

writeXStringSet(
  repseqs,
  filepath = output_asv_fasta
)

# ------------------------------------------------------------------------------
# Write filtering statistics
# ------------------------------------------------------------------------------

write.table(
  stats_df,
  file = output_filter_stats_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ------------------------------------------------------------------------------
# Report summary
# ------------------------------------------------------------------------------

output_n_samples <- nrow(seqtab_filt)
output_n_asvs <- ncol(seqtab_filt)
output_total_reads <- sum(seqtab_filt)

message(
  "ASV filtering completed: ",
  input_n_asvs, " -> ", output_n_asvs, " ASVs; ",
  input_total_reads, " -> ", output_total_reads, " reads; ",
  input_n_samples, " -> ", output_n_samples, " samples."
)

message(
  "Filtering thresholds: ",
  "min_asv_totalfreq=", min_asv_totalfreq, ", ",
  "min_asv_sample=", min_asv_sample, ", ",
  "min_asv_prevalence=", min_asv_prevalence
)

message(
  "ASVs passing individual filters: ",
  "total frequency=", sum(keep_totalfreq), ", ",
  "sample count=", sum(keep_samplecount), ", ",
  "prevalence=", sum(keep_prevalence), ", ",
  "all filters=", sum(keep_asvs)
)

quit(save = "no", status = 0)
