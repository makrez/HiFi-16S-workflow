#!/usr/bin/env Rscript

# Benchmark abundance and prevalence filters on a DADA2
# post-chimera-removal sequence table.
#
# Usage:
#   Rscript benchmark_asv_filters.R <seqtab_nochim.rds> <output_dir>
#
# No additional R packages required.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop(
    "Usage: Rscript benchmark_asv_filters.R ",
    "<seqtab_nochim.rds> <output_dir>"
  )
}

input_file <- args[1]
outdir <- args[2]

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, filename) {
  write.table(
    x,
    file = file.path(outdir, filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = "NA"
  )
}

# =============================================================================
# 1. Load data
# =============================================================================

seqtab <- readRDS(input_file)

if (!is.matrix(seqtab)) {
  stop("Input must be a matrix.")
}

if (nrow(seqtab) == 0 || ncol(seqtab) == 0) {
  stop("Input sequence table is empty.")
}

if (!is.numeric(seqtab) ||
  any(!is.finite(seqtab)) ||
  any(seqtab < 0)) {
  stop("Sequence table must contain finite, non-negative counts.")
}

n_samples <- nrow(seqtab)
n_asvs <- ncol(seqtab)
total_reads <- sum(seqtab)

sample_ids <- rownames(seqtab)

if (is.null(sample_ids)) {
  sample_ids <- paste0("sample_", seq_len(n_samples))
}

asv_ids <- colnames(seqtab)

if (is.null(asv_ids)) {
  stop("Input sequence table has no ASV column names.")
}

# =============================================================================
# 2. Calculate abundance and prevalence
# =============================================================================

# Total reads per ASV
asv_total_reads <- colSums(seqtab)

# Number of samples containing each ASV
asv_sample_count <- colSums(seqtab > 0)

# Fraction of samples containing each ASV
asv_prevalence <- asv_sample_count / n_samples

# Total reads per sample
sample_depth <- rowSums(seqtab)

# =============================================================================
# 3. Describe the unfiltered dataset
# =============================================================================

describe_vector <- function(x, metric) {
  q <- quantile(
    x,
    probs = c(0, 0.10, 0.25, 0.50, 0.75, 0.90, 1),
    names = FALSE
  )

  data.frame(
    metric = metric,
    n = length(x),
    mean = mean(x),
    min = q[1],
    p10 = q[2],
    p25 = q[3],
    median = q[4],
    p75 = q[5],
    p90 = q[6],
    max = q[7]
  )
}

dataset_summary <- do.call(
  rbind,
  list(
    describe_vector(
      asv_total_reads,
      "asv_total_reads"
    ),
    describe_vector(
      asv_sample_count,
      "asv_sample_count"
    ),
    describe_vector(
      asv_prevalence,
      "asv_prevalence"
    ),
    describe_vector(
      sample_depth,
      "sample_read_depth"
    )
  )
)

write_tsv(
  dataset_summary,
  "dataset_summary.tsv"
)

# Individual ASV metrics
asv_metrics <- data.frame(
  FeatureID = asv_ids,
  total_reads = asv_total_reads,
  sample_count = asv_sample_count,
  prevalence = asv_prevalence,
  stringsAsFactors = FALSE
)

write_tsv(
  asv_metrics,
  "asv_metrics.tsv"
)

# Individual sample read depths
sample_depths <- data.frame(
  sample = sample_ids,
  total_reads = sample_depth,
  stringsAsFactors = FALSE
)

write_tsv(
  sample_depths,
  "sample_depths.tsv"
)

# =============================================================================
# 4. Define filtering thresholds
# =============================================================================

# Minimum total read count per ASV
total_freq_thresholds <- c(
  0, 10, 50, 100, 1000, 5000
)

# Minimum absolute number of samples
sample_thresholds <- c(
  0, 2, 5
)

# Minimum fraction of samples
prevalence_thresholds <- c(
  0, 0.01, 0.05, 0.10, 0.20
)

# Generate all combinations
grid <- expand.grid(
  min_asv_total_freq = total_freq_thresholds,
  min_asv_sample = sample_thresholds,
  min_asv_prevalence = prevalence_thresholds
)

# =============================================================================
# 5. Benchmark filtering combinations
# =============================================================================

benchmark_filter <- function(min_total, min_samples, min_prev) {
  # Total abundance
  keep_totalfreq <- if (min_total > 0) {
    asv_total_reads >= min_total
  } else {
    rep(TRUE, n_asvs)
  }

  # Absolute prevalence
  keep_samplecount <- if (min_samples > 0) {
    asv_sample_count >= min_samples
  } else {
    rep(TRUE, n_asvs)
  }

  # Relative prevalence
  keep_prevalence <- if (min_prev > 0) {
    asv_prevalence >= min_prev
  } else {
    rep(TRUE, n_asvs)
  }

  # Combine all enabled filters
  keep_asvs <- (
    keep_totalfreq &
      keep_samplecount &
      keep_prevalence
  )

  # Retained ASVs
  output_asvs <- sum(keep_asvs)

  # Retained reads per sample
  filtered_sample_depth <- rowSums(
    seqtab[, keep_asvs, drop = FALSE]
  )

  # Total retained reads
  output_reads <- sum(filtered_sample_depth)

  # Samples remaining after dropping empty samples
  output_samples <- sum(filtered_sample_depth > 0)

  # Effective minimum number of samples required
  effective_min_samples <- max(
    min_samples,
    ceiling(min_prev * n_samples)
  )

  data.frame(
    min_asv_total_freq = min_total,
    min_asv_sample = min_samples,
    min_asv_prevalence = min_prev,
    effective_min_samples = effective_min_samples,
    input_samples = n_samples,
    output_samples = output_samples,
    removed_samples = n_samples - output_samples,
    pct_samples_retained = 100 * output_samples / n_samples,
    input_asvs = n_asvs,
    output_asvs = output_asvs,
    removed_asvs = n_asvs - output_asvs,
    pct_asvs_retained = 100 * output_asvs / n_asvs,
    input_total_reads = total_reads,
    output_total_reads = output_reads,
    removed_reads = total_reads - output_reads,
    pct_reads_retained = if (total_reads > 0) {
      100 * output_reads / total_reads
    } else {
      NA_real_
    },
    asvs_passing_totalfreq = sum(keep_totalfreq),
    asvs_passing_samplecount = sum(keep_samplecount),
    asvs_passing_prevalence = sum(keep_prevalence)
  )
}

results <- lapply(seq_len(nrow(grid)), function(i) {
  benchmark_filter(
    min_total = grid$min_asv_total_freq[i],
    min_samples = grid$min_asv_sample[i],
    min_prev = grid$min_asv_prevalence[i]
  )
})

benchmark <- do.call(rbind, results)

# Sort by abundance, prevalence, and sample threshold
benchmark <- benchmark[
  order(
    benchmark$min_asv_total_freq,
    benchmark$min_asv_prevalence,
    benchmark$min_asv_sample
  ),
]

rownames(benchmark) <- NULL

write_tsv(
  benchmark,
  "filter_benchmark.tsv"
)

# =============================================================================
# 6. Print overview
# =============================================================================

cat("\n=== Input dataset ===\n")

cat("Samples:     ", n_samples, "\n")
cat("ASVs:        ", n_asvs, "\n")
cat("Total reads: ", total_reads, "\n")

cat(
  "ASVs found in only one sample: ",
  sum(asv_sample_count == 1),
  "\n"
)

cat(
  "ASVs with fewer than 100 total reads: ",
  sum(asv_total_reads < 100),
  "\n"
)

cat(
  "ASVs with fewer than 1000 total reads: ",
  sum(asv_total_reads < 1000),
  "\n"
)

cat("\n=== Benchmark ===\n")

cat("Threshold combinations: ", nrow(benchmark), "\n")

cat("\nOutput directory: ", outdir, "\n", sep = "")

cat("\nFiles written:\n")
cat("  dataset_summary.tsv\n")
cat("  asv_metrics.tsv\n")
cat("  sample_depths.tsv\n")
cat("  filter_benchmark.tsv\n")
