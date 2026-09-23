#!/usr/bin/env Rscript

# Generate an ASV filtering benchmark report.
#
# Usage:
#   Rscript report_asv_filters.R <benchmark_dir> <report_dir>
#
# Requires:
#   dataset_summary.tsv
#   filter_benchmark.tsv
#   asv_metrics.tsv
#   sample_depths.tsv

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  stop("Usage: Rscript report_asv_filters.R <benchmark_dir> <report_dir>")
}

input_dir <- args[1]
outdir <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Read benchmark results
# =============================================================================

read_benchmark <- function(filename) {
  path <- file.path(input_dir, filename)

  if (!file.exists(path)) {
    stop("Missing benchmark file: ", path)
  }

  read.delim(path, check.names = FALSE)
}

summary_df <- read_benchmark("dataset_summary.tsv")
benchmark <- read_benchmark("filter_benchmark.tsv")
asvs <- read_benchmark("asv_metrics.tsv")
samples <- read_benchmark("sample_depths.tsv")

n_samples <- nrow(samples)
n_asvs <- nrow(asvs)
total_reads <- sum(samples$total_reads)

# =============================================================================
# 2. Helper functions
# =============================================================================

fmt <- function(x) {
  format(
    round(x),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}

pct <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.1f%%", x))
}

metric <- function(name, column) {
  x <- summary_df[summary_df$metric == name, column]

  if (length(x) != 1) {
    stop("Missing or duplicate metric: ", name)
  }

  x
}

save_plot <- function(filename, expression) {
  png(
    file.path(outdir, filename),
    width = 1100,
    height = 750,
    res = 130
  )

  tryCatch(
    eval.parent(substitute(expression)),
    finally = dev.off()
  )
}

# =============================================================================
# 3. Describe the dataset
# =============================================================================

single_sample_asvs <- sum(asvs$sample_count == 1)
low_abundance_asvs <- sum(asvs$total_reads < 1000)
low_depth_samples <- sum(samples$total_reads < 1000)

# =============================================================================
# 4. Plot ASV abundance
# =============================================================================

save_plot("asv_abundance.png", {
  hist(
    log10(asvs$total_reads + 1),
    breaks = 60,
    main = "Distribution of ASV abundance",
    xlab = "log10(total reads per ASV + 1)",
    ylab = "Number of ASVs"
  )

  abline(
    v = log10(1000 + 1),
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = "1,000 reads",
    lty = 2,
    lwd = 2,
    bty = "n"
  )
})

# =============================================================================
# 5. Plot ASV prevalence
# =============================================================================

save_plot("asv_prevalence.png", {
  hist(
    log10(asvs$sample_count),
    breaks = 40,
    main = "Distribution of ASV prevalence",
    xlab = "log10(number of samples containing an ASV)",
    ylab = "Number of ASVs"
  )

  abline(
    v = log10(ceiling(0.1 * n_samples)),
    lty = 2,
    lwd = 2
  )

  legend(
    "topright",
    legend = "10% prevalence",
    lty = 2,
    lwd = 2,
    bty = "n"
  )
})

# =============================================================================
# 6. Plot filtering outcomes
# =============================================================================

# Use only scenarios without an additional absolute sample threshold.
comparison <- subset(
  benchmark,
  min_asv_sample == 0
)

save_plot("filter_retention.png", {
  freq_levels <- sort(unique(comparison$min_asv_total_freq))

  point_symbols <- match(
    comparison$min_asv_total_freq,
    freq_levels
  )

  plot(
    comparison$pct_asvs_retained,
    comparison$pct_reads_retained,
    pch = point_symbols,
    xlim = c(0, 100),
    ylim = c(0, 100),
    xlab = "ASVs retained (%)",
    ylab = "Reads retained (%)",
    main = "Effects of abundance and prevalence filtering"
  )

  grid()

  legend(
    "bottomright",
    legend = paste0("Min. reads: ", freq_levels),
    pch = seq_along(freq_levels),
    bty = "n",
    cex = 0.85
  )
})

# =============================================================================
# 7. Select informative filtering scenarios
# =============================================================================

scenarios <- data.frame(
  min_asv_total_freq = c(
    0, 100, 1000, 0, 0, 100, 1000
  ),
  min_asv_prevalence = c(
    0, 0, 0, 0.01, 0.05, 0.01, 0.10
  )
)

select_scenario <- function(total_freq, prevalence) {
  matches <- which(
    benchmark$min_asv_total_freq == total_freq &
      benchmark$min_asv_sample == 0 &
      abs(benchmark$min_asv_prevalence - prevalence) < 1e-9
  )

  if (length(matches) != 1) {
    stop(
      "Missing or duplicate benchmark scenario: ",
      total_freq, " / ", prevalence
    )
  }

  benchmark[matches, , drop = FALSE]
}

selected <- do.call(
  rbind,
  lapply(seq_len(nrow(scenarios)), function(i) {
    select_scenario(
      scenarios$min_asv_total_freq[i],
      scenarios$min_asv_prevalence[i]
    )
  })
)

# =============================================================================
# 8. Build the Markdown report
# =============================================================================

report <- c(
  "# ITS ASV filtering benchmark",
  "",
  "## 1. Dataset overview",
  "",
  paste0(
    "The unfiltered dataset contains **",
    fmt(n_samples), " samples**, **",
    fmt(n_asvs), " ASVs**, and **",
    fmt(total_reads), " reads**."
  ),
  "",
  "| Metric | Value |",
  "|---|---:|",
  paste0("| Samples | ", fmt(n_samples), " |"),
  paste0("| ASVs | ", fmt(n_asvs), " |"),
  paste0("| Total reads | ", fmt(total_reads), " |"),
  paste0(
    "| Median reads per sample | ",
    fmt(metric("sample_read_depth", "median")), " |"
  ),
  paste0(
    "| Median reads per ASV | ",
    fmt(metric("asv_total_reads", "median")), " |"
  ),
  paste0(
    "| Median samples per ASV | ",
    fmt(metric("asv_sample_count", "median")), " |"
  ),
  "",
  "## 2. ASV abundance",
  "",
  paste0(
    "The median ASV contains ",
    fmt(metric("asv_total_reads", "median")),
    " reads. The 90th percentile is ",
    fmt(metric("asv_total_reads", "p90")),
    " reads."
  ),
  "",
  paste0(
    "**", fmt(low_abundance_asvs), " ASVs (",
    pct(100 * low_abundance_asvs / n_asvs),
    ") have fewer than 1,000 total reads.**"
  ),
  "",
  "![ASV abundance distribution](asv_abundance.png)",
  "",
  "## 3. ASV prevalence",
  "",
  paste0(
    "**", fmt(single_sample_asvs), " ASVs (",
    pct(100 * single_sample_asvs / n_asvs),
    ") occur in only one sample.**"
  ),
  "",
  paste0(
    "A prevalence threshold of 0.1 requires presence in at least ",
    ceiling(0.1 * n_samples),
    " of the ", n_samples, " samples."
  ),
  "",
  "![ASV prevalence distribution](asv_prevalence.png)",
  "",
  "## 4. Filtering benchmark",
  "",
  "All scenarios were evaluated against the same unfiltered ",
  "post-chimera-removal sequence table.",
  "",
  "| Min. reads | Min. prevalence | ASVs retained | ASVs retained (%) | Reads retained (%) | Samples retained |",
  "|---:|---:|---:|---:|---:|---:|"
)

for (i in seq_len(nrow(selected))) {
  x <- selected[i, ]

  report <- c(
    report,
    paste0(
      "| ", fmt(x$min_asv_total_freq),
      " | ", x$min_asv_prevalence,
      " | ", fmt(x$output_asvs),
      " | ", pct(x$pct_asvs_retained),
      " | ", pct(x$pct_reads_retained),
      " | ", fmt(x$output_samples),
      " |"
    )
  )
}

report <- c(
  report,
  "",
  "### ASV retention versus read retention",
  "",
  "Each point represents a filtering scenario. ",
  "The plot compares the fraction of ASVs retained ",
  "with the fraction of reads retained.",
  "",
  "![Filtering outcomes](filter_retention.png)",
  "",
  "## 5. Sample-level considerations",
  "",
  paste0(
    "Sample read depth ranges from ",
    fmt(min(samples$total_reads)), " to ",
    fmt(max(samples$total_reads)), " reads."
  ),
  "",
  paste0(
    fmt(low_depth_samples),
    " samples contain fewer than 1,000 reads before ASV filtering."
  ),
  "",
  "Samples with very low read depth should be reviewed separately ",
  "before downstream diversity or differential abundance analyses.",
  "",
  "## 6. Interpretation",
  "",
  "- Total-frequency filtering removes ASVs with low abundance across the dataset.",
  "- Prevalence filtering removes ASVs observed in few samples.",
  "- The two filters capture different properties and can remove different ASVs.",
  "- Read retention and ASV retention should be considered together.",
  "- Prevalence is calculated across all original samples, not separately within biological groups.",
  "- Rare or group-specific ASVs are not necessarily sequencing artifacts.",
  "",
  "The benchmark quantifies filtering effects but does not establish ",
  "which threshold provides the best biological or statistical results.",
  "",
  "## 7. Reproducibility",
  "",
  paste0("Input: `", input_dir, "`"),
  "",
  "All figures and statistics were generated from the original ",
  "benchmark TSV files using base R."
)

# =============================================================================
# 9. Write report
# =============================================================================

report_file <- file.path(outdir, "ASV_filtering_report.md")

writeLines(report, report_file)

message("Report written to: ", report_file)
