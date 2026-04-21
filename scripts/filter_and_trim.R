#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

parse_optional_numeric <- function(x, arg_name) {
  if (is.null(x) || x %in% c("NA", "NULL", "", "null")) {
    return(NULL)
  }
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) {
    err_quit(sprintf("Argument '%s' must be numeric, NA, or NULL. Got: %s", arg_name, x))
  }
  value
}

parse_required_numeric <- function(x, arg_name) {
  if (is.null(x) || x %in% c("NA", "NULL", "", "null")) {
    err_quit(sprintf("Missing required numeric argument: %s", arg_name))
  }
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) {
    err_quit(sprintf("Argument '%s' must be numeric. Got: %s", arg_name, x))
  }
  value
}

parse_optional_integer <- function(x, arg_name, default = NULL) {
  if (is.null(x) || x %in% c("NA", "NULL", "", "null")) {
    return(default)
  }
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) {
    err_quit(sprintf("Argument '%s' must be an integer. Got: %s", arg_name, x))
  }
  value
}

parse_optional_logical <- function(x, arg_name, default = FALSE) {
  if (is.null(x) || x %in% c("NA", "NULL", "", "null")) {
    return(default)
  }

  x_lower <- tolower(x)
  if (x_lower %in% c("true", "t", "1")) {
    return(TRUE)
  }
  if (x_lower %in% c("false", "f", "0")) {
    return(FALSE)
  }

  err_quit(sprintf("Argument '%s' must be TRUE/FALSE. Got: %s", arg_name, x))
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  err_quit(
    paste(
      "Expected exactly 8 positional arguments:",
      "1) input_fastq",
      "2) output_fastq",
      "3) output_stats",
      "4) max_expected_errors",
      "5) min_length",
      "6) max_length",
      "7) min_quality",
      "8) rm_phix",
      sep = "\n"
    )
  )
}

input_fastq <- args[1]
output_fastq <- args[2]
output_stats <- args[3]
max_expected_errors <- parse_required_numeric(args[4], "max_expected_errors")
min_len <- parse_optional_numeric(args[5], "min_length")
max_len <- parse_optional_numeric(args[6], "max_length")
min_quality <- parse_optional_integer(args[7], "min_quality", default = 0L)
rm_phix <- parse_optional_logical(args[8], "rm_phix", default = FALSE)

if (!file.exists(input_fastq)) {
  err_quit(sprintf("Input FASTQ does not exist: %s", input_fastq))
}

if (dir.exists(output_fastq)) {
  err_quit(sprintf("Output FASTQ path is a directory, not a file: %s", output_fastq))
}

if (dir.exists(output_stats)) {
  err_quit(sprintf("Output stats path is a directory, not a file: %s", output_stats))
}

output_fastq_dir <- dirname(output_fastq)
output_stats_dir <- dirname(output_stats)

if (!dir.exists(output_fastq_dir)) {
  dir.create(output_fastq_dir, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(output_stats_dir)) {
  dir.create(output_stats_dir, recursive = TRUE, showWarnings = FALSE)
}

sample_id <- basename(output_fastq)
sample_id <- sub("\\.filtered\\.fastq\\.gz$", "", sample_id)

message("Running dada2::filterAndTrim() on sample: ", sample_id)
message("Input FASTQ: ", input_fastq)
message("Output FASTQ: ", output_fastq)
message("Output stats: ", output_stats)
message("maxEE: ", max_expected_errors)
message("minLen: ", ifelse(is.null(min_len), "NULL", as.character(min_len)))
message("maxLen: ", ifelse(is.null(max_len), "NULL", as.character(max_len)))
message("minQ: ", min_quality)
message("rm.phix: ", rm_phix)

filter_out <- suppressWarnings(
  filterAndTrim(
    fwd = input_fastq,
    filt = output_fastq,
    maxEE = max_expected_errors,
    rm.phix = rm_phix,
    maxLen = max_len,
    minLen = min_len,
    minQ = min_quality,
    multithread = FALSE,
    compress = TRUE,
    verbose = TRUE
  )
)

if (!file.exists(output_fastq)) {
  err_quit("Filtering finished but output FASTQ was not created.")
}

if (!is.matrix(filter_out) && !is.data.frame(filter_out)) {
  err_quit("Unexpected filterAndTrim output type.")
}

if (nrow(filter_out) != 1) {
  err_quit(sprintf("Unexpected filterAndTrim output shape: expected 1 row, got %s.", nrow(filter_out)))
}

required_cols <- c("reads.in", "reads.out")
missing_cols <- setdiff(required_cols, colnames(filter_out))
if (length(missing_cols) > 0) {
  err_quit(sprintf(
    "filterAndTrim output is missing expected columns: %s",
    paste(missing_cols, collapse = ", ")
  ))
}

stats_df <- data.frame(
  sample = sample_id,
  input_reads = filter_out[1, "reads.in"],
  filtered_reads = filter_out[1, "reads.out"],
  stringsAsFactors = FALSE
)

write.table(
  stats_df,
  file = output_stats,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

message("Done.")
quit(save = "no", status = 0)
