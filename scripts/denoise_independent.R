#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

parse_required_numeric <- function(x, arg_name) {
  value <- suppressWarnings(as.numeric(x))
  if (is.na(value)) {
    err_quit(sprintf("Argument '%s' must be numeric. Got: %s", arg_name, x))
  }
  value
}

parse_optional_integer <- function(x, arg_name) {
  if (x %in% c("NA", "NULL", "", "null")) {
    return(NULL)
  }
  value <- suppressWarnings(as.integer(x))
  if (is.na(value)) {
    err_quit(sprintf("Argument '%s' must be an integer or NA. Got: %s", arg_name, x))
  }
  value
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  err_quit(
    paste(
      "Expected 6 positional arguments:",
      "1) filtered_fastq",
      "2) error_model_rds",
      "3) output_dada_rds",
      "4) omega_c",
      "5) band_size",
      "6) homopolymer_gap_penalty",
      sep = "\n"
    )
  )
}

filtered_fastq <- args[1]
error_model_rds <- args[2]
output_dada_rds <- args[3]
omega_c <- parse_required_numeric(args[4], "omega_c")
band_size <- parse_optional_integer(args[5], "band_size")
homopolymer_gap_penalty <- parse_optional_integer(args[6], "homopolymer_gap_penalty")

if (!file.exists(filtered_fastq)) {
  err_quit(sprintf("Filtered FASTQ does not exist: %s", filtered_fastq))
}

if (!file.exists(error_model_rds)) {
  err_quit(sprintf("Error model RDS does not exist: %s", error_model_rds))
}

err <- readRDS(error_model_rds)
drp <- derepFastq(filtered_fastq)

dada_args <- list(
  derep = drp,
  err = err,
  multithread = TRUE,
  verbose = TRUE,
  OMEGA_C = omega_c
)

if (!is.null(band_size)) {
  dada_args$BAND_SIZE <- band_size
}
if (!is.null(homopolymer_gap_penalty)) {
  dada_args$HOMOPOLYMER_GAP_PENALTY <- homopolymer_gap_penalty
}

dd <- do.call(dada, dada_args)

getN <- function(x) sum(getUniques(x))

denoised_reads <- getN(dd)

stats_df <- data.frame(
  sample = sub("\\.dada\\.rds$", "", basename(output_dada_rds)),
  denoised_reads = denoised_reads,
  stringsAsFactors = FALSE
)

write.table(
  stats_df,
  file = sub("\\.dada\\.rds$", ".denoise_stats.tsv", output_dada_rds),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

saveRDS(dd, output_dada_rds)
quit(save = "no", status = 0)
