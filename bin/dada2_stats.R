#!/usr/bin/env Rscript

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 6) {
  err_quit(
    paste(
      "Expected 6 positional arguments:",
      "1) input_seqtab_filtered_rds",
      "2) metadata_tsv",
      "3) output_sample_frequency_detail_tsv",
      "4) output_dada2_qc_tsv",
      "5) output_rarefaction_depth_txt",
      "6) output_alpha_depth_txt",
      sep = "\n"
    )
  )
}

input_seqtab_rds <- args[1]
metadata_tsv <- args[2]
output_freq_tsv <- args[3]
output_qc_tsv <- args[4]
output_rarefaction_txt <- args[5]
output_alpha_txt <- args[6]

if (!file.exists(input_seqtab_rds)) {
  err_quit(sprintf("Input seqtab RDS does not exist: %s", input_seqtab_rds))
}

if (!file.exists(metadata_tsv)) {
  err_quit(sprintf("Metadata TSV does not exist: %s", metadata_tsv))
}

seqtab <- readRDS(input_seqtab_rds)

if (!is.matrix(seqtab)) {
  err_quit("Loaded seqtab is not a matrix.")
}

sample_depths <- rowSums(seqtab)

freq_df <- data.frame(
  sample = names(sample_depths),
  frequency = as.numeric(sample_depths),
  stringsAsFactors = FALSE
)

write.table(
  freq_df,
  file = output_freq_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Read metadata and merge if possible
meta <- tryCatch(
  read.table(metadata_tsv, header = TRUE, sep = "\t", quote = "", comment.char = "", check.names = FALSE),
  error = function(e) NULL
)

qc_df <- freq_df
if (!is.null(meta) && "sample_name" %in% colnames(meta)) {
  qc_df <- merge(freq_df, meta, by.x = "sample", by.y = "sample_name", all.x = TRUE, sort = FALSE)
}

write.table(
  qc_df,
  file = output_qc_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Reproduce old rarefaction heuristic
sorted_depths <- sort(sample_depths, decreasing = TRUE)
number <- length(sorted_depths)

if (number <= 2) {
  rarefaction_d <- sorted_depths[number]
  alpha_d <- rarefaction_d
} else if (number < 5) {
  result <- as.integer(number * 0.8)
  if (result < 1) result <- 1L
  rarefaction_d <- sorted_depths[result]
  alpha_d <- rarefaction_d
} else {
  result <- as.integer(number * 0.8)
  if (result < 1) result <- 1L
  rarefaction_d <- sorted_depths[result]

  alpha_res <- as.integer(number * 0.2)
  if (alpha_res < 1) alpha_res <- 1L
  alpha_d <- sorted_depths[alpha_res]
}

int_rarefaction_d <- as.integer(floor(rarefaction_d))
int_alpha_d <- as.integer(floor(alpha_d))

writeLines(as.character(int_rarefaction_d), con = output_rarefaction_txt)
writeLines(as.character(int_alpha_d), con = output_alpha_txt)

# Export for Nextflow env outputs
cat(sprintf("int_rarefaction_d=%d\n", int_rarefaction_d), file = ".command.env", append = TRUE)
cat(sprintf("int_alpha_d=%d\n", int_alpha_d), file = ".command.env", append = TRUE)
