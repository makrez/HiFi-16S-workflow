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

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  err_quit(
    paste(
      "Expected 5 positional arguments:",
      "1) input_seqtab_rds",
      "2) output_seqtab_nochim_rds",
      "3) output_chimera_stats_tsv",
      "4) chimera_method",
      "5) min_parent_fold",
      sep = "\n"
    )
  )
}

input_seqtab_rds <- args[1]
output_seqtab_nochim_rds <- args[2]
output_chimera_stats_tsv <- args[3]
chimera_method <- args[4]
min_parent_fold <- parse_required_numeric(args[5], "min_parent_fold")

if (!file.exists(input_seqtab_rds)) {
  err_quit(sprintf("Input seqtab RDS does not exist: %s", input_seqtab_rds))
}

valid_methods <- c("none", "pooled", "consensus")
if (!(chimera_method %in% valid_methods)) {
  err_quit(sprintf(
    "chimera_method must be one of: %s. Got: %s",
    paste(valid_methods, collapse = ", "),
    chimera_method
  ))
}

seqtab <- readRDS(input_seqtab_rds)

if (!is.matrix(seqtab)) {
  err_quit("Loaded seqtab is not a matrix.")
}

message("Removing chimeras")
message("Method: ", chimera_method)
message("min_parent_fold: ", min_parent_fold)
message("Input seqtab dimensions: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs")

if (chimera_method == "none") {
  seqtab_nochim <- seqtab
} else {
  seqtab_nochim <- removeBimeraDenovo(
    seqtab,
    method = chimera_method,
    minFoldParentOverAbundance = min_parent_fold,
    multithread = FALSE
  )
}

saveRDS(seqtab_nochim, output_seqtab_nochim_rds)

input_total <- sum(seqtab)
nochim_total <- sum(seqtab_nochim)

stats_df <- data.frame(
  input_total_reads = input_total,
  nonchim_total_reads = nochim_total,
  removed_reads = input_total - nochim_total,
  fraction_nonchim = if (input_total > 0) nochim_total / input_total else NA_real_,
  input_asvs = ncol(seqtab),
  nonchim_asvs = ncol(seqtab_nochim),
  removed_asvs = ncol(seqtab) - ncol(seqtab_nochim),
  stringsAsFactors = FALSE
)

write.table(
  stats_df,
  file = output_chimera_stats_tsv,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

message("Output seqtab dimensions: ", nrow(seqtab_nochim), " samples x ", ncol(seqtab_nochim), " ASVs")
message("Wrote ", output_seqtab_nochim_rds)
message("Wrote ", output_chimera_stats_tsv)

quit(save = "no", status = 0)
