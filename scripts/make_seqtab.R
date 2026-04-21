#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
})

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  err_quit("Expected at least one input .dada.rds file.")
}

dada_files <- args

missing_files <- dada_files[!file.exists(dada_files)]
if (length(missing_files) > 0) {
  err_quit(
    paste(
      "Missing input .dada.rds files:",
      paste(missing_files, collapse = "\n"),
      sep = "\n"
    )
  )
}

sample_ids <- basename(dada_files)
sample_ids <- sub("\\.dada\\.rds$", "", sample_ids)

if (anyDuplicated(sample_ids)) {
  err_quit("Duplicate sample IDs detected after stripping '.dada.rds' suffix.")
}

dds <- lapply(dada_files, readRDS)
names(dds) <- sample_ids

message("Building sequence table from ", length(dds), " denoised samples")
message("Samples:")
for (sid in sample_ids) {
  message("  - ", sid)
}

seqtab <- makeSequenceTable(dds)

saveRDS(seqtab, "seqtab.rds")

message("Wrote seqtab.rds")
quit(save = "no", status = 0)
