#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(Biostrings)
})

err_quit <- function(msg, status = 1) {
  message("Error: ", msg)
  quit(save = "no", status = status)
}

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 8) {
  err_quit(
    paste(
      "Expected 8 positional arguments:",
      "1) asv_fasta",
      "2) db_fasta",
      "3) db_name",
      "4) n_threads",
      "5) min_boot",
      "6) try_rc",
      "7) output_bootstraps",
      "8) tax_levels (comma-separated or NULL)",
      sep = "\n"
    )
  )
}

# -------------------------------------------------------------------------
# Parse arguments
# -------------------------------------------------------------------------

asv_fasta <- args[1]
db_fasta <- args[2]
db_name <- args[3]

n_threads <- suppressWarnings(as.integer(args[4]))
min_boot <- suppressWarnings(as.numeric(args[5]))

parse_bool <- function(x, name) {
  value <- toupper(x)

  if (!value %in% c("TRUE", "FALSE")) {
    err_quit(sprintf("%s must be TRUE or FALSE", name))
  }

  value == "TRUE"
}

try_rc <- parse_bool(args[6], "try_rc")
output_bootstraps <- parse_bool(args[7], "output_bootstraps")

tax_levels <- if (args[8] == "NULL") {
  NULL
} else {
  trimws(strsplit(args[8], ",", fixed = TRUE)[[1]])
}

if (!file.exists(asv_fasta)) {
  err_quit(sprintf("ASV FASTA not found: %s", asv_fasta))
}

if (!file.exists(db_fasta)) {
  err_quit(sprintf("DB FASTA not found: %s", db_fasta))
}

if (is.na(n_threads) || n_threads < 1) {
  err_quit("n_threads must be a positive integer")
}

if (!is.finite(min_boot) || min_boot < 0 || min_boot > 100) {
  err_quit("min_boot must be between 0 and 100")
}

if (!is.null(tax_levels) &&
    (length(tax_levels) == 0 || any(!nzchar(tax_levels)))) {
  err_quit("tax_levels must contain non-empty rank names")
}

# -------------------------------------------------------------------------
# Read ASVs
# -------------------------------------------------------------------------

seqs <- readDNAStringSet(asv_fasta)

seq_vec <- as.character(seqs)
feature_ids <- names(seqs)

if (length(seq_vec) == 0 ||
    is.null(feature_ids) ||
    anyNA(feature_ids) ||
    any(!nzchar(feature_ids)) ||
    anyDuplicated(feature_ids)) {
  err_quit("ASV FASTA must contain unique, non-empty identifiers")
}

# -------------------------------------------------------------------------
# Taxonomy assignment
# -------------------------------------------------------------------------

message("Running assignTaxonomy for DB: ", db_name)
message("  minBoot: ", min_boot)
message("  tryRC: ", try_rc)
message("  threads: ", n_threads)

dada_args <- list(
  seqs = seq_vec,
  refFasta = db_fasta,
  minBoot = min_boot,
  tryRC = try_rc,
  multithread = n_threads,
  outputBootstraps = TRUE,
  verbose = TRUE
)

if (!is.null(tax_levels)) {
  dada_args$taxLevels <- tax_levels
}

tax <- do.call(assignTaxonomy, dada_args)

taxa_mat <- tax$tax
boot_mat <- tax$boot

if (is.null(dim(taxa_mat)) || is.null(dim(boot_mat))) {
  err_quit(sprintf(
    "Unexpected assignTaxonomy output for DB: %s",
    db_name
  ))
}

if (!identical(dim(taxa_mat), dim(boot_mat)) ||
    nrow(taxa_mat) != length(feature_ids)) {
  err_quit("Taxonomy and bootstrap dimensions do not match ASVs")
}

# -------------------------------------------------------------------------
# Collapse taxonomy
# -------------------------------------------------------------------------

collapse_taxonomy <- function(x) {
  x <- x[!is.na(x) & x != ""]

  if (length(x) == 0) {
    return("Unclassified")
  }

  paste(x, collapse = ";")
}

confidence_fun <- function(taxa, boot) {
  assigned <- which(!is.na(taxa) & taxa != "")

  if (length(assigned) == 0) {
    return(NA_real_)
  }

  deepest_rank <- max(assigned)

  suppressWarnings(as.numeric(boot[deepest_rank]))
}

confidences <- vapply(
  seq_len(nrow(taxa_mat)),
  function(i) {
    confidence_fun(
      taxa_mat[i, ],
      boot_mat[i, ]
    )
  },
  numeric(1)
)

taxon_strings <- apply(
  taxa_mat,
  1,
  collapse_taxonomy
)

# -------------------------------------------------------------------------
# Write standard taxonomy output
# -------------------------------------------------------------------------

out <- data.frame(
  FeatureID = feature_ids,
  Taxon = taxon_strings,
  Confidence = confidences,
  stringsAsFactors = FALSE
)

out_file <- sprintf("%s_nb.tsv", db_name)

write.table(
  out,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

message("Wrote: ", out_file)

# -------------------------------------------------------------------------
# Write bootstrap scores
# -------------------------------------------------------------------------

if (output_bootstraps) {

  boot_df <- data.frame(
    FeatureID = feature_ids,
    boot_mat,
    check.names = FALSE
  )

  boot_file <- sprintf(
    "%s_nb_bootstraps.tsv",
    db_name
  )

  write.table(
    boot_df,
    file = boot_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )

  message("Wrote: ", boot_file)
}

quit(save = "no", status = 0)
