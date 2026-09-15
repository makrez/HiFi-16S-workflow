#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Usage: merge_taxonomy_with_table.R <taxonomy.tsv> <asv_table.tsv> <out.tsv>")
}

tax_file <- args[1]
asv_file <- args[2]
out_file <- args[3]

if (!file.exists(tax_file)) {
  stop(sprintf("Taxonomy file not found: %s", tax_file))
}

if (!file.exists(asv_file)) {
  stop(sprintf("ASV table file not found: %s", asv_file))
}

# ------------------------------------------------------------------------------
# Read taxonomy table
# Expected current schema from taxonomy_best:
#   FeatureID / id / Sequence  + Taxon + Confidence
# We normalize the join key to "Sequence"
# ------------------------------------------------------------------------------
tax <- fread(tax_file, check.names = FALSE)

if ("FeatureID" %in% colnames(tax)) {
  setnames(tax, "FeatureID", "Sequence")
} else if ("id" %in% colnames(tax)) {
  # If taxonomy currently uses hashed IDs, this merge script cannot infer Sequence.
  # For the current pipeline state, taxonomy should still use sequence-based IDs.
  stop("Taxonomy file contains 'id' but no 'FeatureID'/'Sequence'. This script expects sequence-based taxonomy keys before MD5 is added.")
}

required_tax_cols <- c("Sequence", "Taxon", "Confidence")
missing_tax_cols <- setdiff(required_tax_cols, colnames(tax))
if (length(missing_tax_cols) > 0) {
  stop(sprintf(
    "Taxonomy file is missing required columns: %s",
    paste(missing_tax_cols, collapse = ", ")
  ))
}

tax <- tax[, .(Sequence, Taxon, Confidence)]

# ------------------------------------------------------------------------------
# Normalize taxonomy to the canonical pipeline format
#
# Native NB taxonomy is positional:
#
#   Cryptista;Cryptophyta;Cryptophyceae;...;pyrenoidifera
#
# Canonical pipeline format:
#
#   d:Cryptista;p:Cryptophyta;c:Cryptophyceae;...;s:pyrenoidifera
# ------------------------------------------------------------------------------

rank_codes <- c("d", "p", "c", "o", "f", "g", "s")

canonicalize_nb_taxonomy <- function(taxon) {

  if (is.na(taxon) || trimws(taxon) == "") {
    return(NA_character_)
  }

  values <- strsplit(
    taxon,
    ";",
    fixed = TRUE
  )[[1]]

  values <- trimws(values)

  if (length(values) > length(rank_codes)) {
    stop(sprintf(
      "NB taxonomy contains %d ranks; expected at most %d: %s",
      length(values),
      length(rank_codes),
      taxon
    ))
  }

  # Assign rank by POSITION before removing empty entries.
  # This is important because an empty intermediate rank must not
  # cause all subsequent ranks to shift.
  ranked <- paste0(
    rank_codes[seq_along(values)],
    ":",
    values
  )

  keep <- (
    !is.na(values) &
      values != "" &
      values != "NA" &
      values != "Unassigned"
  )

  ranked <- ranked[keep]

  if (length(ranked) == 0) {
    return(NA_character_)
  }

  paste(
    ranked,
    collapse = ";"
  )
}

tax[
  ,
  Taxon := vapply(
    Taxon,
    canonicalize_nb_taxonomy,
    character(1)
  )
]

# ------------------------------------------------------------------------------
# Read ASV abundance table
# Expected format:
#   first column = sample name
#   remaining columns = ASV sequences
# ------------------------------------------------------------------------------
asv <- fread(asv_file, check.names = FALSE)

if (ncol(asv) < 2) {
  stop("ASV table must contain at least one sample column and one ASV column")
}

setnames(asv, 1, "sample")

sample_names <- asv[["sample"]]
seq_cols <- setdiff(colnames(asv), "sample")

if (length(seq_cols) == 0) {
  stop("No ASV sequence columns found in ASV table")
}

# Convert sample x ASV table to ASV x sample table
abund_mat <- as.matrix(asv[, ..seq_cols])
rownames(abund_mat) <- sample_names

asv_long <- as.data.table(t(abund_mat), keep.rownames = "Sequence")

# After transpose, columns are sample names
setnames(
  asv_long,
  old = setdiff(colnames(asv_long), "Sequence"),
  new = sample_names
)

# ------------------------------------------------------------------------------
# Merge taxonomy with abundance
# ------------------------------------------------------------------------------
merged <- merge(asv_long, tax, by = "Sequence", all.x = TRUE, sort = FALSE)

# Reorder to old-pipeline style minus the md5 id:
# Sequence, Taxon, Confidence, <sample columns...>
sample_cols <- setdiff(colnames(merged), c("Sequence", "Taxon", "Confidence"))
setcolorder(merged, c("Sequence", "Taxon", "Confidence", sample_cols))

fwrite(merged, out_file, sep = "\t", quote = FALSE, na = "")

message("Merged table written: ", out_file)
