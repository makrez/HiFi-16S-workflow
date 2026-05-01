#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: taxonomy_best.R <nb_files...> <priority_db>")
}

priority_db <- tail(args, 1)
nb_files <- head(args, -1)

# Load all NB results
db_list <- lapply(nb_files, function(f) {
  db_name <- sub("_nb.tsv$", "", basename(f))
  dt <- fread(f)
  dt[, DB := db_name]
  dt
})

all_tax <- rbindlist(db_list)

# Function to count taxonomy depth
tax_depth <- function(t) {
  if (is.na(t) || t == "Unclassified") {
    return(0)
  }
  length(strsplit(t, ";")[[1]])
}

all_tax[, depth := sapply(Taxon, tax_depth)]

# Choose best taxonomy per ASV
best <- all_tax[,
  {
    max_depth <- max(depth, na.rm = TRUE)
    candidates <- .SD[depth == max_depth]

    if (nrow(candidates) > 1) {
      # prefer priority DB
      pref <- candidates[DB == priority_db]
      if (nrow(pref) > 0) {
        candidates <- pref
      }
    }

    candidates[1]
  },
  by = FeatureID
]

# Write outputs
# Write outputs
fwrite(
  best[, .(FeatureID, Taxon, Confidence)],
  "best_taxonomy.tsv",
  sep = "\t"
)

fwrite(
  best[, .(FeatureID, Taxon, Confidence, DB)],
  "best_taxonomy_withDB.tsv",
  sep = "\t"
)
