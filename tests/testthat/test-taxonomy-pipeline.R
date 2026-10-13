library(testthat)
library(data.table)

repo_root <- normalizePath(test_path("..", ".."))

run_cli <- function(script, ...) {
  command <- file.path(R.home("bin"), "Rscript")
  arguments <- shQuote(c(file.path(repo_root, "bin", script), ...))
  output <- suppressWarnings(system2(
    command,
    args = arguments, stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(as.integer(status), 0L, info = paste(output, collapse = "\n"))
}

make_abundances <- function(file) {
  fwrite(data.table(
    sample = c("S1", "S2"),
    ACGT = c(10L, 3L),
    TGCA = c(2L, 5L),
    NNNN = c(1L, 0L)
  ), file, sep = "\t")
}

summarise_merged <- function(merged_file, workdir) {
  merged <- fread(
    merged_file,
    na.strings = c("", "NA")
  )
  merged[, id := paste0("asv", seq_len(.N))]
  setcolorder(merged, c("id", "Sequence", "Taxon", "Confidence", "S1", "S2"))
  with_ids <- file.path(workdir, "merged-with-ids.tsv")
  assignment <- file.path(workdir, "assignment.tsv")
  deepest <- file.path(workdir, "deepest.tsv")
  fwrite(merged, with_ids, sep = "\t", na = "")
  run_cli("summarize_taxonomy.R", with_ids, assignment, deepest)
  list(merged = merged, assignment = fread(assignment), deepest = fread(deepest))
}

test_that("NB harmonisation, merge and summary preserve ranks and reads", {
  workdir <- tempfile("taxonomy-nb-")
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  raw <- file.path(workdir, "raw.tsv")
  harmonised <- file.path(workdir, "harmonised.tsv")
  abundance <- file.path(workdir, "abundance.tsv")
  merged <- file.path(workdir, "merged.tsv")

  fwrite(data.table(
    FeatureID = c("ACGT", "TGCA"),
    Taxon = c("d__Bacteria;c__Actinomycetes", "d__Bacteria;p__Bacillota"),
    Confidence = c(90L, 85L)
  ), raw, sep = "\t")
  make_abundances(abundance)

  run_cli("harmonise_taxonomy.R", raw, "gg2", harmonised)
  expect_identical(
    fread(harmonised)[FeatureID == "ACGT", Taxon],
    "d__Bacteria;p__;c__Actinomycetes"
  )

  run_cli("merge_taxonomy_with_table.R", harmonised, abundance, merged)
  result <- summarise_merged(merged, workdir)

  expect_identical(
    result$merged[Sequence == "ACGT", Taxon],
    "d__Bacteria;p__;c__Actinomycetes"
  )
  expect_identical(
    result$merged[Sequence == "TGCA", Taxon],
    "d__Bacteria;p__Bacillota"
  )
  expect_true(is.na(result$merged[Sequence == "NNNN", Taxon]))
  expect_false(any(grepl("[dpcofgs]:[dpcofgs]__", result$merged$Taxon)))

  # ACGT has 13 reads (Class, missing Phylum), TGCA has 7 (Phylum).
  # NNNN has 1 read and no taxonomy.
  stats <- result$assignment
  expect_equal(stats[rank == "Domain", assigned_asvs], 2L)
  expect_equal(stats[rank == "Phylum", assigned_asvs], 1L)
  expect_equal(stats[rank == "Class", assigned_asvs], 1L)
  expect_equal(stats[rank == "Class", assigned_reads], 13)
  expect_equal(stats[rank == "Phylum", assigned_reads], 7)
  expect_equal(stats[rank == "Domain", total_reads], 21)

  deepest <- result$deepest
  expect_equal(deepest[deepest_rank == "Class", asvs], 1L)
  expect_equal(deepest[deepest_rank == "Phylum", reads], 7)
  expect_equal(deepest[deepest_rank == "Unassigned", reads], 1)
})

test_that("native SINTAX VSEARCH labels become d__ in merged output", {
  workdir <- tempfile("taxonomy-vsearch-")
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  lca <- file.path(workdir, "lca.tsv")
  hits <- file.path(workdir, "hits.tsv")
  abundance <- file.path(workdir, "abundance.tsv")
  merged <- file.path(workdir, "merged.tsv")

  writeLines(c(
    "ACGT\td:Fungi,g:Candida,s:albicans",
    "TGCA\td:Fungi,p:Ascomycota"
  ), lca)
  writeLines(c(
    "ACGT\tref1;tax=d:Fungi,g:Candida,s:albicans;\t99\t100\t4\t4\t4",
    "TGCA\tref2;tax=d:Fungi,p:Ascomycota;\t99\t100\t4\t4\t4"
  ), hits)
  make_abundances(abundance)

  run_cli("merge_vsearch_taxonomy_with_table.R", lca, hits, abundance, merged)
  result <- summarise_merged(merged, workdir)

  expect_identical(
    result$merged[Sequence == "ACGT", Taxon],
    "d__Fungi;p__;c__;o__;f__;g__Candida;s__albicans"
  )
  expect_identical(
    result$merged[Sequence == "TGCA", Taxon],
    "d__Fungi;p__Ascomycota"
  )
  # Native SINTAX species epithets are NOT renamed by this formatting fix.
  expect_equal(result$assignment[rank == "Species", assigned_reads], 13)
  expect_equal(result$deepest[deepest_rank == "Species", reads], 13)
})

test_that("summarisation rejects legacy d: or duplicated rank prefixes", {
  workdir <- tempfile("taxonomy-invalid-")
  dir.create(workdir)
  on.exit(unlink(workdir, recursive = TRUE), add = TRUE)

  infile <- file.path(workdir, "invalid.tsv")
  fwrite(data.table(
    id = "asv1", Sequence = "ACGT",
    Taxon = "d:d__Bacteria;p:p__Bacillota", Confidence = 90L,
    S1 = 10L
  ), infile, sep = "\t")

  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    args = shQuote(c(
      file.path(repo_root, "bin", "summarize_taxonomy.R"),
      infile,
      file.path(workdir, "assignment.tsv"),
      file.path(workdir, "deepest.tsv")
    )),
    stdout = TRUE, stderr = TRUE
  ))
  expect_true(!is.null(attr(output, "status")))
  expect_match(paste(output, collapse = "\n"), "Non-canonical taxonomy")
})
