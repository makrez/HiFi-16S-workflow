library(testthat)

source(test_path("..", "..", "bin", "lib", "taxonomy_harmonisation.R"))

test_that("positional SILVA/GTDB taxonomies receive double-underscore ranks", {
  expect_identical(
    harmonise_silva_taxonomy("Bacteria;Actinomycetota;Actinomycetes"),
    "d__Bacteria;p__Actinomycetota;c__Actinomycetes"
  )
  expect_identical(
    harmonise_gtdb_taxonomy("Bacteria;Bacillota;Clostridia"),
    "d__Bacteria;p__Bacillota;c__Clostridia"
  )
})

test_that("GG2-prefixed lineages do not acquire extra rank codes", {
  x <- paste(
    "d__Bacteria", "p__Actinomycetota", "c__Actinomycetes",
    "o__Propionibacteriales", "f__Propionibacteriaceae",
    "g__Cutibacterium", "s__Cutibacterium_acnes",
    sep = ";"
  )
  expect_identical(harmonise_gg2_taxonomy(x), x)
  expect_identical(normalise_taxonomy(x), x)
})

test_that("missing phylum never shifts a class assignment", {
  expected <- "d__Bacteria;p__;c__Clostridia"
  expect_identical(normalise_taxonomy("Bacteria;;Clostridia"), expected)
  expect_identical(normalise_taxonomy("d__Bacteria;p__;c__Clostridia"), expected)
  expect_identical(normalise_taxonomy("d__Bacteria;c__Clostridia"), expected)
  expect_identical(harmonise_gg2_taxonomy("d__Bacteria;c__Clostridia"), expected)
})

test_that("missing arbitrary internal ranks are padded through the deepest assignment", {
  expect_identical(
    normalise_taxonomy("d__Bacteria;g__Cutibacterium"),
    "d__Bacteria;p__;c__;o__;f__;g__Cutibacterium"
  )
  expect_identical(
    normalise_taxonomy("Bacteria;;;;;Cutibacterium;acnes"),
    "d__Bacteria;p__;c__;o__;f__;g__Cutibacterium;s__acnes"
  )
  expect_identical(
    harmonise_silva_taxonomy("Bacteria;;;;;Cutibacterium;acnes"),
    "d__Bacteria;p__;c__;o__;f__;g__Cutibacterium;s__Cutibacterium_acnes"
  )
  expect_identical(normalise_taxonomy(";Actinomycetota"), "d__;p__Actinomycetota")
})

test_that("trailing unassigned levels are omitted without shifting", {
  expect_identical(normalise_taxonomy("d__Bacteria;p__"), "d__Bacteria")
  expect_identical(normalise_taxonomy("Bacteria;Bacillota;"), "d__Bacteria;p__Bacillota")
  expect_identical(normalise_taxonomy("d__Bacteria;p__;c__"), "d__Bacteria")
})

test_that("NA and unassigned values are handled", {
  expect_identical(
    normalise_taxonomy(c(NA_character_, "", " ", "NA", "Unassigned", "d__")),
    rep(NA_character_, 6L)
  )
  expect_identical(
    normalise_taxonomy("d__Bacteria;p__Unassigned;c__Clostridia"),
    "d__Bacteria;p__;c__Clostridia"
  )
})

test_that("species normalisation retains the previous per-database behaviour", {
  expect_identical(
    harmonise_silva_taxonomy("Bacteria;;;;;Cutibacterium;acnes"),
    "d__Bacteria;p__;c__;o__;f__;g__Cutibacterium;s__Cutibacterium_acnes"
  )
  expect_identical(
    harmonise_gg2_taxonomy("d__Bacteria;g__Cutibacterium;s__sp123"),
    "d__Bacteria;p__;c__;o__;f__;g__Cutibacterium;s__Cutibacterium_sp123"
  )
  expect_identical(
    normalise_taxonomy("d__Fungi;g__Candida;s__albicans"),
    "d__Fungi;p__;c__;o__;f__;g__Candida;s__albicans"
  )
  expect_identical(
    harmonise_gg2_taxonomy("d__Bacteria;s__Cutibacterium_acnes"),
    "d__Bacteria;p__;c__;o__;f__;g__;s__Cutibacterium_acnes"
  )
})

test_that("GTDB genome suffix is removed only from species", {
  expect_identical(
    harmonise_gtdb_taxonomy(
      "Bacteria;Bacillota;Clostridia;Eubacteriales;Lachnospiraceae;Paraclostridium;Paraclostridium_bifermentans(RS_GCF_019916025_1"
    ),
    "d__Bacteria;p__Bacillota;c__Clostridia;o__Eubacteriales;f__Lachnospiraceae;g__Paraclostridium;s__Paraclostridium_bifermentans"
  )
  expect_identical(
    harmonise_gtdb_taxonomy("d__Bacteria;g__Foo;s__Foo_bar(RS_GCF_123"),
    "d__Bacteria;p__;c__;o__;f__;g__Foo;s__Foo_bar"
  )
})

test_that("normalisation is idempotent and vectorised", {
  input <- c(
    "Bacteria;Bacillota;Bacilli", "d__Bacteria;c__Clostridia",
    "d__Bacteria;p__;c__Actinomycetes", NA_character_
  )
  once <- normalise_taxonomy(input)
  expect_identical(normalise_taxonomy(once), once)
  expect_length(once, length(input))
})

test_that("malformed, mixed, or old colon input is rejected", {
  expect_error(normalise_taxonomy("d__Bacteria;p__A;p__B"), "Duplicate")
  expect_error(normalise_taxonomy("p__A;d__Bacteria"), "order")
  expect_error(normalise_taxonomy("d__Bacteria;k__A"), "Unsupported")
  expect_error(normalise_taxonomy("Bacteria;p__Bacillota"), "Mixed")
  expect_error(normalise_taxonomy("d:Bacteria;p:Bacillota"), "colon")
  expect_error(normalise_taxonomy("a;b;c;d;e;f;g;h"), "more than seven")
})

test_that("only padded double-underscore output is canonical", {
  input <- c(
    "d__Bacteria;p__Bacillota",
    "d__Bacteria;p__;c__Clostridia",
    "d__Bacteria;c__Clostridia",
    "d__Bacteria;p__",
    "d:d__Bacteria",
    "Bacteria;Bacillota",
    NA_character_, "Unassigned"
  )
  expect_identical(
    is_canonical_taxonomy(input),
    c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
  )
})
