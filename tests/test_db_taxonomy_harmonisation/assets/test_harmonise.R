library(here)
library(tidyverse)
source(here("scripts/harmonise_taxonomy_function.R"))

data_path <- here("tests/test_db_taxonomy_harmonisation/assets")

gg2 <- read.delim(file.path(data_path, "gg2_nb.tsv"), sep = "\t", check.names = FALSE)
gg2$Taxon_harmonised <- harmonise_gg2_taxonomy(gg2$Taxon)
gg2 |>
  filter(str_detect(Taxon, "s__")) |>
  head()

gtdb <- read.delim(file.path(data_path, "gtdb_nb.tsv"), sep = "\t", check.names = FALSE)
gtdb$Taxon_harmonised <- harmonise_gtdb_taxonomy(gtdb$Taxon)

gtdb |>
  filter(str_detect(Taxon_harmonised, "s__")) 

silva <- read.delim(file.path(data_path, "silva_nb.tsv"), sep = "\t", check.names = FALSE)

silva$Taxon_harmonised <- harmonise_silva_taxonomy(silva$Taxon)
silva |>
  filter(str_detect(Taxon_harmonised, "s__")) 

