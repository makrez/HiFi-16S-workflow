fix_species <- function(parts) {
  if (length(parts) < 7) return(parts)

  genus <- parts[6]
  species <- parts[7]

  if (genus == "" || species == "") return(parts)

  # Case 1: placeholder species (spXXXXX)
  if (grepl("^sp[0-9]+$", species)) {
    parts[7] <- paste0(genus, "_", species)
    return(parts)
  }

  # Case 2: already binomial
  if (startsWith(species, paste0(genus, "_"))) {
    return(parts)
  }

  # Case 3: normal species epithet
  parts[7] <- paste0(genus, "_", species)

  parts
}

pad_taxonomy <- function(parts) {
  ranks <- c("d", "p", "c", "o", "f", "g", "s")
  parts <- trimws(parts)
  parts <- parts[parts != ""]

  # Limit to max number of ranks
  if (length(parts) > length(ranks)) {
    parts <- parts[seq_along(ranks)]
  }

  # Build full taxonomy vector
  full <- rep("", length(ranks))
  full[seq_along(parts)] <- parts
  full <- fix_species(full)

  # Find last non-empty rank
  last_idx <- max(which(full != ""), 0)

  if (last_idx == 0) {
    return(NA_character_)
  }

  # Truncate to last valid rank
  full <- full[seq_len(last_idx)]
  ranks <- ranks[seq_len(last_idx)]

  paste0(ranks, "__", full, collapse = ";")
}

harmonise_gg2_taxonomy <- function(taxon) {
  vapply(taxon, function(x) {
    if (is.na(x) || x == "") {
      return(pad_taxonomy(character(0)))
    }

    parts <- strsplit(x, ";", fixed = TRUE)[[1]]
    parts <- trimws(parts)

    # GG2 already uses rank prefixes like d__, p__, ...
    parts <- sub("^[a-z]__", "", parts)

    pad_taxonomy(parts)
  }, character(1))
}

harmonise_gtdb_taxonomy <- function(taxon) {
  vapply(taxon, function(x) {
    if (is.na(x) || x == "") {
      return(pad_taxonomy(character(0)))
    }

    parts <- strsplit(x, ";", fixed = TRUE)[[1]]
    parts <- trimws(parts)

    # GTDB NB can contain species genome suffixes:
    # Paraclostridium_bifermentans(RS_GCF_019916025_1
    parts <- sub("\\(.*$", "", parts)

    pad_taxonomy(parts)
  }, character(1))
}

harmonise_silva_taxonomy <- function(taxon) {
  vapply(taxon, function(x) {
    if (is.na(x) || x == "") {
      return(pad_taxonomy(character(0)))
    }

    parts <- strsplit(x, ";", fixed = TRUE)[[1]]
    parts <- trimws(parts)

    pad_taxonomy(parts)
  }, character(1))
}
