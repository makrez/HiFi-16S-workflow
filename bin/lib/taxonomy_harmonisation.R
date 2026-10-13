# Canonical taxonomy format: d__Bacteria;p__...;c__...;o__...;f__...;g__...;s__...
# Ranks are positional for unprefixed input and explicit for prefixed input.
# Missing internal ranks are retained as empty x__ components through the
# deepest assigned rank. No missing rank may shift another rank's position.
TAXON_RANKS <- c("d", "p", "c", "o", "f", "g", "s")

fix_species <- function(parts) {
  if (length(parts) < 7L) {
    return(parts)
  }

  genus <- parts[6]
  species <- parts[7]

  if (genus == "" || species == "") {
    return(parts)
  }

  # Preserve existing binomials; expand epithets and sp123 placeholders.
  if (!startsWith(species, paste0(genus, "_"))) {
    parts[7] <- paste0(genus, "_", species)
  }
  parts
}

normalise_taxonomy <- function(taxon, add_genus_to_species = FALSE) {
  vapply(taxon, function(x) {
    if (is.na(x) || trimws(x) %in% c("", "NA", "Unassigned")) {
      return(NA_character_)
    }

    parts <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
    if (length(parts) > length(TAXON_RANKS)) {
      stop("Taxonomy contains more than seven ranks: ", x, call. = FALSE)
    }

    # Old merged d:Bacteria syntax must not be silently turned into d__d:Bacteria.
    if (any(grepl("^[[:alpha:]]:", parts))) {
      stop("Old colon rank prefixes are not supported: ", x, call. = FALSE)
    }

    prefixed <- grepl("^[[:alpha:]]__", parts)
    if (any(prefixed) && !all(prefixed)) {
      stop("Mixed prefixed and positional taxonomy: ", x, call. = FALSE)
    }

    values <- rep("", length(TAXON_RANKS))
    if (all(prefixed)) {
      codes <- sub("__.*$", "", parts)
      positions <- match(codes, TAXON_RANKS)

      if (anyNA(positions)) {
        stop("Unsupported taxonomy rank prefix: ", x, call. = FALSE)
      }
      if (anyDuplicated(positions)) {
        stop("Duplicate taxonomy rank: ", x, call. = FALSE)
      }
      if (is.unsorted(positions, strictly = TRUE)) {
        stop("Incorrect taxonomy rank order: ", x, call. = FALSE)
      }
      values[positions] <- sub("^[[:alpha:]]__", "", parts)
    } else {
      # Empty positional fields stay in place: Bacteria;;Clostridia means
      # d__Bacteria;p__;c__Clostridia, never p__Clostridia.
      values[seq_along(parts)] <- parts
    }

    values[values %in% c("NA", "Unassigned")] <- ""
    if (add_genus_to_species) values <- fix_species(values)

    assigned <- which(nzchar(values))
    if (!length(assigned)) {
      return(NA_character_)
    }

    last <- max(assigned)
    paste0(
      TAXON_RANKS[seq_len(last)], "__", values[seq_len(last)],
      collapse = ";"
    )
  }, character(1), USE.NAMES = FALSE)
}

# Backward-compatible helper interface for callers with positional components.
pad_taxonomy <- function(parts) {
  normalise_taxonomy(paste(parts, collapse = ";"), add_genus_to_species = TRUE)
}

harmonise_gg2_taxonomy <- function(taxon) {
  normalise_taxonomy(taxon, add_genus_to_species = TRUE)
}

harmonise_silva_taxonomy <- function(taxon) {
  normalise_taxonomy(taxon, add_genus_to_species = TRUE)
}

harmonise_gtdb_taxonomy <- function(taxon) {
  cleaned <- vapply(taxon, function(x) {
    if (is.na(x)) {
      return(NA_character_)
    }

    parts <- strsplit(x, ";", fixed = TRUE)[[1]]
    if (!length(parts)) {
      return(x)
    }

    # Only remove GTDB genome identifiers from the species field.
    has_prefix <- grepl("^[[:alpha:]]__", trimws(parts))
    species_idx <- if (all(has_prefix)) {
      which(grepl("^s__", trimws(parts)))
    } else if (!any(has_prefix) && length(parts) >= 7L) {
      7L
    } else {
      integer(0)
    }
    if (length(species_idx)) {
      parts[species_idx] <- sub("\\(.*$", "", parts[species_idx])
    }
    paste(parts, collapse = ";")
  }, character(1), USE.NAMES = FALSE)

  normalise_taxonomy(cleaned, add_genus_to_species = TRUE)
}

# Summaries require an already-canonical lineage, but accept unassigned ASVs.
is_canonical_taxonomy <- function(taxon) {
  vapply(taxon, function(x) {
    if (is.na(x) || trimws(x) %in% c("", "Unassigned")) {
      return(TRUE)
    }
    canonical <- tryCatch(
      normalise_taxonomy(x),
      error = function(e) NA_character_
    )
    !is.na(canonical) && identical(x, canonical)
  }, logical(1), USE.NAMES = FALSE)
}
