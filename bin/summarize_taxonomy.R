#!/usr/bin/env Rscript

suppressPackageStartupMessages({
    library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
    stop(
        paste(
            "Usage: summarize_taxonomy.R",
            "<taxonomy_table.tsv>",
            "<assignment_summary.tsv>",
            "<deepest_rank_summary.tsv>"
        )
    )
}

input_file <- args[1]
assignment_out <- args[2]
deepest_out <- args[3]

# -------------------------------------------------------------------------
# Taxonomic ranks
#
# Canonical taxonomy format expected by this script:
#
#   d:Bacteria;p:Bacillota;c:Bacilli;...;s:species
#
# Ranks may stop at any level, for example:
#
#   d:Cryptista;p:Cryptophyta;c:Cryptophyceae;
#   o:Cryptomonadales;f:Cryptomonadaceae;g:Cryptomonas
#
# Missing or empty taxonomy represents an unassigned ASV.
# -------------------------------------------------------------------------

rank_codes <- c(
    "d",
    "p",
    "c",
    "o",
    "f",
    "g",
    "s"
)

rank_names <- c(
    d = "Domain",
    p = "Phylum",
    c = "Class",
    o = "Order",
    f = "Family",
    g = "Genus",
    s = "Species"
)

rank_order <- setNames(
    seq_along(rank_codes),
    rank_codes
)

# -------------------------------------------------------------------------
# Read input
# -------------------------------------------------------------------------

dat <- fread(
    input_file,
    check.names = FALSE
)

required_cols <- c(
    "id",
    "Sequence",
    "Taxon",
    "Confidence"
)

missing_cols <- setdiff(
    required_cols,
    names(dat)
)

if (length(missing_cols) > 0) {
    stop(
        sprintf(
            "Input table is missing required columns: %s",
            paste(
                missing_cols,
                collapse = ", "
            )
        )
    )
}

# -------------------------------------------------------------------------
# Identify sample abundance columns
#
# Final taxonomy table is expected to contain:
#
#   id
#   Sequence
#   Taxon
#   Confidence
#   <sample columns...>
#
# Everything after Confidence is treated as abundance.
# -------------------------------------------------------------------------

confidence_idx <- match(
    "Confidence",
    names(dat)
)

if (confidence_idx == ncol(dat)) {
    stop(
        "No sample abundance columns found after Confidence."
    )
}

sample_cols <- names(dat)[
    (confidence_idx + 1):ncol(dat)
]

if (length(sample_cols) == 0) {
    stop(
        "No sample abundance columns found."
    )
}

# -------------------------------------------------------------------------
# Validate and convert abundance columns
# -------------------------------------------------------------------------

for (col in sample_cols) {

    original <- dat[[col]]

    numeric_values <- suppressWarnings(
        as.numeric(original)
    )

    invalid <- (
        !is.na(original) &
        trimws(as.character(original)) != "" &
        is.na(numeric_values)
    )

    if (any(invalid)) {
        stop(
            sprintf(
                "Sample abundance column '%s' contains non-numeric values.",
                col
            )
        )
    }

    set(
        dat,
        j = col,
        value = numeric_values
    )
}

# Total reads represented by each ASV across all samples
dat[
    ,
    total_reads := rowSums(
        .SD,
        na.rm = TRUE
    ),
    .SDcols = sample_cols
]

# -------------------------------------------------------------------------
# Validate canonical taxonomy format
# -------------------------------------------------------------------------

validate_taxonomy <- function(taxon) {

    if (
        is.na(taxon) ||
        trimws(taxon) == "" ||
        trimws(taxon) == "Unassigned"
    ) {
        return(TRUE)
    }

    taxon <- trimws(taxon)

    parts <- strsplit(
        taxon,
        ";",
        fixed = TRUE
    )[[1]]

    parts <- trimws(parts)

    parts <- parts[
        parts != ""
    ]

    if (length(parts) == 0) {
        return(TRUE)
    }

    # Every component must have one canonical rank prefix.
    if (!all(
        grepl(
            "^[dpcofgs]:.+$",
            parts
        )
    )) {
        return(FALSE)
    }

    codes <- sub(
        ":.*$",
        "",
        parts
    )

    # Rank must not occur more than once.
    if (anyDuplicated(codes)) {
        return(FALSE)
    }

    # Ranks must occur in taxonomic order.
    positions <- unname(
        rank_order[codes]
    )

    if (is.unsorted(
        positions,
        strictly = TRUE
    )) {
        return(FALSE)
    }

    TRUE
}

valid_taxonomy <- vapply(
    dat$Taxon,
    validate_taxonomy,
    logical(1)
)

if (!all(valid_taxonomy)) {

    bad_idx <- which(
        !valid_taxonomy
    )[1]

    stop(
        sprintf(
            paste(
                "Non-canonical taxonomy detected.",
                "Expected format d:...;p:...;c:...;o:...;f:...;g:...;s:...",
                "Example offending value: %s"
            ),
            dat$Taxon[bad_idx]
        )
    )
}

# -------------------------------------------------------------------------
# Parse taxonomy
# -------------------------------------------------------------------------

extract_rank <- function(taxon, rank) {

    if (
        is.na(taxon) ||
        trimws(taxon) == "" ||
        trimws(taxon) == "Unassigned"
    ) {
        return(NA_character_)
    }

    parts <- strsplit(
        taxon,
        ";",
        fixed = TRUE
    )[[1]]

    parts <- trimws(parts)

    pattern <- paste0(
        "^",
        rank,
        ":"
    )

    hit <- grep(
        pattern,
        parts,
        value = TRUE
    )

    if (length(hit) == 0) {
        return(NA_character_)
    }

    value <- sub(
        pattern,
        "",
        hit[1]
    )

    value <- trimws(value)

    if (
        value == "" ||
        value == "NA" ||
        value == "Unassigned"
    ) {
        return(NA_character_)
    }

    value
}

for (rank in rank_codes) {

    dat[
        ,
        (rank) := vapply(
            Taxon,
            extract_rank,
            character(1),
            rank = rank
        )
    ]
}

# -------------------------------------------------------------------------
# Assignment summary per taxonomic rank
#
# This table is cumulative.
#
# An ASV assigned to genus contributes to:
#
#   Domain
#   Phylum
#   Class
#   Order
#   Family
#   Genus
#
# but not Species.
# -------------------------------------------------------------------------

total_asvs <- nrow(dat)

total_reads <- sum(
    dat$total_reads,
    na.rm = TRUE
)

assignment_summary <- rbindlist(
    lapply(
        rank_codes,
        function(rank) {

            assigned <- !is.na(
                dat[[rank]]
            )

            assigned_asvs <- sum(
                assigned
            )

            unassigned_asvs <- (
                total_asvs -
                assigned_asvs
            )

            assigned_reads <- sum(
                dat$total_reads[
                    assigned
                ],
                na.rm = TRUE
            )

            unassigned_reads <- (
                total_reads -
                assigned_reads
            )

            data.table(
                rank = unname(
                    rank_names[rank]
                ),
                total_asvs = total_asvs,
                assigned_asvs = assigned_asvs,
                unassigned_asvs = unassigned_asvs,
                pct_assigned_asvs = if (
                    total_asvs > 0
                ) {
                    100 *
                        assigned_asvs /
                        total_asvs
                } else {
                    NA_real_
                },
                total_reads = total_reads,
                assigned_reads = assigned_reads,
                unassigned_reads = unassigned_reads,
                pct_assigned_reads = if (
                    total_reads > 0
                ) {
                    100 *
                        assigned_reads /
                        total_reads
                } else {
                    NA_real_
                }
            )
        }
    )
)

assignment_summary[
    ,
    `:=`(
        pct_assigned_asvs = round(
            pct_assigned_asvs,
            2
        ),
        pct_assigned_reads = round(
            pct_assigned_reads,
            2
        )
    )
]

fwrite(
    assignment_summary,
    assignment_out,
    sep = "\t",
    quote = FALSE,
    na = ""
)

# -------------------------------------------------------------------------
# Determine deepest assigned taxonomic rank
#
# Examples:
#
#   d:...;p:...;c:...;o:...;f:...;g:...
#       -> Genus
#
#   d:...;p:...
#       -> Phylum
#
#   empty taxonomy
#       -> Unassigned
# -------------------------------------------------------------------------

get_deepest_rank <- function(...) {

    values <- list(...)

    assigned <- which(
        !vapply(
            values,
            function(x) {
                is.na(x) ||
                    trimws(x) == ""
            },
            logical(1)
        )
    )

    if (length(assigned) == 0) {
        return(
            "Unassigned"
        )
    }

    unname(
        rank_names[
            rank_codes[
                max(assigned)
            ]
        ]
    )
}

dat[
    ,
    deepest_rank := mapply(
        get_deepest_rank,
        d,
        p,
        c,
        o,
        f,
        g,
        s
    )
]

# -------------------------------------------------------------------------
# Summarize deepest assigned rank
# -------------------------------------------------------------------------

deepest_levels <- c(
    "Unassigned",
    "Domain",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
)

deepest_summary <- dat[
    ,
    .(
        asvs = .N,
        reads = sum(
            total_reads,
            na.rm = TRUE
        )
    ),
    by = deepest_rank
]

# Add zero-count ranks so the output schema is always stable.
deepest_summary <- merge(
    data.table(
        deepest_rank = deepest_levels
    ),
    deepest_summary,
    by = "deepest_rank",
    all.x = TRUE
)

deepest_summary[
    is.na(asvs),
    asvs := 0L
]

deepest_summary[
    is.na(reads),
    reads := 0
]

deepest_summary[
    ,
    `:=`(
        pct_asvs = if (
            total_asvs > 0
        ) {
            round(
                100 *
                    asvs /
                    total_asvs,
                2
            )
        } else {
            NA_real_
        },
        pct_reads = if (
            total_reads > 0
        ) {
            round(
                100 *
                    reads /
                    total_reads,
                2
            )
        } else {
            NA_real_
        }
    )
]

# Preserve biologically meaningful rank order.
deepest_summary[
    ,
    deepest_rank := factor(
        deepest_rank,
        levels = deepest_levels
    )
]

setorder(
    deepest_summary,
    deepest_rank
)

deepest_summary[
    ,
    deepest_rank := as.character(
        deepest_rank
    )
]

setcolorder(
    deepest_summary,
    c(
        "deepest_rank",
        "asvs",
        "pct_asvs",
        "reads",
        "pct_reads"
    )
)

fwrite(
    deepest_summary,
    deepest_out,
    sep = "\t",
    quote = FALSE,
    na = ""
)

message(
    "Taxonomy assignment summary written: ",
    assignment_out
)

message(
    "Deepest-rank summary written: ",
    deepest_out
)
