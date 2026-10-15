#!/usr/bin/env bash
#
# export_app_results.sh - copy only the files the HiFi Amplicon Results
# Explorer (companion Shiny app) reads from a HiFi-16S-workflow results
# directory (--outdir), preserving the relative directory structure.
#
# Intended location in the pipeline repository: utils/export_app_results.sh
# (not bin/, which Nextflow puts on the PATH of every task).
#
# The point is to let users share a small subset of a completed run instead of
# the full tree (the example 16S run is ~8.6 GB, almost all of it FASTQs, RDS
# objects and per-read stats that the app never touches).
#
# Usage:
#   utils/export_app_results.sh [--dry-run] [--force] [-v|--verbose] \
#       <results_dir> <output_dir>
#   utils/export_app_results.sh -h|--help
#
# The file list mirrors the app's readers (R/load_results.R and the R/io_*.R
# modules of the app repository). See EXPORT_REQUIRED / EXPORT_OPTIONAL_FILES /
# EXPORT_GLOBS below; each entry is annotated with the reader that consumes it.
# If the app starts reading new files, update this list.
#
# POSIX tooling only (find, cp, mkdir, du, sort, date); no rsync.

set -euo pipefail

SCRIPT_NAME="export_app_results.sh"
SCRIPT_VERSION="1.0.0"

usage() {
  cat <<'EOF'
Usage: export_app_results.sh [--dry-run] [--force] [-v|--verbose] <results_dir> <output_dir>

Copy only the files the HiFi Amplicon Results Explorer reads from a pipeline
results directory, preserving relative paths.

Options:
  --dry-run      List what would be copied; write nothing.
  --force        Allow writing into an existing non-empty output_dir.
  -v, --verbose  List every copied file.
  -h, --help     Show this help.

Exit status:
  0  success
  1  missing/invalid results directory, or a required file is missing
  2  invalid usage or an unsafe/occupied output directory

Required files (the app cannot load without these):
  final/dada2_tracking.tsv
  final/frequency.tsv
  at least one of:
    final/best_nb_tax_merged_freq_tax.tsv
    final/vsearch_tax_merged_freq_tax.tsv

Optional files/globs are copied when present and skipped silently otherwise.
EOF
}

log()  { printf '%s\n' "$*"; }
warn() { printf '%s\n' "$*" >&2; }
die()  { local code="$1"; shift; warn "Error: $*"; exit "$code"; }

# --- required files (mirror REQUIRED_RESULT_FILES / FINAL_TAXONOMY_FILES) ----
EXPORT_REQUIRED_FILES=(
  "final/dada2_tracking.tsv"
  "final/frequency.tsv"
)
EXPORT_REQUIRED_ANY=(
  "final/best_nb_tax_merged_freq_tax.tsv"
  "final/vsearch_tax_merged_freq_tax.tsv"
)

# --- optional literal files (skipped silently when absent) ------------------
# annotated with the reader that consumes each one
EXPORT_OPTIONAL_FILES=(
  # read_asv_table() - load_results.R / io_tables.R
  "final/best_nb_tax_merged_freq_tax.tsv"
  "final/vsearch_tax_merged_freq_tax.tsv"
  # read_final_summary() - load_results.R + .final_pipeline_summaries()
  "final/best_nb_taxonomy_assignment_summary.tsv"
  "final/best_nb_taxonomy_deepest_rank_summary.tsv"
  "final/vsearch_taxonomy_assignment_summary.tsv"
  "final/vsearch_taxonomy_deepest_rank_summary.tsv"
  # read_rarefaction_depth() - load_results.R / io_tables.R
  "final/rarefaction_depth_suggested.txt"
  # read_taxonomy_db() - .read_winning_db() / io_tables.R
  "taxonomy/nb_tax/best_nb_taxonomy_withDB.tsv"
  # read_chimera_stats() - load_results.R / io_tables.R
  "dada2/chimera_removal/chimera_stats.tsv"
  # read_asv_filter_stats() - load_results.R / io_tables.R
  "dada2/filtered_asvs/filter_asv_stats.tsv"
  # read_seqkit_summary() - .read_qc_seqkit() (primary raw-read source)
  "preprocessing/reads_QC/all_samples_seqkit.summarystats.tsv"
  # read_length_quality() - load_results.R
  "preprocessing/reads_QC/seqkit.summarised_stats.group_by_samples.tsv"
  # read_execution_trace() - load_results.R / io_tables.R
  "pipeline_info/execution_trace.txt"
  # read_run_info() - load_results.R / io_execution_report.R
  "pipeline_info/execution_report.html"
)

# --- optional globs: "dir|name-pattern|exclude-regex|label" -----------------
# The exclude regex is matched against the basename with bash [[ =~ ]].
EXPORT_GLOBS=(
  # list_nb_database_files() / read_per_db_taxonomy() / compare_db_assignments()
  "taxonomy/nb_tax|*_nb.tsv|(_bootstraps\.tsv$|^best_nb_taxonomy)|per-database naive-Bayes tables"
  # .vsearch_databases() / detect_taxonomy_methods(): only the file *names* are
  # read (to recover database names); contents are not parsed.
  "taxonomy/vsearch_tax|*.vsearch_lca.tsv||per-database VSEARCH LCA tables"
  # .read_qc_seqkit(): fallback raw-read source when the combined file is absent
  "preprocessing/qc_raw|*.seqkit.summarystats.tsv||per-sample raw seqkit summaries (fallback)"
  # .read_qc_seqkit(): post-trim reads always come from per-sample files
  "preprocessing/qc_post_trim|*.post_trim.seqkit.summarystats.tsv||per-sample post-trim seqkit summaries"
  # read_cutadapt_logs() - load_results.R / io_cutadapt.R
  "preprocessing/cutadapt|*.cutadapt.log||cutadapt logs"
)

# ---------------------------------------------------------------------------

dry_run=0
force=0
verbose=0
positional=()

while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) dry_run=1 ;;
    --force)   force=1 ;;
    -v|--verbose) verbose=1 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; while [ $# -gt 0 ]; do positional+=("$1"); shift; done; break ;;
    -*) die 2 "unknown option: $1 (try --help)" ;;
    *) positional+=("$1") ;;
  esac
  shift
done

if [ "${#positional[@]}" -ne 2 ]; then
  usage >&2
  die 2 "expected exactly 2 arguments: <results_dir> <output_dir>"
fi

SRC="${positional[0]}"
OUT="${positional[1]}"

# --- canonicalise paths (works for a not-yet-existing output leaf) ----------
abs_path() {
  local target="$1" dir base
  # strip trailing slashes (but keep a lone "/")
  while [ "${target%/}" != "$target" ] && [ "$target" != "/" ]; do
    target="${target%/}"
  done
  if [ -e "$target" ]; then
    (cd -- "$target" && pwd -P)
    return 0
  fi
  dir=$(dirname -- "$target")
  base=$(basename -- "$target")
  if [ -d "$dir" ]; then
    printf '%s/%s\n' "$(cd -- "$dir" && pwd -P)" "$base"
  else
    printf '%s/%s\n' "$(abs_path "$dir")" "$base"
  fi
}

[ -d "$SRC" ] || die 1 "results directory does not exist or is not a directory: $SRC"

SRC_ABS="$(abs_path "$SRC")"
OUT_ABS="$(abs_path "$OUT")"

[ "$OUT_ABS" != "$SRC_ABS" ] || die 2 "output_dir must differ from results_dir"
case "$OUT_ABS/" in
  "$SRC_ABS"/*) die 2 "output_dir must not be inside results_dir ($SRC_ABS)" ;;
esac

# --- required-file validation ----------------------------------------------
missing_required=()
for rel in "${EXPORT_REQUIRED_FILES[@]}"; do
  [ -f "$SRC/$rel" ] || missing_required+=("$rel")
done

any_final=0
for rel in "${EXPORT_REQUIRED_ANY[@]}"; do
  if [ -f "$SRC/$rel" ]; then any_final=1; fi
done

if [ "${#missing_required[@]}" -gt 0 ] || [ "$any_final" -eq 0 ]; then
  warn "Error: $SRC_ABS does not look like a loadable results directory."
  if [ "${#missing_required[@]}" -gt 0 ]; then
    warn "  Missing required file(s):"
    for rel in "${missing_required[@]}"; do warn "    - $rel"; done
  fi
  if [ "$any_final" -eq 0 ]; then
    warn "  Missing required final taxonomy table (need at least one of):"
    for rel in "${EXPORT_REQUIRED_ANY[@]}"; do warn "    - $rel"; done
  fi
  exit 1
fi

# --- output-directory checks -----------------------------------------------
if [ -e "$OUT" ] && [ ! -d "$OUT" ]; then
  die 2 "output path exists and is not a directory: $OUT"
fi
if [ -d "$OUT" ] && [ -n "$(ls -A -- "$OUT" 2>/dev/null)" ] && [ "$force" -ne 1 ]; then
  die 2 "output_dir exists and is non-empty: $OUT (use --force to write into it)"
fi

# --- collect the files to copy ---------------------------------------------
rel_files=()
missing_optional=()
missing_groups=()

add_literal() {
  local rel="$1"
  if [ -f "$SRC/$rel" ]; then
    rel_files+=("$rel")
  else
    missing_optional+=("$rel")
  fi
}

add_glob() {
  local spec="$1"
  local dir pat excl label
  dir="${spec%%|*}"; spec="${spec#*|}"
  pat="${spec%%|*}"; spec="${spec#*|}"
  excl="${spec%%|*}"; label="${spec#*|}"

  local d="$SRC/$dir"
  if [ ! -d "$d" ]; then
    missing_groups+=("$dir/ ($label)")
    return 0
  fi
  local found=0 f base rel
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    base="${f##*/}"
    if [ -n "$excl" ] && [[ "$base" =~ $excl ]]; then continue; fi
    rel="$dir/$base"
    rel_files+=("$rel")
    found=$((found + 1))
  done < <(find "$d" -maxdepth 1 -type f -name "$pat" | LC_ALL=C sort)
  if [ "$found" -eq 0 ]; then
    missing_groups+=("$dir/$pat ($label)")
  fi
}

for rel in "${EXPORT_REQUIRED_FILES[@]}"; do rel_files+=("$rel"); done
for rel in "${EXPORT_OPTIONAL_FILES[@]}"; do add_literal "$rel"; done
for spec in "${EXPORT_GLOBS[@]}"; do add_glob "$spec"; done

# dedupe + sort
unique_files=()
if [ "${#rel_files[@]}" -gt 0 ]; then
  while IFS= read -r line; do
    [ -n "$line" ] && unique_files+=("$line")
  done < <(printf '%s\n' "${rel_files[@]}" | LC_ALL=C sort -u)
fi

n_files="${#unique_files[@]}"

# --- dry-run ----------------------------------------------------------------
if [ "$dry_run" -eq 1 ]; then
  log "Dry run: $SRC_ABS -> $OUT_ABS"
  log "Would copy $n_files file(s):"
  for rel in "${unique_files[@]}"; do log "  $rel"; done
  if [ "${#missing_optional[@]}" -gt 0 ]; then
    log ""
    log "Optional files not present (skipped):"
    for rel in "${missing_optional[@]}"; do log "  $rel"; done
  fi
  if [ "${#missing_groups[@]}" -gt 0 ]; then
    log ""
    log "Optional groups with no matching files (skipped):"
    for g in "${missing_groups[@]}"; do log "  $g"; done
  fi
  log ""
  log "Nothing written (--dry-run)."
  exit 0
fi

# --- copy -------------------------------------------------------------------
mkdir -p -- "$OUT"
for rel in "${unique_files[@]}"; do
  dest="$OUT/$rel"
  mkdir -p -- "$(dirname -- "$dest")"
  cp -p -- "$SRC/$rel" "$dest"
  if [ "$verbose" -eq 1 ]; then log "  copied: $rel"; fi
done

# --- manifest ---------------------------------------------------------------
manifest="$OUT/EXPORT_MANIFEST.txt"
{
  printf '# HiFi Amplicon Results Explorer export manifest\n'
  printf 'script: %s\n' "$SCRIPT_NAME"
  printf 'script_version: %s\n' "$SCRIPT_VERSION"
  printf 'date_utc: %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'source: %s\n' "$SRC_ABS"
  printf 'export: %s\n' "$OUT_ABS"
  printf 'files: %s\n' "$n_files"
  printf '\n# relative_path\tsource_path\n'
  for rel in "${unique_files[@]}"; do
    printf '%s\t%s\n' "$rel" "$SRC_ABS/$rel"
  done
} > "$manifest"

# --- summary ----------------------------------------------------------------
src_size="$(du -sh -- "$SRC_ABS" 2>/dev/null | cut -f1)"
out_size="$(du -sh -- "$OUT_ABS" 2>/dev/null | cut -f1)"

log ""
log "Export complete."
log "  source:        $SRC_ABS ($src_size)"
log "  export:        $OUT_ABS ($out_size)"
log "  files copied:  $n_files"
log "  manifest:      $manifest"

if [ "${#missing_optional[@]}" -gt 0 ]; then
  log ""
  log "Optional files not present (skipped, ${#missing_optional[@]}):"
  for rel in "${missing_optional[@]}"; do log "  $rel"; done
fi
if [ "${#missing_groups[@]}" -gt 0 ]; then
  log ""
  log "Optional groups with no matching files (skipped, ${#missing_groups[@]}):"
  for g in "${missing_groups[@]}"; do log "  $g"; done
fi

exit 0
