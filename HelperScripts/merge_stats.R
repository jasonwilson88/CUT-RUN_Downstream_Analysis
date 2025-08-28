#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

# --------------------------------------------------------------------
# Usage: merge_stats.R <annotations.tsv> <stats.tsv> <out.tsv>
# Env:
#   MERGE_STATS_SHIFT_START (default "1")  # +1 to stats start (HOMER 0-based -> anno 1-based)
#   MERGE_STATS_VERBOSE     (default "1")
# --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: merge_stats.R <annotations.tsv> <stats.tsv> <out.tsv>")
anno_path  <- args[[1]]
stats_path <- args[[2]]
out_path   <- args[[3]]

SHIFT <- as.integer(Sys.getenv("MERGE_STATS_SHIFT_START", "1"))
VERB  <- as.integer(Sys.getenv("MERGE_STATS_VERBOSE", "1")) == 1
logm  <- function(...) if (VERB) message(...)

# --- helper: normalized id (match only; NEVER mutate names) ---
norm_id <- function(x) tolower(gsub("[^a-z0-9]+", "", x))

find_idx <- function(raw_names, candidates) {
  nid <- vapply(raw_names, norm_id, character(1))
  for (cand in candidates) {
    hit <- which(nid == cand)
    if (length(hit)) return(hit[1])
  }
  NA_integer_
}

# ================== Read annotation (KEEP ALL COLUMNS) ==================
anno <- suppressMessages(readr::read_tsv(anno_path, show_col_types = FALSE))
if (!nrow(anno)) stop("Empty annotations: ", anno_path)
an <- names(anno)

# Keys
ai_chr   <- find_idx(an, c("chr","chrom","chromosome"))
ai_start <- find_idx(an, c("start","startbp","startpos","startposition"))
ai_end   <- find_idx(an, c("end","endbp","endpos","stop","endposition"))

# Fallback when first header is "PeakID (cmd=...)" and Chr/Start/End are cols 2/3/4
if (any(is.na(c(ai_chr, ai_start, ai_end))) &&
    grepl("^\\s*PeakID\\s*\\(cmd=", an[1], ignore.case = TRUE) &&
    ncol(anno) >= 4) {
  ai_chr <- 2L; ai_start <- 3L; ai_end <- 4L
  logm("[merge_stats] Fallback: using columns 2/3/4 as Chr/Start/End (PeakID cmd-first header)")
}

if (any(is.na(c(ai_chr, ai_start, ai_end)))) {
  stop("Annotations missing chr/start/end. Found: ", paste(an, collapse = ", "))
}
logm(sprintf("[merge_stats] anno keys -> %s / %s / %s", an[ai_chr], an[ai_start], an[ai_end]))

# Rename ONLY the annotation gene label: "Gene Name" -> "gene_name"
gene_idx_exact <- which(names(anno) == "Gene Name")
if (length(gene_idx_exact) == 0) {
  # be slightly tolerant to case/underscore
  gene_idx_exact <- which(tolower(names(anno)) %in% c("gene name","gene_name"))
}
if (length(gene_idx_exact) == 0) {
  stop('Annotation file does not contain a "Gene Name" column.')
}
names(anno)[gene_idx_exact[1]] <- "gene_name"

# Build canonical join keys (DO NOT rename user headers beyond 'gene_name')
anno_keys <- anno %>%
  transmute(.a_chr   = .data[[an[ai_chr]]],
            .a_start = suppressWarnings(as.integer(.data[[an[ai_start]]])),
            .a_end   = suppressWarnings(as.integer(.data[[an[ai_end]]]))) %>%
  mutate(.a_row = row_number())
anno_aug <- bind_cols(anno_keys, anno)

# ================== Read stats (KEEP ALL COLUMNS) ==================
stats <- suppressMessages(readr::read_tsv(stats_path, show_col_types = FALSE))
if (!nrow(stats)) stop("Empty stats: ", stats_path)
sn <- names(stats)

si_chr   <- find_idx(sn, c("chr","chrom","chromosome"))
si_start <- find_idx(sn, c("start","startbp","startpos","startposition"))
si_end   <- find_idx(sn, c("end","endbp","endpos","stop","endposition"))

if (any(is.na(c(si_chr, si_start, si_end)))) {
  stop("Stats missing chr/start/end. Found: ", paste(sn, collapse = ", "))
}
logm(sprintf("[merge_stats] stats keys -> %s / %s / %s", sn[si_chr], sn[si_start], sn[si_end]))

# If stats also has a 'gene_name' column, rename it to avoid conflicting with annotation's gene_name
if ("gene_name" %in% names(stats)) {
  names(stats)[which(names(stats) == "gene_name")] <- "gene_name_stats"
}

stats_keys <- stats %>%
  transmute(.s_chr   = .data[[sn[si_chr]]],
            .s_start = suppressWarnings(as.integer(.data[[sn[si_start]]])),
            .s_end   = suppressWarnings(as.integer(.data[[sn[si_end]]]))) %>%
  mutate(.s_row = row_number())
if (SHIFT != 0L) stats_keys$.s_start <- stats_keys$.s_start + SHIFT

stats_aug <- bind_cols(stats_keys, stats)

logm(sprintf("[merge_stats] rows: anno=%s, stats=%s, shift=%d", nrow(anno), nrow(stats), SHIFT))

# ================== Inner join on canonical keys; KEEP ALL COLUMNS ==================
joined <- dplyr::inner_join(
  stats_aug, anno_aug,
  by = c(".s_chr" = ".a_chr", ".s_start" = ".a_start", ".s_end" = ".a_end"),
  suffix = c(".stat", ".anno")
)

# Retry without shift if empty
if (nrow(joined) == 0 && SHIFT != 0L) {
  logm("[merge_stats] 0 rows after join; retrying with SHIFT=0")
  stats_keys2 <- stats_keys; stats_keys2$.s_start <- stats_keys2$.s_start - SHIFT
  stats_aug2  <- bind_cols(stats_keys2, stats)
  joined <- dplyr::inner_join(
    stats_aug2, anno_aug,
    by = c(".s_chr" = ".a_chr", ".s_start" = ".a_start", ".s_end" = ".a_end"),
    suffix = c(".stat",".anno")
  )
}

logm(sprintf("[merge_stats] rows after join: %s", nrow(joined)))

# ================== Finalize: front keys + annotation gene_name ==================
# Clean key columns (provide canonical copies without touching originals)
joined$chr   <- joined$.s_chr
joined$start <- joined$.s_start
joined$end   <- joined$.s_end

# Ensure there is a single authoritative 'gene_name' from annotation
# (we already renamed annotation's "Gene Name" -> gene_name, so it should be present unsuffixed)
if (!"gene_name" %in% names(joined)) {
  # If dplyr suffixed it due to a stray stats 'gene_name', recover from .anno suffix
  if ("gene_name.anno" %in% names(joined)) {
    joined$gene_name <- joined[["gene_name.anno"]]
  } else {
    # last resort: create NA column (but this should not happen with the rename above)
    joined$gene_name <- NA_character_
  }
}

# Reorder: keys first, then EVERYTHING else (unaltered)
front <- c("chr","start","end","gene_name")
drop_tmp <- c(".a_row",".s_row",".s_chr",".s_start",".s_end",".a_chr",".a_start",".a_end")
rest  <- setdiff(names(joined), c(front, drop_tmp))
out <- joined %>%
  select(all_of(front), all_of(rest)) %>%
  mutate(start = as.integer(start), end = as.integer(end))

readr::write_tsv(out, out_path)

