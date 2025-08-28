#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: homer_clean_and_bed.R in_homer_txt out_clean_tsv out_bed blacklist_path [--schema tf|histone|auto]")
}
in_txt  <- args[[1]]
out_tsv <- args[[2]]
out_bed <- args[[3]]
bl_path <- args[[4]]
schema  <- if (length(args) >= 5) tolower(sub("^--schema\\s*", "", args[[5]])) else "auto"

# ---------- Robust HOMER reader ----------
read_homer_findpeaks <- function(path){
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) stop("Empty HOMER file: ", path)
  
  # find the header line that starts with PeakID (comments optional)
  hdr_idx <- which(grepl("^\\s*#?\\s*PeakID(\\t|\\s|$)", lines))[1]
  if (is.na(hdr_idx)) {
    # sometimes there is a "# Column Headers:" line just before
    ch <- which(grepl("^\\s*#?\\s*Column Headers", lines))[1]
    if (!is.na(ch) && ch < length(lines) &&
        grepl("^\\s*#?\\s*PeakID(\\t|\\s|$)", lines[ch + 1])) {
      hdr_idx <- ch + 1
    } else {
      stop("Could not find 'PeakID' header in HOMER file: ", path)
    }
  }
  
  header_line <- sub("^\\s*#\\s*", "", lines[hdr_idx])
  data_lines  <- lines[(hdr_idx + 1):length(lines)]
  # drop comments and blank lines after header
  data_lines  <- data_lines[!grepl("^\\s*#", data_lines)]
  data_lines  <- data_lines[nzchar(trimws(data_lines))]
  
  tsv <- paste0(header_line, "\n", paste(data_lines, collapse = "\n"))
  df  <- suppressMessages(readr::read_tsv(I(tsv), show_col_types = FALSE, guess_max = 100000, progress = FALSE))
  
  # normalize names: spaces/hyphens -> underscore, tolower
  names(df) <- names(df) |>
    stringr::str_replace_all(c(" " = "_", "," = "", "-" = "_")) |>
    tolower()
  
  df
}

# ---------- Blacklist reader (CSV or BED) ----------
read_blacklist <- function(path){
  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("bed","")) {
    bl <- suppressMessages(readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE, progress = FALSE))
    if (ncol(bl) < 3) stop("Blacklist BED must have at least 3 columns.")
    bl <- bl[, 1:3]; names(bl) <- c("Chr","Start","End")
  } else {
    bl <- suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
    nms <- tolower(names(bl))
    cidx <- which(nms %in% c("chr","chrom","chromosome"))[1]
    sidx <- which(nms %in% c("start","start_bp","startpos"))[1]
    eidx <- which(nms %in% c("end","end_bp","endpos","stop"))[1]
    if (any(is.na(c(cidx,sidx,eidx)))) stop("CSV blacklist must contain chr/start/end columns.")
    bl <- tibble(Chr = bl[[cidx]], Start = bl[[sidx]], End = bl[[eidx]])
  }
  bl %>% mutate(Start = as.integer(Start), End = as.integer(End))
}

# ---------- Overlap + filter ----------
overlap_any <- function(a_start, a_end, b_start, b_end) (a_end >= b_start & a_start <= b_end)

filter_blacklist <- function(df, bl){
  req <- c("chr","start","end")
  if (!all(req %in% names(df))) stop("Input must contain chr/start/end columns (HOMER findPeaks output).")
  if (nrow(df) == 0) return(df)
  
  df <- df %>% mutate(start = as.integer(start), end = as.integer(end))
  keep <- logical(nrow(df))
  idx_by_chr <- split(seq_len(nrow(df)), df$chr)
  
  for (chr in names(idx_by_chr)) {
    idx <- idx_by_chr[[chr]]
    blc <- bl %>% filter(Chr == chr)
    if (nrow(blc) == 0) { keep[idx] <- TRUE; next }
    si <- df$start[idx]; ei <- df$end[idx]
    no_overlap <- rep(TRUE, length(idx))
    for (k in seq_len(nrow(blc))) {
      no_overlap <- no_overlap & !overlap_any(si, ei, blc$Start[k], blc$End[k])
      if (!any(no_overlap)) break
    }
    keep[idx] <- no_overlap
  }
  df[keep, , drop = FALSE]
}

# ---------- Main ----------
schema <- if (schema %in% c("tf","histone")) schema else "auto"
homer  <- read_homer_findpeaks(in_txt)

# auto-detect schema (we keep it for compatibility; behavior is same)
has_focus  <- "focus_ratio" %in% names(homer) || "p_value_vs_local" %in% names(homer)
has_region <- "region_size" %in% names(homer) && !"focus_ratio" %in% names(homer)
if (schema == "auto") schema <- if (has_region && !has_focus) "histone" else "tf"

bl     <- read_blacklist(bl_path)
clean  <- filter_blacklist(homer, bl)

# write cleaned TSV (full table with normalized names)
readr::write_tsv(clean, out_tsv)

# write 3-col BED
bed <- clean %>% select(chr, start, end)
readr::write_tsv(bed, out_bed, col_names = FALSE)


