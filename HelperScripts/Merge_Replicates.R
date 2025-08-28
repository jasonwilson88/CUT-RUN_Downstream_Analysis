#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
})

# --------------------------------------------------------------------
# Usage: merge_replicates.R out_overlap_tsv rep1.tsv rep2.tsv
# Env:
#   OVERLAP_PAD_BP             (default 0)      # expand intervals ±bp for overlap test
#   OVERLAP_REQUIRE_SAME_GENE  (default 0)      # 1 to require same gene label across reps
#   MERGE_VERBOSE              (default 1)      # 1 to print diagnostics
# --------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: merge_replicates.R out_overlap_tsv rep1.tsv rep2.tsv")
out_overlap <- args[1]; rep1_path <- args[2]; rep2_path <- args[3]

PAD_BP   <- as.integer(Sys.getenv("OVERLAP_PAD_BP", "0"))
REQ_SAME <- as.integer(Sys.getenv("OVERLAP_REQUIRE_SAME_GENE", "0")) == 1
VERB     <- as.integer(Sys.getenv("MERGE_VERBOSE", "1")) == 1
logmsg   <- function(...) if (VERB) message(...)

# ---------- helpers (for matching only; never rename user headers) ----------
norm_soft <- function(x) tolower(gsub("[\\s,\\-]+","_", x))
norm_flat <- function(x) tolower(gsub("[^a-z0-9]+","", x))

pick_idx <- function(nflat, patterns) {
  hit <- integer(0)
  for (p in patterns) hit <- c(hit, grep(p, nflat, perl = TRUE))
  u <- unique(hit)
  if (length(u)) u[1] else NA_integer_
}

# FC/P: accept unsuffixed or suffixed (...control1/_1), or compact forms (fc, pvalue)
find_fc_idx <- function(nflat) {
  idx <- pick_idx(nflat, c("^foldchangevscontrol$", "^foldchangevscontrol[0-9]+$"))
  if (is.na(idx)) idx <- pick_idx(nflat, c("^fc[0-9]*$"))
  idx
}
find_p_idx <- function(nflat) {
  idx <- pick_idx(nflat, c("^pvaluevscontrol$", "^pvaluevscontrol[0-9]+$"))
  if (is.na(idx)) idx <- pick_idx(nflat, c("^pvalue[0-9]*$"))
  idx
}

# Choose a gene label column by coalescing across likely candidates (preserve original case)
extract_gene_vector <- function(df) {
  nm <- names(df)
  nf <- norm_flat(nm)
  
  # candidates in priority order
  cand_sets <- list(
    grep("^genename(\\.|$)", nf),
    grep("^genesymbol(\\.|$)", nf),
    grep("^symbol(\\.|$)", nf),
    grep("\\bgene(\\.|$)", nf)
  )
  # exclude bad types from each set
  is_bad <- function(ix) {
    any(grepl("alias|nearest|unigene|refseq|ensembl|description|type|promoter|distance|tss",
              norm_soft(nm[ix])))
  }
  
  cols <- integer(0)
  for (s in cand_sets) {
    s <- s[!sapply(s, is_bad)]
    cols <- c(cols, s)
  }
  cols <- unique(cols)
  if (!length(cols)) return(rep(NA_character_, nrow(df)))
  
  # Coalesce across all candidate columns, row-wise
  out <- rep(NA_character_, nrow(df))
  for (ci in cols) {
    v <- df[[nm[ci]]]
    v <- if (is.factor(v)) as.character(v) else as.character(v)
    v <- trimws(v)
    out[is.na(out) | out==""] <- v[is.na(out) | out==""]
    if (all(!is.na(out) & out != "")) break
  }
  # final clean
  out[out==""] <- NA_character_
  out
}

read_rep <- function(path) {
  df <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE))
  if (!nrow(df)) return(NULL)
  nm <- names(df); nf <- norm_flat(nm)
  
  i_chr   <- pick_idx(nf, c("^chr$","^chrom$","^chromosome$"))
  i_start <- pick_idx(nf, c("^start$","^startbp$","^startpos$","^startposition$"))
  i_end   <- pick_idx(nf, c("^end$","^endbp$","^endpos$","^stop$","^endposition$"))
  if (any(is.na(c(i_chr, i_start, i_end))))
    stop("Missing chr/start/end in ", path, " names=", paste(nm, collapse=", "))
  
  i_fc <- find_fc_idx(nf)
  i_p  <- find_p_idx(nf)
  if (is.na(i_fc) || is.na(i_p))
    stop("Missing fold_change_vs_controlN / p_value_vs_controlN in ", path,
         " names=", paste(nm, collapse=", "))
  
  gene_vec <- extract_gene_vector(df)
  
  tibble::tibble(
    chr   = trimws(as.character(df[[nm[i_chr]]])),  # chr as character
    start = suppressWarnings(as.integer(df[[nm[i_start]]])),
    end   = suppressWarnings(as.integer(df[[nm[i_end]]])),
    gene_label = gene_vec,
    gene_key   = tolower(trimws(coalesce(gene_vec, NA_character_))),
    fc = suppressWarnings(as.numeric(df[[nm[i_fc]]])),
    p  = suppressWarnings(as.numeric(df[[nm[i_p]]])),
    # keep original coords (per-rep)
    start_orig = suppressWarnings(as.integer(df[[nm[i_start]]])),
    end_orig   = suppressWarnings(as.integer(df[[nm[i_end]]]))
  ) %>% filter(!is.na(chr), !is.na(start), !is.na(end))
}

# Expand intervals by PAD_BP for overlap test (but preserve originals)
expand_for_overlap <- function(d) {
  if (PAD_BP > 0) {
    d$start_pad <- pmax(0L, d$start - PAD_BP)
    d$end_pad   <- d$end + PAD_BP
  } else {
    d$start_pad <- d$start
    d$end_pad   <- d$end
  }
  d
}

# Emit ALL overlapping pairs in pure R (two sorted lists per chromosome)
emit_overlaps <- function(A, B) {
  out <- vector("list", 0)
  achrs <- unique(A$chr); bchrs <- unique(B$chr)
  for (chr in intersect(achrs, bchrs)) {
    a <- A[A$chr == chr, , drop = FALSE]
    b <- B[B$chr == chr, , drop = FALSE]
    if (!nrow(a) || !nrow(b)) next
    ia <- 1; ib <- 1
    while (ia <= nrow(a) && ib <= nrow(b)) {
      if (a$end_pad[ia] < b$start_pad[ib]) { ia <- ia + 1; next }
      if (b$end_pad[ib] < a$start_pad[ia]) { ib <- ib + 1; next }
      # overlap region exists
      jb <- ib
      while (jb <= nrow(b) && b$start_pad[jb] <= a$end_pad[ia]) {
        if (b$end_pad[jb] >= a$start_pad[ia]) {
          # record one pair (keep per-rep originals)
          out[[length(out)+1]] <- data.frame(
            chr = chr,
            start1 = a$start_orig[ia],
            end1   = a$end_orig[ia],
            start2 = b$start_orig[jb],
            end2   = b$end_orig[jb],
            start  = min(a$start_orig[ia], b$start_orig[jb], na.rm = TRUE),
            end    = max(a$end_orig[ia],   b$end_orig[jb],   na.rm = TRUE),
            gene_name = dplyr::coalesce(
              ifelse(is.na(a$gene_label[ia]) || a$gene_label[ia]=="", NA_character_, as.character(a$gene_label[ia])),
              ifelse(is.na(b$gene_label[jb]) || b$gene_label[jb]=="", NA_character_, as.character(b$gene_label[jb]))
            ),
            fold_change_vs_control1 = a$fc[ia],
            p_value_vs_control1     = a$p[ia],
            fold_change_vs_control2 = b$fc[jb],
            p_value_vs_control2     = b$p[jb],
            stringsAsFactors = FALSE
          )
        }
        jb <- jb + 1
      }
      # advance the one that ends first
      if (a$end_pad[ia] <= b$end_pad[ib]) ia <- ia + 1 else ib <- ib + 1
    }
  }
  if (!length(out)) NULL else dplyr::bind_rows(out)
}

# =============================== Main ===============================
rep1 <- read_rep(rep1_path)
rep2 <- read_rep(rep2_path)
logmsg(sprintf("[merge] %s rows=%s", basename(rep1_path), nrow(rep1)))
logmsg(sprintf("[merge] %s rows=%s", basename(rep2_path), nrow(rep2)))

# Optional same-gene prefilter
if (REQ_SAME) {
  k1 <- unique(rep1$gene_key[!is.na(rep1$gene_key) & nzchar(rep1$gene_key)])
  k2 <- unique(rep2$gene_key[!is.na(rep2$gene_key) & nzchar(rep2$gene_key)])
  keep_keys <- intersect(k1, k2)
  rep1 <- rep1[rep1$gene_key %in% keep_keys, , drop = FALSE]
  rep2 <- rep2[rep2$gene_key %in% keep_keys, , drop = FALSE]
  logmsg(sprintf("[merge] same-gene filter -> rep1=%s, rep2=%s", nrow(rep1), nrow(rep2)))
}

# Prepare & overlap
rep1 <- rep1 %>% arrange(chr, start, end) %>% expand_for_overlap()
rep2 <- rep2 %>% arrange(chr, start, end) %>% expand_for_overlap()
ov   <- emit_overlaps(rep1, rep2)

if (is.null(ov) || !nrow(ov)) {
  # write empty with header (plotter can still run gracefully)
  cols <- c("chr","start","end","gene_name",
            "start1","end1","start2","end2",
            "fold_change_vs_control1","p_value_vs_control1",
            "fold_change_vs_control2","p_value_vs_control2")
  empty <- as.data.frame(matrix(ncol=length(cols), nrow=0)); names(empty) <- cols
  readr::write_tsv(empty, out_overlap); quit(save="no", status=0)
}

# Deduplicate exact duplicate pairs if any; order nicely
ov <- ov %>%
  distinct(chr, start1, end1, start2, end2, .keep_all = TRUE) %>%
  arrange(chr, start, end)

readr::write_tsv(ov, out_overlap)

