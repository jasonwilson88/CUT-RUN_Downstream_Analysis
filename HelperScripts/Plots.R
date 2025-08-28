#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

# ---------------- args & square size ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Plots.R overlap_tsv pval_pdf fc_pdf volcano1_pdf volcano2_pdf [--highlight G1,G2,...]")
}
overlap_path <- args[[1]]
pval_pdf     <- args[[2]]
fc_pdf       <- args[[3]]
vol1_pdf     <- args[[4]]
vol2_pdf     <- args[[5]]

# Parse --highlight as "--highlight X" or "--highlight=X"
parse_highlight <- function(argv) {
  hl <- character(0)
  if (length(argv) >= 6) {
    for (i in 6:length(argv)) {
      a <- argv[[i]]
      if (grepl("^--highlight(=|$)", a)) {
        if (a == "--highlight" && i + 1 <= length(argv)) {
          hl <- argv[[i + 1]]
        } else {
          hl <- sub("^--highlight=", "", a)
        }
        break
      }
    }
  }
  if (length(hl) && nzchar(hl)) trimws(unlist(strsplit(hl, ","))) else character(0)
}
highlight <- parse_highlight(args)

# Square output size (inches); override with env var PLOT_SIZE_IN
size <- as.numeric(Sys.getenv("PLOT_SIZE_IN", "6"))
if (!is.finite(size) || size <= 0) size <- 6

# ---------------- helpers ----------------
norm_names <- function(x) tolower(gsub("[\\s,\\-]+", "_", x))
alnum_key  <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
rep_num    <- function(nm) suppressWarnings(as.integer(sub(".*?([0-9]+)$", "\\1", nm)))

# Forcing identical numeric limits on x/y to make a true square plot
equal_limits <- function(x, y, zero_floor = FALSE, pad_frac = 0.05) {
  rng <- range(c(x, y), na.rm = TRUE)
  if (!all(is.finite(rng))) rng <- c(0, 1)
  if (zero_floor) rng[1] <- min(0, rng[1])
  span <- diff(rng)
  if (!is.finite(span) || span == 0) span <- 1
  pad <- span * pad_frac
  c(rng[1] - pad, rng[2] + pad)
}

topn_by <- function(tbl, col, n) {
  v <- tbl[[col]]
  ord <- order(replace(v, is.na(v), -Inf), decreasing = TRUE)
  head(tbl[ord, , drop = FALSE], n)
}

pick_gene_col <- function(df, norm, orig) {
  exact_priority <- c("gene_name","gene_symbol","symbol","gene","geneid","gene_id")
  for (p in exact_priority) if (p %in% norm) return(orig[match(p, norm)])
  cand <- grep("gene", norm)
  if (length(cand)) {
    excl <- grepl("description|type|nearest|detailed|annotation|promoter", norm[cand])
    cand <- cand[!excl]
  }
  if (!length(cand)) stop("No suitable gene column found.")
  lens <- sapply(cand, function(i) {
    v <- df[[orig[i]]]; v <- v[!is.na(v)]
    if (!length(v)) Inf else median(nchar(as.character(v)))
  })
  orig[cand[which.min(lens)]]
}

# Robust finder for replicate columns (FC & p)
find_rep_cols <- function(norm, kind = c("fc","pv")) {
  kind <- match.arg(kind)
  parts <- if (kind == "fc") c("fold","change","vs","control") else c("^p","value","vs","control")
  has_all <- function(nm, parts) all(sapply(parts, function(p) grepl(p, nm)))
  idx <- which(sapply(norm, function(nm) has_all(nm, parts) && grepl("[0-9]+$", nm)))
  if (length(idx) < 2) {
    stem <- if (kind == "fc") "fold.*change.*control" else "p.*value.*control"
    idx <- which(grepl(stem, norm) & grepl("[0-9]+$", norm))
  }
  if (length(idx) >= 2) idx[order(rep_num(norm[idx]))] else integer(0)
}

# ---------------- read & normalize ----------------
df <- suppressMessages(readr::read_tsv(overlap_path, show_col_types = FALSE))
if (!nrow(df)) stop("Empty overlap file: ", overlap_path)

orig_names <- names(df)
norm <- norm_names(orig_names)

gene_col <- pick_gene_col(df, norm, orig_names)
fc_idx   <- find_rep_cols(norm, "fc")
pv_idx   <- find_rep_cols(norm, "pv")

message("[plots] gene column chosen: ", gene_col)
message("[plots] found FC cols: ", if (length(fc_idx)) paste(orig_names[fc_idx], collapse=", ") else "<none>")
message("[plots] found P  cols: ", if (length(pv_idx)) paste(orig_names[pv_idx], collapse=", ") else "<none>")

if (length(fc_idx) < 2 || length(pv_idx) < 2) {
  stop(
    "Need >= 2 columns each for Fold_Change_vs_ControlN and p_value_vs_ControlN (N=1,2,...). ",
    "Saw FC=", paste(orig_names[fc_idx], collapse=", "),
    " ; P=",  paste(orig_names[pv_idx], collapse=", ")
  )
}

FC1 <- orig_names[fc_idx[1]]; FC2 <- orig_names[fc_idx[2]]
P1  <- orig_names[pv_idx[1]];  P2  <- orig_names[pv_idx[2]]

# ---------------- compute stats & keys ----------------
eps <- 1e-300
df2 <- df %>%
  mutate(
    neglog10p1 = -log10(pmax(.data[[P1]], eps)),
    neglog10p2 = -log10(pmax(.data[[P2]], eps)),
    gene_key   = tolower(trimws(.data[[gene_col]])),
    gene_key2  = alnum_key(.data[[gene_col]]),
    label      = .data[[gene_col]]          # preserve original casing for labels
  )

# alias support (optional)
alias_candidates <- c("gene_alias","gene_aliases","alias","aliases")
alias_col <- NA_character_
for (cand in alias_candidates) {
  idx <- which(norm == cand)
  if (length(idx)) { alias_col <- orig_names[idx[1]]; break }
}
if (!is.na(alias_col)) df2$alias_key2 <- alnum_key(df[[alias_col]]) else df2$alias_key2 <- NA_character_

# ---------------- highlights (case-insensitive match; labels keep original case) ----------------
hl_raw0 <- highlight
hl_raw  <- tolower(trimws(hl_raw0))
hl_key2 <- alnum_key(hl_raw0)

match_exact       <- if (length(hl_raw)) df2$gene_key  %in% hl_raw  else rep(FALSE, nrow(df2))
match_gene_sanit  <- if (length(hl_key2)) df2$gene_key2 %in% alnum_key(hl_raw0) else rep(FALSE, nrow(df2))
match_alias_sanit <- if (!is.na(alias_col) && length(hl_key2)) df2$alias_key2 %in% alnum_key(hl_raw0) else rep(FALSE, nrow(df2))

match_regex <- rep(FALSE, nrow(df2))
if (length(hl_raw0)) {
  for (h in hl_raw0) {
    if (!nzchar(h)) next
    rx <- paste0("\\b", stringr::str_replace_all(h, "([\\W_])", "\\\\\\1"), "\\b")
    match_regex <- match_regex | grepl(rx, df[[gene_col]], ignore.case = TRUE)
  }
}
hl_mask <- match_exact | match_gene_sanit | match_alias_sanit | match_regex
hl_df   <- df2[hl_mask, , drop = FALSE]

if (length(hl_raw0)) {
  found_names <- unique(df2$label[hl_mask])
  missing     <- setdiff(tolower(unique(hl_raw0)), unique(df2$gene_key[hl_mask]))
  message(sprintf("[plots] highlight asked : %s", if (length(hl_raw0)) paste(unique(hl_raw0), collapse=", ") else "<none>"))
  message(sprintf("[plots] highlight found : %s", if (length(found_names)) paste(found_names, collapse=", ") else "<none>"))
  if (length(missing)) message(sprintf("[plots] highlight NOT found: %s", paste(missing, collapse=", ")))
}

# ---------------- label sets (top-N + always include highlights) ----------------
p1_tbl <- data.frame(label = df2$label, key = df2$gene_key, val = df2$neglog10p1, stringsAsFactors = FALSE)
p2_tbl <- data.frame(label = df2$label, key = df2$gene_key, val = df2$neglog10p2, stringsAsFactors = FALSE)
f1_tbl <- data.frame(label = df2$label, key = df2$gene_key, val = df[[FC1]],    stringsAsFactors = FALSE)
f2_tbl <- data.frame(label = df2$label, key = df2$gene_key, val = df[[FC2]],    stringsAsFactors = FALSE)

pval_top_keys <- unique(c(topn_by(p1_tbl, "val", 7)$key, topn_by(p2_tbl, "val", 7)$key, df2$gene_key[hl_mask]))
fc_top_keys   <- unique(c(topn_by(f1_tbl, "val", 7)$key, topn_by(f2_tbl, "val", 7)$key, df2$gene_key[hl_mask]))
v1_top_keys   <- unique(c(topn_by(f1_tbl, "val", 5)$key, topn_by(p1_tbl, "val", 5)$key, df2$gene_key[hl_mask]))
v2_top_keys   <- unique(c(topn_by(f2_tbl, "val", 5)$key, topn_by(p2_tbl, "val", 5)$key, df2$gene_key[hl_mask]))

pval_df <- df2[df2$gene_key %in% pval_top_keys, , drop = FALSE]
fc_df   <- df2[df2$gene_key %in% fc_top_keys,   , drop = FALSE]
v1_df   <- df2[df2$gene_key %in% v1_top_keys,   , drop = FALSE]
v2_df   <- df2[df2$gene_key %in% v2_top_keys,   , drop = FALSE]

# ---------------- compute unified limits for square axes ----------------
lim_pval <- equal_limits(df2$neglog10p1, df2$neglog10p2, zero_floor = TRUE,  pad_frac = 0.05)
lim_fc   <- equal_limits(df[[FC1]],      df[[FC2]],      zero_floor = TRUE,  pad_frac = 0.05)
lim_v1   <- equal_limits(df[[FC1]],      df2$neglog10p1, zero_floor = TRUE,  pad_frac = 0.05)
lim_v2   <- equal_limits(df[[FC2]],      df2$neglog10p2, zero_floor = TRUE,  pad_frac = 0.05)

# ---------------- plots (square geometry & square ranges) ----------------
# P-value vs P-value
gp <- ggplot(df2, aes(neglog10p1, neglog10p2)) +
  geom_point() +
  ggrepel::geom_text_repel(data = pval_df, aes(label = label), size = 3) +
  geom_point(data = hl_df, aes(neglog10p1, neglog10p2), inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = hl_df, aes(neglog10p1, neglog10p2, label = label),
                           inherit.aes = FALSE, size = 3, min.segment.length = 0, max.overlaps = Inf) +
  scale_x_continuous(limits = lim_pval) +
  scale_y_continuous(limits = lim_pval) +
  coord_fixed() + theme_bw() +
  labs(x = "Rep 1 -log10 p-value", y = "Rep 2 -log10 p-value")
ggsave(pval_pdf, gp, width = size, height = size, units = "in")

# FC vs FC
gf <- ggplot(df2, aes(.data[[FC1]], .data[[FC2]])) +
  geom_point() +
  ggrepel::geom_text_repel(data = fc_df, aes(label = label), size = 3) +
  geom_point(data = hl_df, aes(.data[[FC1]], .data[[FC2]]), inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = hl_df, aes(.data[[FC1]], .data[[FC2]], label = label),
                           inherit.aes = FALSE, size = 3, min.segment.length = 0, max.overlaps = Inf) +
  scale_x_continuous(limits = lim_fc) +
  scale_y_continuous(limits = lim_fc) +
  coord_fixed() + theme_bw() +
  labs(x = "Rep 1 Fold Change", y = "Rep 2 Fold Change")
ggsave(fc_pdf, gf, width = size, height = size, units = "in")

# Volcano 1 (force equal numeric ranges across X/Y)
gv1 <- ggplot(df2, aes(.data[[FC1]], neglog10p1)) +
  geom_point() +
  ggrepel::geom_text_repel(data = v1_df, aes(label = label), size = 3) +
  geom_point(data = hl_df, aes(.data[[FC1]], neglog10p1), inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = hl_df, aes(.data[[FC1]], neglog10p1, label = label),
                           inherit.aes = FALSE, size = 3, min.segment.length = 0, max.overlaps = Inf) +
  scale_x_continuous(limits = lim_v1) +
  scale_y_continuous(limits = lim_v1) +
  coord_fixed() + theme_bw() +
  labs(x = "Rep 1 Fold Change", y = "Rep 1 -log10 p-value")
ggsave(vol1_pdf, gv1, width = size, height = size, units = "in")

# Volcano 2 (force equal numeric ranges across X/Y)
gv2 <- ggplot(df2, aes(.data[[FC2]], neglog10p2)) +
  geom_point() +
  ggrepel::geom_text_repel(data = v2_df, aes(label = label), size = 3) +
  geom_point(data = hl_df, aes(.data[[FC2]], neglog10p2), inherit.aes = FALSE) +
  ggrepel::geom_text_repel(data = hl_df, aes(.data[[FC2]], neglog10p2, label = label),
                           inherit.aes = FALSE, size = 3, min.segment.length = 0, max.overlaps = Inf) +
  scale_x_continuous(limits = lim_v2) +
  scale_y_continuous(limits = lim_v2) +
  coord_fixed() + theme_bw() +
  labs(x = "Rep 2 Fold Change", y = "Rep 2 -log10 p-value")
ggsave(vol2_pdf, gv2, width = size, height = size, units = "in")
