#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(dplyr); library(readr) })

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: center_bed.R in_bed out_bed")

# Try with/without header automatically
guess <- suppressWarnings(readr::read_tsv(args[1], col_names = FALSE, show_col_types = FALSE))
if (ncol(guess) < 3) stop("Input BED must have >= 3 columns")
bed <- guess[,1:3]
names(bed) <- c("chr","start","end")

bed2 <- bed %>%
  mutate(start = as.integer(start), end = as.integer(end)) %>%
  mutate(center = as.integer(ceiling((start + end)/2))) %>%
  transmute(chr, start = center, end = center + 1)

readr::write_tsv(bed2, args[2], col_names = FALSE)
