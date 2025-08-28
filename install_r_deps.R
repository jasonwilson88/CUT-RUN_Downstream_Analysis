#!/usr/bin/env Rscript
pkgs <- c("dplyr","readr","stringr","ggplot2","ggrepel")
inst <- installed.packages()[,"Package"]
to_install <- setdiff(pkgs, inst)
if (length(to_install)) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
} else {
  message("All R deps already installed: ", paste(pkgs, collapse=", "))
}
