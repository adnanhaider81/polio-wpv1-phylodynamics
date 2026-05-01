#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
install_optional <- "--optional" %in% args

repos <- "https://cloud.r-project.org"
core <- c(
  "optparse", "readr", "dplyr", "stringr", "lubridate",
  "ape", "phangorn", "treedater", "phytools",
  "ggplot2", "scales", "tidyr", "igraph"
)
optional_cran <- c("readxl", "circlize", "sf", "geodata", "geosphere")
optional_bioc <- c("ggtree")

install_missing_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) install.packages(missing, repos = repos)
}

install_missing_cran(core)

if (install_optional) {
  install_missing_cran(optional_cran)
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = repos)
  }
  missing_bioc <- optional_bioc[!vapply(optional_bioc, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_bioc) > 0) {
    BiocManager::install(missing_bioc, ask = FALSE, update = TRUE)
  }
}

message("Package installation check complete.")
