# scripts/download_boundaries.R
# One-time helper to pre-download GADM boundaries to a project root.
# This makes map plotting faster and repeatable.

suppressPackageStartupMessages({
  library(geodata)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 1) stop("Usage: Rscript scripts/download_boundaries.R ROOT_PATH")

ROOT <- args[[1]]
dir.create(file.path(ROOT, "data/gis"), recursive = TRUE, showWarnings = FALSE)

message("Downloading Pakistan GADM level 1 and 2...")
geodata::gadm(country="PAK", level=1, path=file.path(ROOT, "data/gis"))
geodata::gadm(country="PAK", level=2, path=file.path(ROOT, "data/gis"))

message("Downloading Afghanistan GADM level 1...")
geodata::gadm(country="AFG", level=1, path=file.path(ROOT, "data/gis"))

message("Done. Files saved under: ", file.path(ROOT, "data/gis"))
