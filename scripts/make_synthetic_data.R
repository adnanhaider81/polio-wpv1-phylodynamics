#!/usr/bin/env Rscript

set.seed(20260501)

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
root <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
out_dir <- file.path(root, "data", "synthetic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

make_base <- function(n = 906) {
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}

mutate_sequence <- function(seq, positions, bases) {
  x <- strsplit(seq, "")[[1]]
  x[positions] <- bases
  paste(x, collapse = "")
}

regions <- data.frame(
  country = c("PAK", "PAK", "PAK", "PAK", "PAK", "AFG", "AFG"),
  state = c("SINDH", "SINDH", "BALOCHISTAN", "KHYBER PAKHTUNKHWA", "PUNJAB", "HILMAND", "KANDAHAR"),
  locality = c("KHI_DEMO", "HYDERABAD_DEMO", "QUETTA_DEMO", "PESHAWAR_DEMO", "LAHORE_DEMO", "LASHKARGAH_DEMO", "KANDAHAR_DEMO"),
  cluster = c("C1", "C1", "C2", "C2", "C3", "C2", "C3"),
  n = c(7, 5, 6, 5, 5, 4, 4),
  stringsAsFactors = FALSE
)

base <- make_base(906)
records <- list()
metadata <- list()
idx <- 1

for (r in seq_len(nrow(regions))) {
  lineage_positions <- sample(seq_len(906), 12)
  lineage_bases <- sample(c("A", "C", "G", "T"), 12, replace = TRUE)
  lineage_base <- mutate_sequence(base, lineage_positions, lineage_bases)
  for (j in seq_len(regions$n[[r]])) {
    strain <- sprintf("SYN_WPV1_%03d", idx)
    extra_positions <- sample(seq_len(906), 2 + (idx %% 4))
    extra_bases <- sample(c("A", "C", "G", "T"), length(extra_positions), replace = TRUE)
    seq <- mutate_sequence(lineage_base, extra_positions, extra_bases)
    date <- as.Date("2025-01-15") + (idx - 1) * 8 + sample(0:3, 1)
    records[[strain]] <- seq
    metadata[[idx]] <- data.frame(
      FILENAME = strain,
      Country = regions$country[[r]],
      StateProv = regions$state[[r]],
      locality = regions$locality[[r]],
      onsetdate = format(date, "%Y-%m-%d"),
      IsENV = ifelse(idx %% 5 == 0, "FALSE", "TRUE"),
      Cluster = regions$cluster[[r]],
      stringsAsFactors = FALSE
    )
    idx <- idx + 1
  }
}

fasta_path <- file.path(out_dir, "synthetic_alignment.fasta")
con <- file(fasta_path, open = "w")
on.exit(close(con), add = TRUE)
for (name in names(records)) {
  writeLines(paste0(">", name), con)
  writeLines(records[[name]], con)
}

meta <- do.call(rbind, metadata)
utils::write.csv(meta, file.path(out_dir, "synthetic_metadata.csv"), row.names = FALSE, quote = TRUE)
message("Wrote synthetic data to: ", out_dir)
