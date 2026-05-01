#!/usr/bin/env Rscript

# WPV1 phylogeographic analysis pipeline.
# The public repository ships with synthetic example data only.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]])))
} else {
  getwd()
}
setwd(script_dir)

msg <- function(...) cat(sprintf(...), "\n")

require_or_install_cran <- function(pkgs, install = FALSE, repos = "https://cloud.r-project.org") {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  if (!install) {
    stop(paste0(
      "Missing R packages: ", paste(missing, collapse = ", "), "\n",
      "Run: Rscript scripts/install_packages.R\n",
      "Or rerun this pipeline with --install_missing_pkgs"
    ), call. = FALSE)
  }
  install.packages(missing, repos = repos)
  invisible(TRUE)
}

require_or_install_bioc <- function(pkgs, install = FALSE) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) == 0) return(invisible(TRUE))
  if (!install) {
    stop(paste0(
      "Missing Bioconductor packages: ", paste(missing, collapse = ", "), "\n",
      "Run: Rscript scripts/install_packages.R --optional"
    ), call. = FALSE)
  }
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(pkgs, ask = FALSE, update = TRUE)
  invisible(TRUE)
}

check_exe <- function(candidates) {
  hits <- Sys.which(candidates)
  hits <- hits[nzchar(hits)]
  if (length(hits) == 0) return("")
  hits[[1]]
}

bootstrap_args <- commandArgs(trailingOnly = TRUE)
bootstrap_install <- "--install_missing_pkgs" %in% bootstrap_args
if (!requireNamespace("optparse", quietly = TRUE)) {
  if (bootstrap_install) {
    install.packages("optparse", repos = "https://cloud.r-project.org")
  } else {
    stop("Missing R package: optparse. Run: Rscript scripts/install_packages.R", call. = FALSE)
  }
}

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--root"), type = "character", default = script_dir, help = "Project root [default: repository root]"),
  make_option(c("--mode"), type = "character", default = "distance", help = "distance, fasttree, or iqtree [default: distance]"),
  make_option(c("--fasta"), type = "character", default = "data/synthetic/synthetic_alignment.fasta", help = "FASTA path relative to --root"),
  make_option(c("--metadata"), type = "character", default = "data/synthetic/synthetic_metadata.csv", help = "Metadata CSV/TSV/XLSX path relative to --root"),
  make_option(c("--region_scheme"), type = "character", default = "region10", help = "Geographic grouping scheme [default: region10]"),
  make_option(c("--max_per_region_month"), type = "integer", default = 0, help = "0 disables downsampling"),
  make_option(c("--nsim"), type = "integer", default = 25, help = "Number of stochastic-map replicates"),
  make_option(c("--seed"), type = "integer", default = 1, help = "Random seed"),
  make_option(c("--seq_len"), type = "integer", default = 0, help = "Alignment length for treedater; 0 auto-detects"),
  make_option(c("--force_tree"), action = "store_true", default = FALSE, help = "Rebuild tree and downstream tree products"),
  make_option(c("--install_missing_pkgs"), action = "store_true", default = FALSE, help = "Install missing R packages"),
  make_option(c("--make_maps"), action = "store_true", default = FALSE, help = "Build GIS maps when boundary files are available"),
  make_option(c("--download_boundaries"), action = "store_true", default = FALSE, help = "Download GADM boundaries for map generation"),
  make_option(c("--skip_figures"), action = "store_true", default = FALSE, help = "Skip PNG figure generation"),
  make_option(c("--ufboot"), type = "integer", default = 1000, help = "IQ-TREE ultrafast bootstrap replicates"),
  make_option(c("--established_months"), type = "integer", default = 6, help = "Cluster persistence threshold in months"),
  make_option(c("--program_recent_months"), type = "integer", default = 3, help = "Recent detection window in months")
)
opt <- parse_args(OptionParser(option_list = option_list))

ROOT <- normalizePath(opt$root, mustWork = FALSE)

pkgs_core <- c(
  "readr", "dplyr", "stringr", "lubridate",
  "ape", "phangorn", "treedater", "phytools",
  "ggplot2", "scales", "tidyr", "igraph"
)
pkgs_optional <- character(0)
if (grepl("\\.xlsx?$", opt$metadata, ignore.case = TRUE)) {
  pkgs_optional <- c(pkgs_optional, "readxl")
}
if (opt$make_maps || opt$download_boundaries) {
  pkgs_optional <- c(pkgs_optional, "sf", "geodata", "geosphere")
}

require_or_install_cran(unique(c(pkgs_core, pkgs_optional)), install = opt$install_missing_pkgs)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

need_files <- c(
  "R/utils.R", "R/regions.R", "R/prep.R", "R/tree.R",
  "R/ancestral.R", "R/ltl.R", "R/figures.R"
)
if (opt$make_maps || opt$download_boundaries) {
  need_files <- c(need_files, "R/map.R")
}
missing_local <- need_files[!file.exists(file.path(script_dir, need_files))]
if (length(missing_local) > 0) {
  stop("Missing pipeline files:\n- ", paste(missing_local, collapse = "\n- "), call. = FALSE)
}

source("R/utils.R")
source("R/regions.R")
source("R/prep.R")
source("R/tree.R")
source("R/ancestral.R")
source("R/ltl.R")
source("R/figures.R")
if (opt$make_maps || opt$download_boundaries) {
  source("R/map.R")
}

dir_create(ROOT)
dirs <- c("data/processed", "data/gis", "results/trees", "results/tables", "results/figures", "logs", "tmp")
invisible(lapply(file.path(ROOT, dirs), dir_create))

fasta_in <- file.path(ROOT, opt$fasta)
meta_in <- file.path(ROOT, opt$metadata)
stop_if_missing(fasta_in)
stop_if_missing(meta_in)

if (opt$mode == "fasttree") {
  if (!nzchar(check_exe(c("FastTree", "fasttree")))) {
    stop("FastTree was not found on PATH. Use --mode distance for the synthetic demo.", call. = FALSE)
  }
} else if (opt$mode == "iqtree") {
  if (!nzchar(check_exe(c("iqtree2", "iqtree3", "iqtree")))) {
    stop("IQ-TREE was not found on PATH. Use --mode distance for the synthetic demo.", call. = FALSE)
  }
  if (opt$ufboot < 1000) stop("--ufboot must be at least 1000 for IQ-TREE ultrafast bootstrap.", call. = FALSE)
} else if (opt$mode != "distance") {
  stop("Unknown --mode: ", opt$mode, call. = FALSE)
}

msg("Reading metadata...")
meta <- read_metadata(meta_in)
meta <- apply_region_scheme(meta, scheme = opt$region_scheme)
validate_alignment_metadata(fasta_in, meta)
meta_all <- meta

if (opt$max_per_region_month > 0) {
  msg("Downsampling by region-month to cap=%d", opt$max_per_region_month)
  meta <- downsample_by_region_month(meta, opt$max_per_region_month, seed = opt$seed)
}

ns_path <- file.path(ROOT, "data/processed/metadata_nextstrain.tsv")
write_nextstrain_metadata(meta, ns_path)

aln_path <- file.path(ROOT, "data/processed/alignment.fasta")
file.copy(fasta_in, aln_path, overwrite = TRUE)
alignment_bp <- if (opt$seq_len > 0) opt$seq_len else get_alignment_length(aln_path)

if (opt$mode == "distance") {
  tree_path <- file.path(ROOT, "results/trees/distance.tree")
  if (!opt$force_tree && file.exists(tree_path)) {
    msg("Using existing distance tree: %s", tree_path)
  } else {
    run_distance_tree(aln_path, tree_path)
  }
} else if (opt$mode == "fasttree") {
  tree_path <- file.path(ROOT, "results/trees/ml.tree")
  if (!opt$force_tree && file.exists(tree_path)) {
    msg("Using existing FastTree tree: %s", tree_path)
  } else {
    run_fasttree(aln_path, tree_path)
  }
} else {
  out_prefix <- file.path(ROOT, "results/trees/iqtree_run")
  tree_path <- paste0(out_prefix, ".treefile")
  if (!opt$force_tree && file.exists(tree_path)) {
    msg("Using existing IQ-TREE tree: %s", tree_path)
  } else {
    Sys.setenv(POLIO_UFBOOT = as.character(opt$ufboot))
    run_iqtree(aln_path, out_prefix)
  }
}

diagnostics <- root_to_tip_diagnostics(ape::read.tree(tree_path), meta)
write_csv(diagnostics, file.path(ROOT, "results/tables/root_to_tip_diagnostics.csv"))

dated_rds <- file.path(ROOT, "results/trees/dated_tree.rds")
dated_nwk <- file.path(ROOT, "results/trees/dated_tree.nwk")
if (!opt$force_tree && file.exists(dated_rds) && file.exists(dated_nwk)) {
  msg("Using existing dated tree: %s", dated_rds)
  dated_tree <- readRDS(dated_rds)
  meta2 <- meta %>% filter(strain %in% dated_tree$tip.label)
  meta2 <- meta2[match(dated_tree$tip.label, meta2$strain), ]
} else {
  msg("Time-scaling tree with treedater...")
  td <- time_scale_tree(tree_path, meta, seq_len = alignment_bp)
  dated_tree <- td$tree
  meta2 <- td$meta
  saveRDS(dated_tree, dated_rds)
  ape::write.tree(dated_tree, dated_nwk)
}

simmap_rds <- file.path(ROOT, "results/trees/simmap_list.rds")
if (!opt$force_tree && file.exists(simmap_rds)) {
  msg("Using existing simmap results: %s", simmap_rds)
  sims <- readRDS(simmap_rds)
} else {
  msg("Running simmap with nsim=%d", opt$nsim)
  tip_states <- setNames(meta2$region, meta2$strain)
  sims <- run_simmap(dated_tree, tip_states, nsim = opt$nsim, seed = opt$seed)
  saveRDS(sims, simmap_rds)
}

tip_states <- setNames(meta2$region, meta2$strain)
mats <- transition_counts(sims, tip_states)
summ <- summarize_transition_matrix(mats)
support <- transition_support(mats)

write_matrix_csv(summ$mean, file.path(ROOT, "results/tables/transition_matrix_mean.csv"))
write_matrix_csv(summ$lo, file.path(ROOT, "results/tables/transition_matrix_lo95.csv"))
write_matrix_csv(summ$hi, file.path(ROOT, "results/tables/transition_matrix_hi95.csv"))
write_matrix_csv(support, file.path(ROOT, "results/tables/transition_support.csv"))
write_csv(top_corridors(summ$mean, support), file.path(ROOT, "results/tables/top_corridors.csv"))
write_csv(flow_partition(summ$mean), file.path(ROOT, "results/tables/flow_partition.csv"))

ie_mean <- import_export_from_matrix(summ$mean)
ie_lo <- import_export_from_matrix(summ$lo)
ie_hi <- import_export_from_matrix(summ$hi)
ie <- ie_mean %>%
  rename(imports_mean = imports, exports_mean = exports) %>%
  left_join(rename(ie_lo, imports_lo = imports, exports_lo = exports), by = "region") %>%
  left_join(rename(ie_hi, imports_hi = imports, exports_hi = exports), by = "region") %>%
  arrange(desc(exports_mean))
write_csv(ie, file.path(ROOT, "results/tables/import_export.csv"))

ltl_obj <- extract_lt_lines(sims[[1]], meta2, censor_months = 6)
ltl_class <- classify_lt_lines(ltl_obj$ltl)
write_csv(ltl_class, file.path(ROOT, "results/tables/ltl_summary.csv"))
write_csv(ltl_replicate_summary(sims, meta2), file.path(ROOT, "results/tables/ltl_replicate_summary.csv"))

cluster_status <- summarize_cluster_district_status(
  meta_all,
  sustain_months = opt$established_months,
  recent_months = opt$program_recent_months
)
cluster_status_by_region <- summarize_cluster_status_by_region(cluster_status)
cluster_status_overall <- cluster_status %>%
  count(status, status_label, name = "n") %>%
  arrange(match(status, c(
    "Established_active", "Established_quiet",
    "Recent_under_6m", "Died_off_pre_established"
  )))

write_csv(cluster_status, file.path(ROOT, "results/tables/program_cluster_status.csv"))
write_csv(cluster_status_by_region, file.path(ROOT, "results/tables/program_cluster_status_by_region.csv"))
write_csv(cluster_status_overall, file.path(ROOT, "results/tables/program_cluster_status_overall.csv"))

if (!opt$skip_figures) {
  plot_tree(dated_tree, meta2, file.path(ROOT, "results/figures/fig_tree.png"))
  plot_ltt(dated_tree, censor_days = 180, out_png = file.path(ROOT, "results/figures/fig_ltt.png"))
  plot_import_export(ie, file.path(ROOT, "results/figures/fig_import_export.png"))
  plot_chord(summ$mean, support_mat = support, min_mean = 0.1, out_png = file.path(ROOT, "results/figures/fig_chord.png"))
  plot_ltl_facets(ltl_class, file.path(ROOT, "results/figures/fig_ltl.png"))
  plot_program_cluster_status(
    cluster_status,
    file.path(ROOT, "results/figures/fig_program_cluster_status.png"),
    sustain_months = opt$established_months,
    recent_months = opt$program_recent_months
  )
}

if (opt$make_maps || opt$download_boundaries) {
  have_gis <- length(list.files(file.path(ROOT, "data/gis"), recursive = TRUE, full.names = TRUE)) > 0
  if (!have_gis && opt$download_boundaries) {
    helper <- file.path(script_dir, "scripts", "download_boundaries.R")
    rscript <- Sys.which("Rscript")
    if (!nzchar(rscript)) rscript <- file.path(R.home("bin"), "Rscript.exe")
    system(sprintf('"%s" "%s" "%s"', rscript, helper, ROOT))
    have_gis <- length(list.files(file.path(ROOT, "data/gis"), recursive = TRUE, full.names = TRUE)) > 0
  }
  if (!have_gis) {
    msg("Skipping maps: boundary files are not available under %s", file.path(ROOT, "data/gis"))
  } else {
    region_sf <- build_region10_polygons(ROOT)
    saveRDS(region_sf, file.path(ROOT, "results/trees/region_polygons.rds"))
    plot_region_map(region_sf, file.path(ROOT, "results/figures/fig_region_map.png"))
    cent <- region_centroids(region_sf)
    write_csv(cent, file.path(ROOT, "results/tables/region_centroids.csv"))
    plot_flow_map(region_sf, cent, summ$mean, min_mean = 0.1, out_png = file.path(ROOT, "results/figures/fig_flow_map.png"))
  }
}

msg("Done. Outputs are in: %s", file.path(ROOT, "results"))
