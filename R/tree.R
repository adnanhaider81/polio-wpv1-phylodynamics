suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(phangorn)
  library(treedater)
})

get_alignment_length <- function(aln_fasta) {
  dna <- ape::read.dna(aln_fasta, format = "fasta")
  ncol(as.matrix(dna))
}

run_distance_tree <- function(aln_fasta, out_tree, model = "raw") {
  dna <- ape::read.dna(aln_fasta, format = "fasta")
  d <- ape::dist.dna(dna, model = model, pairwise.deletion = TRUE, as.matrix = FALSE)
  tr <- ape::nj(d)
  tr <- tryCatch(phangorn::midpoint(tr), error = function(e) tr)
  ape::write.tree(tr, out_tree)
  out_tree
}

run_fasttree <- function(aln_fasta, out_tree) {
  exe <- Sys.which(c("FastTree", "fasttree"))
  exe <- exe[nzchar(exe)][[1]]
  cmd <- sprintf('%s -nt %s > %s', shQuote(exe), shQuote(aln_fasta), shQuote(out_tree))
  message("Running: ", cmd)
  status <- system(cmd)
  if (status != 0) stop("FastTree failed.", call. = FALSE)
  out_tree
}

run_iqtree <- function(aln_fasta, out_prefix, threads = "AUTO", model = "HKY+G", ufboot = 1000) {
  env_ufboot <- suppressWarnings(as.integer(Sys.getenv("POLIO_UFBOOT", "")))
  if (!is.na(env_ufboot) && env_ufboot >= 1000) ufboot <- env_ufboot
  exe <- Sys.which(c("iqtree2", "iqtree3", "iqtree"))
  exe <- exe[nzchar(exe)][[1]]
  cmd <- sprintf('%s -s %s -m %s -B %s -nt %s -pre %s',
                 shQuote(exe), shQuote(aln_fasta), model, ufboot, threads, shQuote(out_prefix))
  message("Running: ", cmd)
  status <- system(cmd)
  if (status != 0) stop("IQ-TREE failed.", call. = FALSE)
  paste0(out_prefix, ".treefile")
}

time_scale_tree <- function(tree_path,
                            meta,
                            seq_len,
                            omega0 = 0.01,
                            clock = "strict",
                            ncpu = 1) {
  tr <- ape::read.tree(tree_path)
  meta2 <- meta %>%
    filter(strain %in% tr$tip.label) %>%
    filter(!is.na(decdate))
  keep <- tr$tip.label %in% meta2$strain
  tr <- ape::drop.tip(tr, tr$tip.label[!keep])
  meta2 <- meta2[match(tr$tip.label, meta2$strain), ]
  sts <- setNames(meta2$decdate, meta2$strain)
  td <- treedater::dater(
    tre = tr,
    sts = sts,
    s = seq_len,
    omega0 = omega0,
    clock = clock,
    ncpu = ncpu,
    quiet = TRUE
  )
  list(tree = td, treedater = td, meta = meta2)
}

root_to_tip_diagnostics <- function(tree, meta) {
  meta <- meta[match(tree$tip.label, meta$strain), ]
  depths <- ape::node.depth.edgelength(tree)[seq_along(tree$tip.label)]
  df <- data.frame(
    strain = tree$tip.label,
    decdate = meta$decdate,
    root_to_tip = depths,
    region = meta$region,
    stringsAsFactors = FALSE
  )
  fit <- stats::lm(root_to_tip ~ decdate, data = df)
  df$slope <- unname(stats::coef(fit)[["decdate"]])
  df$r_squared <- summary(fit)$r.squared
  df
}
