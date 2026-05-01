suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ape)
  library(scales)
  library(tidyr)
})

.pretty_region_labels <- c(
  AFG_HELMAND = "AFG Helmand",
  AFG_KANDAHAR = "AFG Kandahar",
  AFG_OTHER = "AFG Other",
  PAK_BALOCHISTAN = "PAK Balochistan",
  PAK_ISLAMABAD = "PAK Islamabad",
  PAK_KARACHI = "PAK Karachi",
  PAK_KP = "PAK KP",
  PAK_NORTH = "PAK North",
  PAK_PUNJAB = "PAK Punjab",
  PAK_SINDH = "PAK Sindh"
)

pretty_region_label <- function(x) {
  x <- as.character(x)
  out <- .pretty_region_labels[x]
  miss <- is.na(out)
  out[miss] <- gsub("_", " ", x[miss])
  unname(out)
}

plot_tree <- function(tree, meta, out_png, width = 10, height = 7) {
  meta <- meta[match(tree$tip.label, meta$strain), ]
  if (requireNamespace("ggtree", quietly = TRUE)) {
    suppressPackageStartupMessages(library(ggtree))
    df <- meta %>% select(strain, region, sample_type)
    p <- ggtree(tree, linewidth = 0.3) %<+% df +
      geom_tippoint(aes(color = region, shape = sample_type), size = 1.2, alpha = 0.9) +
      theme_tree2() +
      theme(legend.position = "right")
    ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
    return(p)
  }

  png(out_png, width = width, height = height, units = "in", res = 300, bg = "white")
  par(mar = c(3, 1, 2, 1))
  ape::plot.phylo(tree, show.tip.label = FALSE, cex = 0.45, no.margin = TRUE)
  title("Time-scaled WPV1 tree")
  dev.off()
  invisible(TRUE)
}

plot_ltt <- function(tree, censor_days = 180, out_png, width = 8, height = 4) {
  bt <- ape::branching.times(tree)
  times <- sort(bt)
  n <- length(tree$tip.label)
  ltt <- data.frame(
    time = c(0, times, max(ape::node.depth.edgelength(tree))),
    lineages = c(1, seq_len(length(times)) + 1, n)
  )
  max_t <- max(ltt$time)
  ltt <- ltt %>% filter(time <= (max_t - censor_days / 365.25))
  p <- ggplot(ltt, aes(x = time, y = lineages)) +
    geom_step(linewidth = 0.7) +
    labs(x = "Time from root (years)", y = "Lineages") +
    theme_bw()
  ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
  p
}

plot_import_export <- function(summ_ie, out_png, width = 10, height = 4) {
  long <- summ_ie %>%
    pivot_longer(cols = c(imports_mean, exports_mean), names_to = "type", values_to = "mean") %>%
    mutate(
      lo = ifelse(type == "imports_mean", imports_lo, exports_lo),
      hi = ifelse(type == "imports_mean", imports_hi, exports_hi),
      type = ifelse(type == "imports_mean", "Importations", "Exportations")
    )
  p <- ggplot(long, aes(x = reorder(pretty_region_label(region), mean), y = mean, fill = type)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(width = 0.8), width = 0.2) +
    coord_flip() +
    labs(x = NULL, y = "Estimated transition count", fill = NULL) +
    theme_bw()
  ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
  p
}

plot_transition_heatmap <- function(mat_mean, out_png, width = 8, height = 7) {
  diag(mat_mean) <- 0
  df <- as.data.frame(as.table(mat_mean), stringsAsFactors = FALSE) %>%
    rename(from = Var1, to = Var2, mean = Freq) %>%
    filter(from != to)
  p <- ggplot(df, aes(x = pretty_region_label(to), y = pretty_region_label(from), fill = mean)) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_viridis_c(option = "C", direction = -1) +
    labs(x = "Destination", y = "Origin", fill = "Mean") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
  p
}

plot_chord <- function(mat_mean, support_mat = NULL, min_mean = 2, max_links = 24, out_png, width = 9, height = 9) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    return(plot_transition_heatmap(mat_mean, out_png, width = width, height = height))
  }
  suppressPackageStartupMessages(library(circlize))

  m <- mat_mean
  diag(m) <- 0
  m[m < min_mean] <- 0
  if (!is.null(support_mat)) m[support_mat < 0.5] <- 0

  links <- as.data.frame(as.table(m), stringsAsFactors = FALSE) %>%
    rename(from = Var1, to = Var2, weight = Freq) %>%
    filter(from != to, weight > 0) %>%
    arrange(desc(weight))

  if (nrow(links) == 0) {
    return(plot_transition_heatmap(mat_mean, out_png, width = width, height = height))
  }
  if (!is.null(max_links) && nrow(links) > max_links) links <- links %>% slice_head(n = max_links)

  regions <- rownames(m)[rownames(m) %in% unique(c(links$from, links$to))]
  region_cols <- c(
    AFG_HELMAND = "#e76f51",
    AFG_KANDAHAR = "#f4a261",
    AFG_OTHER = "#e9c46a",
    PAK_BALOCHISTAN = "#264653",
    PAK_ISLAMABAD = "#2a9d8f",
    PAK_KARACHI = "#1d3557",
    PAK_KP = "#457b9d",
    PAK_NORTH = "#4d908e",
    PAK_PUNJAB = "#90be6d",
    PAK_SINDH = "#f94144"
  )
  miss_cols <- setdiff(regions, names(region_cols))
  if (length(miss_cols) > 0) {
    extra <- grDevices::hcl(h = seq(15, 375, length.out = length(miss_cols) + 1)[-1], c = 80, l = 55)
    region_cols <- c(region_cols, setNames(extra, miss_cols))
  }
  region_cols <- region_cols[regions]

  png(out_png, width = width, height = height, units = "in", res = 300, bg = "white")
  par(mar = c(1, 1, 3, 1), bg = "white")
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.after = rep(5, length(regions)))
  circlize::chordDiagramFromDataFrame(
    links %>% select(from, to, weight),
    order = regions,
    grid.col = region_cols,
    directional = 1,
    direction.type = "arrows",
    link.arr.type = "big.arrow",
    annotationTrack = "grid",
    transparency = 0.35
  )
  title(main = sprintf("Inferred movement corridors (top %d directed links)", nrow(links)), cex.main = 1.05)
  circlize::circos.clear()
  dev.off()
  invisible(TRUE)
}

plot_ltl_facets <- function(ltl_classified, out_png, width = 11, height = 7) {
  df <- ltl_classified %>%
    filter(class != "Unclassified_censored") %>%
    mutate(
      class = recode(
        class,
        Dead_end = "Dead-end lineage",
        Persistent = "Persistent lineage",
        Other = "Other lineage"
      ),
      class = factor(class, levels = c("Dead-end lineage", "Persistent lineage", "Other lineage")),
      month = as.Date(format(first_date, "%Y-%m-01"))
    ) %>%
    count(region, month, class)
  if (nrow(df) == 0) {
    png(out_png, width = width, height = height, units = "in", res = 300, bg = "white")
    plot.new()
    text(0.5, 0.5, "No classified lineages to plot.")
    dev.off()
    return(invisible(FALSE))
  }
  p <- ggplot(df, aes(x = month, y = n, fill = class)) +
    geom_col() +
    facet_wrap(~region, scales = "free_y", ncol = 4, labeller = as_labeller(pretty_region_label)) +
    scale_fill_manual(
      values = c(
        "Dead-end lineage" = "#d95f02",
        "Persistent lineage" = "#1b9e77",
        "Other lineage" = "#6c757d"
      ),
      drop = FALSE
    ) +
    labs(x = NULL, y = "Number of lineages", fill = "Lineage class") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
  p
}

plot_program_cluster_status <- function(cluster_status,
                                        out_png,
                                        sustain_months = 6,
                                        recent_months = 3,
                                        width = 10,
                                        height = 5.5) {
  if (nrow(cluster_status) == 0) {
    png(out_png, width = width, height = height, units = "in", res = 300, bg = "white")
    plot.new()
    text(0.5, 0.5, "No cluster status records available.")
    dev.off()
    return(invisible(FALSE))
  }
  status_levels <- c(
    "Established active",
    "Established quiet",
    "Recent <6 months",
    "Died off before establishment"
  )
  df <- cluster_status %>%
    mutate(status_label = factor(status_label, levels = status_levels)) %>%
    count(region, status_label, name = "n") %>%
    group_by(region) %>%
    mutate(total = sum(n)) %>%
    ungroup()
  p <- ggplot(df, aes(x = reorder(pretty_region_label(region), total), y = n, fill = status_label)) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Established active" = "#d62828",
        "Established quiet" = "#6c757d",
        "Recent <6 months" = "#fcbf49",
        "Died off before establishment" = "#2a9d8f"
      ),
      drop = FALSE
    ) +
    labs(
      x = NULL,
      y = "Cluster-location records",
      fill = "Status",
      subtitle = paste0("Established >= ", sustain_months, " months; recent window = ", recent_months, " months.")
    ) +
    theme_bw(base_size = 11)
  ggsave(out_png, p, width = width, height = height, dpi = 300, bg = "white")
  p
}
