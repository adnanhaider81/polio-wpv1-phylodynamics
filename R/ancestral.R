suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
})

run_simmap <- function(tree, tip_states, nsim = 25, model = "ARD", seed = 1) {
  set.seed(seed)
  tip_states <- tip_states[tree$tip.label]
  if (any(is.na(tip_states))) stop("Missing tip states for some tips in the tree.", call. = FALSE)
  phytools::make.simmap(tree, tip_states, model = model, nsim = nsim, message = TRUE)
}

transition_counts <- function(simmap_list, tip_states) {
  if (length(simmap_list) == 0) stop("Empty simmap list.", call. = FALSE)
  map_states <- unique(unlist(lapply(simmap_list, function(sm) unlist(lapply(sm$maps, names)))))
  all_states <- sort(unique(c(as.character(tip_states), as.character(map_states))))
  all_states <- all_states[!is.na(all_states) & nzchar(all_states)]
  if (length(all_states) < 2) stop("Need at least two states for transition counts.", call. = FALSE)

  lapply(simmap_list, function(sm) {
    m <- matrix(0, nrow = length(all_states), ncol = length(all_states),
                dimnames = list(all_states, all_states))
    for (edge_map in sm$maps) {
      states <- names(edge_map)
      if (length(states) < 2) next
      for (i in 2:length(states)) {
        from <- states[[i - 1]]
        to <- states[[i]]
        if (!is.na(from) && !is.na(to) && nzchar(from) && nzchar(to) && from != to) {
          m[from, to] <- m[from, to] + 1
        }
      }
    }
    m
  })
}

stack_transition_mats <- function(mats) {
  if (length(mats) == 0) stop("No transition matrices to stack.", call. = FALSE)
  n <- nrow(mats[[1]])
  k <- length(mats)
  arr <- array(0, dim = c(n, n, k),
               dimnames = list(rownames(mats[[1]]), colnames(mats[[1]]), NULL))
  for (i in seq_len(k)) arr[, , i] <- mats[[i]]
  arr
}

summarize_transition_matrix <- function(mats) {
  arr <- stack_transition_mats(mats)
  list(
    mean = apply(arr, c(1, 2), mean),
    lo = apply(arr, c(1, 2), quantile, probs = 0.025, na.rm = TRUE),
    hi = apply(arr, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
  )
}

transition_support <- function(mats) {
  arr <- stack_transition_mats(mats)
  apply(arr, c(1, 2), function(x) mean(x > 0))
}

top_corridors <- function(mat_mean, support_mat = NULL) {
  diag(mat_mean) <- 0
  out <- as.data.frame(as.table(mat_mean), stringsAsFactors = FALSE) %>%
    rename(from = Var1, to = Var2, mean_transitions = Freq) %>%
    filter(from != to, mean_transitions > 0) %>%
    arrange(desc(mean_transitions))
  if (!is.null(support_mat)) {
    sup <- as.data.frame(as.table(support_mat), stringsAsFactors = FALSE) %>%
      rename(from = Var1, to = Var2, support = Freq)
    out <- left_join(out, sup, by = c("from", "to"))
  }
  out
}

flow_partition <- function(mat_mean) {
  states <- rownames(mat_mean)
  diag(mat_mean) <- 0
  pak <- grepl("^PAK_", states)
  afg <- grepl("^AFG_", states)
  data.frame(
    partition = c("within_pak", "within_afg", "pak_to_afg", "afg_to_pak", "other"),
    mean_transitions = c(
      sum(mat_mean[pak, pak, drop = FALSE]),
      sum(mat_mean[afg, afg, drop = FALSE]),
      sum(mat_mean[pak, afg, drop = FALSE]),
      sum(mat_mean[afg, pak, drop = FALSE]),
      sum(mat_mean[!(pak | afg), , drop = FALSE]) + sum(mat_mean[, !(pak | afg), drop = FALSE])
    )
  )
}

import_export_from_matrix <- function(mat) {
  states <- rownames(mat)
  imports <- sapply(states, function(s) sum(mat[setdiff(states, s), s]))
  exports <- sapply(states, function(s) sum(mat[s, setdiff(states, s)]))
  data.frame(region = states, imports = imports, exports = exports, stringsAsFactors = FALSE)
}
