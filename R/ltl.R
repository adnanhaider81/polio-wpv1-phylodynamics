# R/ltl.R
suppressPackageStartupMessages({
  library(ape)
  library(phytools)
  library(dplyr)
  library(lubridate)
  library(igraph)
})

# Helper: coerce names like "Node123" or "n123" to integer node ids
parse_node_ids <- function(x) {
  if (is.null(x)) return(integer(0))
  x <- as.character(x)
  x <- gsub("[^0-9]+", "", x)
  suppressWarnings(as.integer(x))
}

# Build local transmission lineages (LTLs) from a single simmap tree.
# We cut edges where the inferred discrete state changes between parent and child node.
# IMPORTANT: some phytools versions do not store simmap_tree$states; we derive node states from:
# 1) phytools::getStates (if available)
# 2) describe.simmap()$ace (posterior probs at nodes), taking max-prob state
extract_lt_lines <- function(simmap_tree, meta, censor_months = 6) {
  ntip <- Ntip(simmap_tree)
  nnode <- simmap_tree$Nnode
  total_nodes <- ntip + nnode

  # Ensure metadata is aligned to tips
  meta <- meta[match(simmap_tree$tip.label, meta$strain), ]

  # Tip states: use observed region labels from metadata (robust)
  tip_states <- as.character(meta$region)
  if (length(tip_states) != ntip) stop("Tip states length mismatch with number of tips.")

  # Internal node states: try getStates first, then describe.simmap ace
  node_states_internal <- rep(NA_character_, nnode)
  internal_node_nums <- (ntip + 1):total_nodes

  got_internal <- FALSE

  # 1) getStates
  gs_ok <- FALSE
  gs_nodes <- NULL
  try({
    if ("getStates" %in% getNamespaceExports("phytools")) {
      gs_nodes <- phytools::getStates(simmap_tree, "nodes")
      gs_ok <- TRUE
    }
  }, silent = TRUE)

  if (gs_ok && !is.null(gs_nodes)) {
    # gs_nodes may be named; try to align by node numbers
    if (!is.null(names(gs_nodes))) {
      ids <- parse_node_ids(names(gs_nodes))
      keep <- !is.na(ids) & ids %in% internal_node_nums
      if (any(keep)) {
        node_states_internal[match(ids[keep], internal_node_nums)] <- as.character(gs_nodes[keep])
        got_internal <- TRUE
      }
    } else if (length(gs_nodes) == nnode) {
      node_states_internal <- as.character(gs_nodes)
      got_internal <- TRUE
    }
  }

  # 2) describe.simmap ace
  if (!got_internal) {
    desc <- describe.simmap(simmap_tree)
    ace <- desc$ace
    if (!is.null(ace) && is.matrix(ace) && nrow(ace) > 0) {
      # ace rows correspond to internal nodes
      st <- apply(ace, 1, function(r) colnames(ace)[which.max(r)])
      if (!is.null(rownames(ace))) {
        ids <- parse_node_ids(rownames(ace))
        keep <- !is.na(ids) & ids %in% internal_node_nums
        if (any(keep)) {
          node_states_internal[match(ids[keep], internal_node_nums)] <- st[keep]
          got_internal <- TRUE
        }
      }
      if (!got_internal && length(st) == nnode) {
        node_states_internal <- st
        got_internal <- TRUE
      }
    }
  }

  if (!got_internal) {
    stop("Could not derive internal node states from simmap tree. Try updating phytools or rerunning simmap.")
  }

  # Build full node state vector indexed by node number
  node_states <- rep(NA_character_, total_nodes)
  node_states[1:ntip] <- tip_states
  node_states[internal_node_nums] <- node_states_internal

  # Cut edges where state changes, including within-edge stochastic-map transitions.
  edge <- simmap_tree$edge
  parent_state <- node_states[edge[,1]]
  child_state  <- node_states[edge[,2]]
  cuts <- which(parent_state != child_state)
  if (!is.null(simmap_tree$maps)) {
    within_edge_cuts <- which(vapply(simmap_tree$maps, function(edge_map) {
      states <- unique(names(edge_map))
      length(states) > 1
    }, logical(1)))
    cuts <- sort(unique(c(cuts, within_edge_cuts)))
  }

  keep_idx <- setdiff(seq_len(nrow(edge)), cuts)
  keep_edges <- edge[keep_idx, , drop = FALSE]

  # Graph components for LTLs, including isolated tips
  vnames <- as.character(seq_len(total_nodes))
  if (nrow(keep_edges) == 0) {
    g <- make_empty_graph(n = total_nodes, directed = TRUE)
    V(g)$name <- vnames
  } else {
    keep_char <- apply(keep_edges, 2, as.character)
    g <- graph_from_edgelist(keep_char, directed = TRUE)
    missing_v <- setdiff(vnames, V(g)$name)
    if (length(missing_v) > 0) g <- add_vertices(g, nv = length(missing_v), name = missing_v)
  }

  comps <- components(as_undirected(g), mode = "weak")
  memb <- comps$membership
  if (is.null(names(memb))) names(memb) <- V(g)$name

  ltl_id <- memb[as.character(seq_len(ntip))]
  if (length(ltl_id) != ntip) ltl_id <- seq_len(ntip)

  df <- data.frame(
    strain = simmap_tree$tip.label,
    ltl = as.integer(ltl_id),
    region = tip_states,
    sample_date = meta$onsetdate,
    sample_type = meta$sample_type,
    stringsAsFactors = FALSE
  )

  ltl_sum <- df %>%
    group_by(ltl, region) %>%
    summarise(
      n = n(),
      first_date = min(sample_date, na.rm = TRUE),
      last_date  = max(sample_date, na.rm = TRUE),
      duration_days = as.numeric(last_date - first_date),
      n_es = sum(sample_type == "ES"),
      n_afp = sum(sample_type == "AFP"),
      .groups = "drop"
    )

  max_date <- max(df$sample_date, na.rm = TRUE)
  ltl_sum <- ltl_sum %>%
    mutate(
      censored = last_date > (max_date %m-% months(censor_months))
    )

  list(tips = df, ltl = ltl_sum, max_date = max_date)
}

classify_lt_lines <- function(ltl_sum,
                              dead_end_max_months = 6,
                              dead_end_max_n = 5,
                              persistent_min_months = 12,
                              persistent_min_n = 10) {
  ltl_sum %>%
    mutate(
      duration_months = duration_days / 30.44,
      class = dplyr::case_when(
        censored ~ "Unclassified_censored",
        duration_months < dead_end_max_months & n < dead_end_max_n ~ "Dead_end",
        duration_months >= persistent_min_months & n >= persistent_min_n ~ "Persistent",
        TRUE ~ "Other"
      )
    )
}

summarize_cluster_district_status <- function(meta,
                                              sustain_months = 6,
                                              recent_months = 3) {
  req <- c("country", "division_norm", "locality", "region", "cluster", "onsetdate", "sample_type")
  miss <- setdiff(req, names(meta))
  if (length(miss) > 0) {
    stop("Metadata is missing required columns for cluster status summary: ", paste(miss, collapse = ", "))
  }

  df <- meta %>%
    transmute(
      country = as.character(country),
      division = as.character(division_norm),
      district = toupper(trimws(as.character(locality))),
      region = as.character(region),
      cluster = as.character(cluster),
      onsetdate = as.Date(onsetdate),
      sample_type = as.character(sample_type)
    ) %>%
    filter(
      !is.na(onsetdate),
      !is.na(cluster) & nzchar(cluster),
      !is.na(district) & nzchar(district)
    )

  if (nrow(df) == 0) {
    return(data.frame(
      country = character(0),
      division = character(0),
      district = character(0),
      region = character(0),
      cluster = character(0),
      n = integer(0),
      first_date = as.Date(character(0)),
      last_date = as.Date(character(0)),
      duration_days = numeric(0),
      duration_months = numeric(0),
      n_es = integer(0),
      n_afp = integer(0),
      recent = logical(0),
      established = logical(0),
      status = character(0),
      status_label = character(0),
      stringsAsFactors = FALSE
    ))
  }

  max_date <- max(df$onsetdate, na.rm = TRUE)
  recent_cutoff <- max_date %m-% months(recent_months)

  df %>%
    group_by(country, division, district, region, cluster) %>%
    summarise(
      n = n(),
      first_date = min(onsetdate, na.rm = TRUE),
      last_date = max(onsetdate, na.rm = TRUE),
      duration_days = as.numeric(last_date - first_date),
      n_es = sum(sample_type == "ES", na.rm = TRUE),
      n_afp = sum(sample_type == "AFP", na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      duration_months = duration_days / 30.44,
      recent = last_date >= recent_cutoff,
      established = duration_months >= sustain_months,
      status = case_when(
        established & recent ~ "Established_active",
        established & !recent ~ "Established_quiet",
        !established & recent ~ "Recent_under_6m",
        TRUE ~ "Died_off_pre_established"
      ),
      status_label = recode(
        status,
        Established_active = "Established active",
        Established_quiet = "Established quiet",
        Recent_under_6m = "Recent <6 months",
        Died_off_pre_established = "Died off before establishment"
      )
    ) %>%
    arrange(desc(established), desc(recent), desc(n), cluster, district)
}

summarize_cluster_status_by_region <- function(cluster_status) {
  if (nrow(cluster_status) == 0) {
    return(data.frame(
      region = character(0),
      status = character(0),
      status_label = character(0),
      n = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  cluster_status %>%
    count(region, status, status_label, name = "n") %>%
    arrange(region, status)
}

ltl_replicate_summary <- function(simmap_list, meta, censor_months = 6) {
  out <- lapply(seq_along(simmap_list), function(i) {
    obj <- extract_lt_lines(simmap_list[[i]], meta, censor_months = censor_months)
    classified <- classify_lt_lines(obj$ltl)
    classified %>%
      count(class, name = "n_lineages") %>%
      mutate(replicate = i, .before = 1)
  })
  dplyr::bind_rows(out)
}
