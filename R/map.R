# R/map.R
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(stringr)
  library(geosphere)
})

.map_region_labels <- c(
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

map_region_label <- function(x) {
  x <- as.character(x)
  out <- .map_region_labels[x]
  miss <- is.na(out)
  out[miss] <- gsub("_", " ", x[miss])
  unname(out)
}

get_gadm <- function(root, iso3, level) {
  dir.create(file.path(root, "data/gis"), recursive = TRUE, showWarnings = FALSE)
  suppressPackageStartupMessages(library(geodata))
  f <- geodata::gadm(country = iso3, level = level, path = file.path(root, "data/gis"))
  st_as_sf(f)
}

build_region10_polygons <- function(root) {
  pak1 <- get_gadm(root, "PAK", 1) %>%
    transmute(admin0="PAK", name1=NAME_1, geom=geometry) %>%
    st_as_sf()
  pak2 <- get_gadm(root, "PAK", 2) %>%
    transmute(admin0="PAK", name1=NAME_1, name2=NAME_2, geom=geometry) %>%
    st_as_sf()
  afg1 <- get_gadm(root, "AFG", 1) %>%
    transmute(admin0="AFG", name1=NAME_1, geom=geometry) %>%
    st_as_sf()

  karachi <- pak2 %>%
    filter(str_to_upper(name1) == "SINDH") %>%
    filter(str_detect(str_to_upper(name2), "KARACHI"))

  sindh_rest <- pak2 %>%
    filter(str_to_upper(name1) == "SINDH") %>%
    filter(!str_detect(str_to_upper(name2), "KARACHI")) %>%
    summarise(admin0="PAK", region="PAK_SINDH", geometry=st_union(geom), .groups="drop") %>%
    st_as_sf()

  karachi_poly <- karachi %>%
    summarise(admin0="PAK", region="PAK_KARACHI", geometry=st_union(geom), .groups="drop") %>%
    st_as_sf()

  pak_other <- pak1 %>%
    mutate(name1_upper = str_to_upper(name1)) %>%
    mutate(region = case_when(
      (str_detect(name1_upper, "KHYBER") & str_detect(name1_upper, "PAKHT")) ~ "PAK_KP",
      str_detect(name1_upper, "FEDERALLY ADMINISTERED TRIBAL") ~ "PAK_KP",
      name1_upper == "PUNJAB" ~ "PAK_PUNJAB",
      name1_upper == "BALOCHISTAN" ~ "PAK_BALOCHISTAN",
      name1_upper == "ISLAMABAD" ~ "PAK_ISLAMABAD",
      name1_upper %in% c("GILGIT-BALTISTAN","AZAD KASHMIR") ~ "PAK_NORTH",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(region)) %>%
    group_by(admin0, region) %>%
    summarise(geometry=st_union(geom), .groups="drop") %>%
    st_as_sf()

  pak_regions <- bind_rows(karachi_poly, sindh_rest, pak_other)

  afg_regions <- afg1 %>%
    mutate(region = case_when(
      str_to_upper(name1) %in% c("HELMAND","HILMAND") ~ "AFG_HELMAND",
      str_to_upper(name1) == "KANDAHAR" ~ "AFG_KANDAHAR",
      TRUE ~ "AFG_OTHER"
    )) %>%
    group_by(admin0, region) %>%
    summarise(geometry=st_union(geom), .groups="drop") %>%
    st_as_sf()

  bind_rows(pak_regions, afg_regions)
}

# Compute stable centroids: project to EPSG:3857 first, then back to lon/lat.
region_centroids <- function(region_sf) {
  region_ll <- st_make_valid(region_sf) %>% st_transform(4326)
  region_proj <- st_transform(region_ll, 3857)
  cent_proj <- st_centroid(region_proj)
  cent_ll <- st_transform(cent_proj, 4326)

  coords <- st_coordinates(cent_ll)
  data.frame(
    region = region_ll$region,
    lon = coords[,1],
    lat = coords[,2],
    stringsAsFactors = FALSE
  )
}

plot_region_map <- function(region_sf, out_png, width=7, height=7) {
  p <- ggplot(region_sf) +
    geom_sf(aes(fill=region), color="grey35", linewidth=0.2, alpha=0.85) +
    coord_sf() +
    theme_void() +
    theme(legend.position="right")
  ggsave(out_png, p, width=width, height=height, dpi=300)
  p
}

plot_sample_points <- function(region_sf, meta, out_png, width=7, height=7) {
  pts <- meta %>% filter(!is.na(latitude) & !is.na(longitude))
  if (nrow(pts) == 0) stop("No lat/long found in metadata. Cannot plot sample points.")
  sf_pts <- st_as_sf(pts, coords=c("longitude","latitude"), crs=4326, remove=FALSE)
  p <- ggplot() +
    geom_sf(data=region_sf, fill=NA, color="grey35", linewidth=0.2) +
    geom_sf(data=sf_pts, aes(color=region, shape=sample_type), size=1.0, alpha=0.7) +
    coord_sf() +
    theme_void() +
    theme(legend.position="right")
  ggsave(out_png, p, width=width, height=height, dpi=300)
  p
}

# Plot movement flows as great-circle arcs between region centroids.
# If no flows pass filter, create a placeholder figure (do not error).
plot_flow_map <- function(region_sf, centroids_df, mat_mean, min_mean = 2, max_corridors = 28, label_top_n = 8, out_png, width=9, height=8) {
  diag(mat_mean) <- 0
  flows <- as.data.frame(as.table(mat_mean), stringsAsFactors = FALSE) %>%
    rename(from=Var1, to=Var2, weight=Freq) %>%
    filter(from != to, weight >= min_mean) %>%
    arrange(desc(weight))

  if (nrow(flows) == 0) {
    png(out_png, width=width, height=height, units="in", res=300, bg="white")
    par(mar = c(0,0,0,0))
    plot.new()
    text(0.5, 0.5,
         labels = paste0("No flows above min_mean = ", min_mean, "\n",
                         "Lower min_mean in run_pipeline.R (plot_flow_map call)\n",
                         "or increase nsim to stabilize counts."),
         cex = 1.1)
    dev.off()
    return(invisible(FALSE))
  }

  coords <- centroids_df
  flows <- flows %>%
    left_join(coords, by=c("from"="region")) %>%
    rename(lon_from=lon, lat_from=lat) %>%
    left_join(coords, by=c("to"="region")) %>%
    rename(lon_to=lon, lat_to=lat)

  good <- with(flows, is.finite(lon_from) & is.finite(lat_from) & is.finite(lon_to) & is.finite(lat_to))
  if (any(!good)) {
    bad_regions <- sort(unique(c(flows$from[!good], flows$to[!good])))
    warning(
      "Dropping ", sum(!good), " flow(s) due to missing/invalid centroid coordinates. Regions: ",
      paste(bad_regions, collapse = ", ")
    )
    flows <- flows[good, , drop = FALSE]
  }

  if (nrow(flows) == 0) {
    png(out_png, width=width, height=height, units="in", res=300)
    par(mar = c(0,0,0,0))
    plot.new()
    text(0.5, 0.5,
         labels = "No plottable flows after centroid matching.\nCheck region naming between metadata and map polygons.",
         cex = 1.1)
    dev.off()
    return(invisible(FALSE))
  }

  if (!is.null(max_corridors) && nrow(flows) > max_corridors) {
    flows <- flows %>% slice_head(n = max_corridors)
  }

  flows <- flows %>%
    mutate(
      seg = row_number(),
      corridor_label = paste0(map_region_label(from), " -> ", map_region_label(to))
    )

  make_arc <- function(lon1, lat1, lon2, lat2, n=50) {
    if (!all(is.finite(c(lon1, lat1, lon2, lat2)))) return(NULL)
    gc <- tryCatch(
      geosphere::gcIntermediate(c(lon1,lat1), c(lon2,lat2), n=n, addStartEnd=TRUE, breakAtDateLine=FALSE),
      error = function(e) NULL
    )
    if (is.null(gc)) {
      return(data.frame(lon = c(lon1, lon2), lat = c(lat1, lat2)))
    }
    gc_df <- as.data.frame(gc)
    if (nrow(gc_df) == 0 || ncol(gc_df) < 2) {
      return(data.frame(lon = c(lon1, lon2), lat = c(lat1, lat2)))
    }
    names(gc_df)[1:2] <- c("lon", "lat")
    gc_df %>% transmute(lon=lon, lat=lat)
  }

  arcs <- lapply(seq_len(nrow(flows)), function(i) {
    f <- flows[i,]
    arc <- make_arc(f$lon_from, f$lat_from, f$lon_to, f$lat_to, n=60)
    if (is.null(arc) || nrow(arc) == 0) return(NULL)
    arc$weight <- f$weight
    arc$seg <- f$seg
    arc$from <- f$from
    arc$to <- f$to
    arc
  }) %>% bind_rows()

  if (nrow(arcs) == 0) {
    png(out_png, width=width, height=height, units="in", res=300, bg="white")
    par(mar = c(0,0,0,0))
    plot.new()
    text(0.5, 0.5,
         labels = "No arcs could be created from available flow coordinates.",
         cex = 1.1)
    dev.off()
    return(invisible(FALSE))
  }

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
  regions_used <- unique(c(as.character(centroids_df$region), as.character(flows$from), as.character(flows$to)))
  miss_cols <- setdiff(regions_used, names(region_cols))
  if (length(miss_cols) > 0) {
    extra <- grDevices::hcl(h = seq(15, 375, length.out = length(miss_cols) + 1)[-1], c = 80, l = 55)
    region_cols <- c(region_cols, setNames(extra, miss_cols))
  }
  region_cols <- region_cols[regions_used]

  nodes <- centroids_df %>%
    mutate(
      region = as.character(region),
      label = map_region_label(region),
      nudge_x = case_when(
        region == "PAK_ISLAMABAD" ~ 0.45,
        region == "PAK_KP" ~ -0.15,
        region == "PAK_PUNJAB" ~ 0.55,
        region == "PAK_NORTH" ~ 0.45,
        region == "PAK_KARACHI" ~ 0.25,
        region == "AFG_KANDAHAR" ~ -0.20,
        TRUE ~ 0.15
      ),
      nudge_y = case_when(
        region == "PAK_KARACHI" ~ -0.55,
        region == "PAK_ISLAMABAD" ~ 0.45,
        region == "PAK_KP" ~ 0.55,
        region == "PAK_PUNJAB" ~ -0.25,
        region == "PAK_NORTH" ~ 0.45,
        region == "PAK_BALOCHISTAN" ~ 0.35,
        TRUE ~ 0.20
      )
    )

  mids <- arcs %>%
    group_by(seg) %>%
    mutate(row_id = row_number(), n_row = n()) %>%
    filter(row_id == floor((n_row + 1) / 2)) %>%
    ungroup() %>%
    select(seg, lon_mid = lon, lat_mid = lat)

  label_n <- min(label_top_n, nrow(flows))
  corridor_labels <- flows %>%
    arrange(desc(weight)) %>%
    slice_head(n = label_n) %>%
    left_join(mids, by = "seg")

  p <- ggplot() +
    geom_sf(data=region_sf, aes(fill=region), color="white", linewidth=0.35, alpha=0.25, show.legend=FALSE) +
    geom_path(data=arcs, aes(x=lon, y=lat, group=seg, color=from, linewidth=weight), alpha=0.72, lineend="round") +
    geom_point(data=centroids_df, aes(x=lon, y=lat, fill=region), shape=21, color="white", stroke=0.5, size=2.8, show.legend=FALSE) +
    geom_segment(
      data=nodes,
      aes(x=lon, y=lat, xend=lon + nudge_x, yend=lat + nudge_y, color=region),
      linewidth=0.25,
      alpha=0.7,
      show.legend=FALSE
    ) +
    geom_label(
      data=nodes,
      aes(x=lon + nudge_x, y=lat + nudge_y, label=label, fill=region),
      color="black",
      size=2.9,
      linewidth=0.15,
      label.padding=grid::unit(0.12, "lines"),
      alpha=0.92,
      show.legend=FALSE
    ) +
    geom_text(
      data=corridor_labels,
      aes(x=lon_mid, y=lat_mid, label=corridor_label, color=from),
      size=2.5,
      fontface="bold",
      check_overlap=TRUE,
      alpha=0.92,
      show.legend=FALSE
    ) +
    coord_sf(crs=st_crs(4326), expand=FALSE) +
    scale_color_manual(values=region_cols, labels=map_region_label(names(region_cols))) +
    scale_fill_manual(values=grDevices::adjustcolor(region_cols, alpha.f=0.35), labels=map_region_label(names(region_cols))) +
    scale_linewidth_continuous(range=c(0.5, 2.8), breaks=scales::pretty_breaks(4)) +
    labs(
      title = "Inferred WPV1 Movement Corridors",
      subtitle = paste0(
        "Top ", nrow(flows), " corridors with mean transitions >= ", min_mean,
        "; labels show top ", label_n, " corridors."
      ),
      color = "Origin region",
      linewidth = "Mean transitions"
    ) +
    theme_minimal(base_size=11) +
    theme(
      panel.grid.major = element_line(color="#d9e2ec", linewidth=0.25),
      panel.grid.minor = element_blank(),
      legend.position="right",
      plot.title = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=9.5),
      panel.background = element_rect(fill="#f6fbff", color=NA),
      plot.background = element_rect(fill="white", color=NA)
    )

  ggsave(out_png, p, width=width, height=height, dpi=320, bg="white")
  p
}
