suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

read_metadata_table <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(readr::read_csv(path, show_col_types = FALSE))
  if (ext %in% c("tsv", "txt")) return(readr::read_tsv(path, show_col_types = FALSE))
  if (ext %in% c("xlsx", "xls")) {
    if (!requireNamespace("readxl", quietly = TRUE)) {
      stop("Install readxl to read Excel metadata: install.packages('readxl')", call. = FALSE)
    }
    return(readxl::read_excel(path))
  }
  stop("Unsupported metadata file extension: ", ext, call. = FALSE)
}

parse_is_env <- function(x) {
  if (is.logical(x)) return(ifelse(is.na(x), FALSE, x))
  y <- toupper(trimws(as.character(x)))
  y %in% c("TRUE", "T", "YES", "Y", "1", "ES", "ENV", "ENVIRONMENTAL")
}

read_metadata <- function(metadata_path) {
  meta <- read_metadata_table(metadata_path)
  req <- c("FILENAME", "Country", "StateProv", "locality", "onsetdate", "IsENV", "Cluster")
  miss <- setdiff(req, names(meta))
  if (length(miss) > 0) {
    stop("Metadata is missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  meta %>%
    mutate(
      strain = str_trim(as.character(FILENAME)),
      country = toupper(str_trim(as.character(Country))),
      StateProv = str_trim(as.character(StateProv)),
      locality = str_trim(as.character(locality)),
      onsetdate = as.Date(onsetdate),
      decdate = vapply(onsetdate, decimal_date, numeric(1)),
      sample_type = ifelse(parse_is_env(IsENV), "ES", "AFP"),
      cluster = as.character(Cluster)
    ) %>%
    filter(!is.na(strain), nzchar(strain))
}

downsample_by_region_month <- function(meta, max_per_region_month = 120, seed = 1) {
  set.seed(seed)
  meta %>%
    mutate(month = format(onsetdate, "%Y-%m")) %>%
    group_by(region, month) %>%
    slice_sample(n = min(n(), max_per_region_month)) %>%
    ungroup() %>%
    select(-month)
}

write_nextstrain_metadata <- function(meta, out_tsv) {
  ns <- meta %>%
    transmute(
      strain = strain,
      date = format(onsetdate, "%Y-%m-%d"),
      country = country,
      division = division_norm,
      location = as.character(locality),
      region = region,
      is_env = sample_type == "ES",
      cluster = cluster
    )
  readr::write_tsv(ns, out_tsv)
  invisible(ns)
}

read_fasta_headers <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sub("^>", "", lines[grepl("^>", lines)])
}

validate_alignment_metadata <- function(fasta_path, meta) {
  headers <- read_fasta_headers(fasta_path)
  if (length(headers) == 0) stop("No FASTA records found in: ", fasta_path, call. = FALSE)
  dup_headers <- headers[duplicated(headers)]
  if (length(dup_headers) > 0) {
    stop("Duplicate FASTA headers found: ", paste(unique(dup_headers), collapse = ", "), call. = FALSE)
  }
  missing_meta <- setdiff(headers, meta$strain)
  missing_fasta <- setdiff(meta$strain, headers)
  if (length(missing_meta) > 0) {
    stop("FASTA records missing from metadata: ", paste(head(missing_meta, 20), collapse = ", "), call. = FALSE)
  }
  if (length(missing_fasta) > 0) {
    warning("Metadata rows without FASTA records will be ignored: ", paste(head(missing_fasta, 20), collapse = ", "))
  }
  invisible(TRUE)
}
