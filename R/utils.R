suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(lubridate)
})

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  invisible(TRUE)
}

dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

as_date_safe <- function(x) {
  if (inherits(x, "Date")) return(x)
  suppressWarnings(as.Date(x))
}

decimal_date <- function(d) {
  d <- as_date_safe(d)
  y <- lubridate::year(d)
  start <- as.Date(paste0(y, "-01-01"))
  end <- as.Date(paste0(y + 1, "-01-01"))
  y + as.numeric(d - start) / as.numeric(end - start)
}

write_tsv <- function(df, path) {
  readr::write_tsv(df, path)
}

write_csv <- function(df, path) {
  readr::write_csv(df, path)
}

write_matrix_csv <- function(mat, path) {
  out <- data.frame(from = rownames(mat), mat, check.names = FALSE)
  readr::write_csv(out, path)
}
