# R/regions.R
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

normalize_state <- function(x) {
  x <- toupper(str_trim(as.character(x)))
  x <- str_replace_all(x, "\\s+", " ")
  if (x %in% c("KPK","K. PAKHTUNKHWA","KHYBER-PAKHTUNKHWA","KHYBER PAKHTOONKHW","KHYBER PAKHTUNKHW","KHYBER PAKHTUNKHWA")) return("KHYBER PAKHTUNKHWA")
  if (x %in% c("GILGIT BALTISTAN","GBALTISTAN","GB")) return("GILGIT-BALTISTAN")
  if (x %in% c("AZAD JAMMU AND KASHMIR","AJK")) return("AJK")
  x
}

# Proof-of-concept region scheme (10 states)
# This is intentionally coarse because 2025-only data and uneven sampling can make fine spatial DTA unstable.
assign_region10 <- function(country, division, locality) {
  c <- toupper(as.character(country))
  div <- normalize_state(division)
  loc <- toupper(as.character(locality))
  if (c == "PAK") {
    if (div == "SINDH") {
      if (str_starts(loc, "KHI")) return("PAK_KARACHI")
      return("PAK_SINDH")
    }
    if (div == "KHYBER PAKHTUNKHWA") return("PAK_KP")
    if (div == "PUNJAB") return("PAK_PUNJAB")
    if (div == "BALOCHISTAN") return("PAK_BALOCHISTAN")
    if (div == "ISLAMABAD") return("PAK_ISLAMABAD")
    if (div %in% c("GILGIT-BALTISTAN","AJK")) return("PAK_NORTH")
    return(paste0("PAK_", div))
  }
  if (c == "AFG") {
    if (div == "HILMAND") return("AFG_HELMAND")
    if (div == "KANDAHAR") return("AFG_KANDAHAR")
    return("AFG_OTHER")
  }
  paste0(c, "_", div)
}

apply_region_scheme <- function(meta, scheme = "region10") {
  if (scheme != "region10") stop("Only region10 is implemented in this template. Add more in R/regions.R.")
  meta %>%
    mutate(
      division_norm = vapply(StateProv, normalize_state, character(1)),
      region = mapply(assign_region10, Country, division_norm, locality) %>% as.character()
    )
}
