#!/usr/bin/env Rscript

# AIMSLTMP-CalibrationData.R
# -------------------------------------------------------------------
# R translation of the Julia script "AIMSLTMP-CalibrationData.jl".
# It pulls reef + sector monitoring data from the AIMS Data API and
# combines them into tidy data frames.
#
# Outputs created in the global environment:
#   - reef_info
#   - reef_photo_df
#   - reef_manta_df
#   - reef_disturbance_df
#   - reef_cots_df
#   - sector_photo_df
#   - sector_manta_df
#   - sector_cots_df
#
# Optional: set `out_dir` (below) to write CSVs and an RDS bundle.
# -------------------------------------------------------------------

# ---- packages ----
required_pkgs <- c("httr", "jsonlite", "dplyr", "tibble", "purrr", "readr")
missing <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop(
    "Missing required packages: ", paste(missing, collapse = ", "),
    "\nInstall them first, e.g.: install.packages(c(", paste(sprintf('"%s"', missing), collapse = ", "), "))"
  )
}

library(httr)
library(jsonlite)
library(dplyr)
library(tibble)
library(purrr)
library(readr)

# ---- base URL + endpoint builders ----
BASE_URL <- "https://api.aims.gov.au/data-v2.0/10.25845/5c09bc4ff315c/"

# simple "joinpath" equivalent for URLs
url_join <- function(base, path) {
  base <- sub("/+$", "/", base)
  path <- sub("^/+", "", path)
  paste0(base, path)
}

# URL encoding similar to HTTP.escapeuri in Julia
uenc <- function(x) utils::URLencode(x, reserved = TRUE)

get_reef_info_url <- url_join(BASE_URL, "reef")

get_photo_transect <- function(name, domain_category = "reef") {
  paste0(
    url_join(BASE_URL, "data"),
    "?domain_name=", uenc(name),
    "&domain_category=", uenc(domain_category),
    "&data_type=photo-transect"
  )
}

get_manta_tow <- function(name, domain_category = "reef") {
  paste0(
    url_join(BASE_URL, "data"),
    "?domain_name=", uenc(name),
    "&domain_category=", uenc(domain_category),
    "&data_type=manta"
  )
}

get_disturbances <- function(name, aggregation = "reef") {
  paste0(
    url_join(BASE_URL, "disturbance"),
    "?reef=", uenc(name),
    "&aggregation=", uenc(aggregation),
    "&zone=_"
  )
}

get_cots <- function(name, domain_category = "reef") {
  paste0(
    url_join(BASE_URL, "cots-by-domain"),
    "?domain_category=", uenc(domain_category),
    "&domain_name=", uenc(name)
  )
}

# ---- HTTP + JSON helpers ----
fetch_json_data <- function(url) {
  resp <- tryCatch(httr::GET(url), error = function(e) e)
  if (inherits(resp, "error")) {
    message("Request error for: ", url, "\n  ", resp$message)
    return(NULL)
  }
  
  status <- httr::status_code(resp)
  if (status != 200) {
    message("Failed (HTTP ", status, ") for: ", url)
    return(NULL)
  }
  
  txt <- httr::content(resp, as = "text", encoding = "UTF-8")
  
  # API sometimes returns [] for no data
  if (identical(trimws(txt), "[]") || nchar(trimws(txt)) == 0) return(NULL)
  
  # flatten=TRUE produces a data.frame for list-of-objects JSON
  out <- jsonlite::fromJSON(txt, flatten = TRUE)
  
  out
}

is_empty_payload <- function(x) {
  is.null(x) ||
    (is.data.frame(x) && nrow(x) == 0) ||
    (is.list(x) && length(x) == 0)
}

# Generic fetch-and-bind for many domains (reefs or sectors)
fetch_and_bind <- function(ids, url_fun, id_col, label = "data") {
  ids <- unique(ids)
  ids <- ids[!is.na(ids) & nzchar(ids)]
  
  if (length(ids) == 0) {
    message("No IDs provided for ", label, ".")
    return(tibble())
  }
  
  message("Fetching ", label, " for ", length(ids), " IDs ...")
  
  pb <- utils::txtProgressBar(min = 0, max = length(ids), style = 3)
  on.exit(try(close(pb), silent = TRUE), add = TRUE)
  
  out_list <- vector("list", length(ids))
  names(out_list) <- ids
  
  for (i in seq_along(ids)) {
    id <- ids[[i]]
    utils::setTxtProgressBar(pb, i)
    
    url <- url_fun(id)
    dat <- fetch_json_data(url)
    
    if (!is_empty_payload(dat)) {
      df <- tibble::as_tibble(dat)
      df[[id_col]] <- id
      out_list[[i]] <- df
    } else {
      out_list[[i]] <- NULL
    }
  }
  
  bind_rows(out_list)
}

# ---- Step 1: list of reefs ----
reef_info <- fetch_json_data(get_reef_info_url) %>% tibble::as_tibble()

if (nrow(reef_info) == 0) {
  stop("No reef info returned. Check API availability and BASE_URL.")
}

reef_names <- unique(reef_info$aims_reef_name)

# ---- Reef-level datasets ----
reef_photo_df <- fetch_and_bind(
  reef_names,
  url_fun = function(r) get_photo_transect(r, domain_category = "reef"),
  id_col  = "reef_name",
  label   = "photo-transect (reef)"
)

reef_manta_df <- fetch_and_bind(
  reef_names,
  url_fun = function(r) get_manta_tow(r, domain_category = "reef"),
  id_col  = "reef_name",
  label   = "manta tow (reef)"
)

reef_disturbance_df <- fetch_and_bind(
  reef_names,
  url_fun = function(r) get_disturbances(r, aggregation = "reef"),
  id_col  = "reef_name",
  label   = "disturbances (reef)"
)
# NOTE: In Julia you had to union columns manually; in R, bind_rows() fills missing
# columns automatically with NA, so no extra work is needed.

reef_cots_df <- fetch_and_bind(
  reef_names,
  url_fun = function(r) get_cots(r, domain_category = "reef"),
  id_col  = "reef_name",
  label   = "COTS (reef)"
)

# ---- Sector-level datasets ----
sectors <- unique(reef_info$a_sector)

sector_photo_df <- fetch_and_bind(
  sectors,
  url_fun = function(s) get_photo_transect(s, domain_category = "sector"),
  id_col  = "sector_name",
  label   = "photo-transect (sector)"
)

sector_manta_df <- fetch_and_bind(
  sectors,
  url_fun = function(s) get_manta_tow(s, domain_category = "sector"),
  id_col  = "sector_name",
  label   = "manta tow (sector)"
)

sector_cots_df <- fetch_and_bind(
  sectors,
  url_fun = function(s) get_cots(s, domain_category = "sector"),
  id_col  = "sector_name",
  label   = "COTS (sector)"
)

# ---- Optional: write outputs ----
out_dir <- "data/aims_ltmp"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
#
readr::write_csv(reef_info,             file.path(out_dir, "reef_info.csv"))
readr::write_csv(reef_photo_df,         file.path(out_dir, "reef_photo_transect.csv"))
readr::write_csv(reef_manta_df,         file.path(out_dir, "reef_manta.csv"))
readr::write_csv(reef_disturbance_df,   file.path(out_dir, "reef_disturbance.csv"))
readr::write_csv(reef_cots_df,          file.path(out_dir, "reef_cots.csv"))
readr::write_csv(sector_photo_df,       file.path(out_dir, "sector_photo_transect.csv"))
readr::write_csv(sector_manta_df,       file.path(out_dir, "sector_manta.csv"))
readr::write_csv(sector_cots_df,        file.path(out_dir, "sector_cots.csv"))
#
save(
  reef_info,
  reef_photo_df,
  reef_manta_df,
  reef_disturbance_df,
  reef_cots_df,
  sector_photo_df,
  sector_manta_df,
  sector_cots_df,
  file = file.path(out_dir, "aims_ltmp.RData")
)

message("Done. Objects available: reef_info, reef_photo_df, reef_manta_df, reef_disturbance_df, reef_cots_df, sector_photo_df, sector_manta_df, sector_cots_df")
