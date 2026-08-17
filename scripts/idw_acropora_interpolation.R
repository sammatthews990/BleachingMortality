# scripts/idw_acropora_interpolation.R
# Fast Vectorized Spatial-Temporal Inverse Distance Weighting (IDW) for % Acropora

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(sf)
  library(stringr)
})

cat("=== Step 1: Loading AIMS LTMP Photo-Transect & Reef Coordinates ===\n")

load("data/aims_ltmp/aims_ltmp.RData")

dat_ref <- read.csv("data/AIMS-Reef_Reference.csv") %>%
  rename(Reef_Name = AIMS_REEF_NAME) %>%
  distinct(ReefName, .keep_all = TRUE)

dat_reef_info <- reef_info %>%
  select(aims_reef_name, LATITUDE = latitude, LONGITUDE = longitude, SECTOR = a_sector) %>%
  mutate(AIMS_REEF_NAME_cap = toupper(aims_reef_name)) %>%
  left_join(dat_ref %>% select(AIMS_REEF_NAME_cap, ReefName), by = "AIMS_REEF_NAME_cap", relationship = "many-to-many") %>%
  filter(!is.na(ReefName))

dat_cheung_coords <- read.csv("data/processed/cheung_recreated_gbr_full.csv") %>%
  distinct(LOC_NAME_S, .keep_all = TRUE) %>%
  select(ReefName = LOC_NAME_S, LATITUDE = lat, LONGITUDE = lon)

coords_df <- dat_reef_info %>%
  select(ReefName, LATITUDE, LONGITUDE, SECTOR) %>%
  bind_rows(dat_cheung_coords) %>%
  filter(!is.na(LATITUDE), !is.na(LONGITUDE)) %>%
  distinct(ReefName, .keep_all = TRUE)

cat("Coordinate lookup built for", nrow(coords_df), "reefs.\n")

cat("\n=== Step 2: Extracting Observed Prop. Acropora from LTMP Benthic Transects ===\n")

acrop_obs <- reef_photo_df %>%
  filter(!is.na(mean)) %>%
  group_by(domain_name, report_year, depth, reefpage_category) %>%
  summarise(mean_cover = mean(mean, na.rm = TRUE), .groups = "drop") %>%
  group_by(domain_name, report_year) %>%
  summarise(
    total_hc = sum(mean_cover[reefpage_category %in% c("Hard Coral", "HARD CORAL", "Hard coral") | grepl("Acropora|Montipora|Porites|Fav|Pocillopora", reefpage_category, ignore.case = TRUE)], na.rm = TRUE),
    acrop_cover = sum(mean_cover[grepl("Acropora", reefpage_category, ignore.case = TRUE)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(total_hc > 0) %>%
  mutate(prop_acropora_obs = acrop_cover / total_hc) %>%
  mutate(Reef_Name_cap = toupper(domain_name)) %>%
  left_join(dat_ref %>% select(AIMS_REEF_NAME_cap, ReefName) %>% distinct(AIMS_REEF_NAME_cap, .keep_all = TRUE), by = c("Reef_Name_cap" = "AIMS_REEF_NAME_cap")) %>%
  filter(!is.na(ReefName)) %>%
  select(ReefName, report_year, prop_acropora_obs) %>%
  group_by(ReefName, report_year) %>%
  summarise(prop_acropora_obs = mean(prop_acropora_obs, na.rm = TRUE), .groups = "drop") %>%
  left_join(coords_df, by = "ReefName") %>%
  filter(!is.na(LATITUDE), !is.na(LONGITUDE))

cat("Extracted", nrow(acrop_obs), "observed reef-year prop_acropora records across", length(unique(acrop_obs$ReefName)), "benthic reefs.\n")

cat("\n=== Step 3: Fast Vectorized Matrix IDW Interpolation ===\n")

dat_cheung <- read.csv("data/processed/cheung_recreated_gbr_full.csv")

target_reefs <- dat_cheung %>%
  distinct(LOC_NAME_S, .keep_all = TRUE) %>%
  select(ReefName = LOC_NAME_S, LATITUDE = lat, LONGITUDE = lon) %>%
  left_join(coords_df %>% select(ReefName, SECTOR), by = "ReefName") %>%
  distinct(ReefName, .keep_all = TRUE)

years_all <- 1985:2025

target_grid <- expand.grid(
  ReefName = target_reefs$ReefName,
  report_year = years_all,
  stringsAsFactors = FALSE
) %>%
  left_join(target_reefs, by = "ReefName")

# Match observed directly first
target_grid <- target_grid %>%
  left_join(acrop_obs %>% select(ReefName, report_year, prop_acropora_obs), by = c("ReefName", "report_year"))

# Pre-calculate Distance Matrix
unique_target_sf <- target_reefs %>% st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
unique_obs_reefs <- acrop_obs %>% distinct(ReefName, LONGITUDE, LATITUDE) %>% st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

dist_matrix <- st_distance(unique_target_sf, unique_obs_reefs)
dist_km <- matrix(as.numeric(dist_matrix) / 1000, nrow = nrow(unique_target_sf), ncol = nrow(unique_obs_reefs))
rownames(dist_km) <- unique_target_sf$ReefName
colnames(dist_km) <- unique_obs_reefs$ReefName

# Local reef mean lookup
reef_means <- acrop_obs %>%
  group_by(ReefName) %>%
  summarise(reef_mean_acrop = mean(prop_acropora_obs, na.rm = TRUE), .groups = "drop")

# Sector annual mean lookup
sector_annual_means <- acrop_obs %>%
  group_by(SECTOR, report_year) %>%
  summarise(sec_year_mean = mean(prop_acropora_obs, na.rm = TRUE), .groups = "drop")

overall_mean <- mean(acrop_obs$prop_acropora_obs, na.rm = TRUE)

target_grid$idw_val <- NA_real_

# Loop by Year (40 iterations instead of 280,000!)
for (yr in years_all) {
  sub_obs <- acrop_obs %>% filter(abs(report_year - yr) <= 2)
  if (nrow(sub_obs) == 0) next
  
  obs_reefs_yr <- sub_obs$ReefName
  obs_vals_yr <- sub_obs$prop_acropora_obs
  
  # Distance matrix subset for this year's observations
  D_yr <- dist_km[, obs_reefs_yr, drop = FALSE]
  W <- 1 / (D_yr + 0.1)^2
  W[D_yr > 150] <- 0 # zero weight outside 150 km
  
  w_sum <- rowSums(W)
  valid_rows <- which(w_sum > 0)
  
  if (length(valid_rows) > 0) {
    idw_calc <- (W[valid_rows, , drop = FALSE] %*% obs_vals_yr) / w_sum[valid_rows]
    reef_names_valid <- rownames(dist_km)[valid_rows]
    
    grid_idx <- which(target_grid$report_year == yr & target_grid$ReefName %in% reef_names_valid)
    # Match order
    m_idx <- match(target_grid$ReefName[grid_idx], reef_names_valid)
    target_grid$idw_val[grid_idx] <- idw_calc[m_idx]
  }
}

# Cascade fallback
target_grid <- target_grid %>%
  left_join(reef_means, by = "ReefName") %>%
  left_join(sector_annual_means, by = c("SECTOR", "report_year")) %>%
  mutate(
    prop_acropora_idw = case_when(
      !is.na(prop_acropora_obs) ~ prop_acropora_obs,
      !is.na(idw_val)            ~ idw_val,
      !is.na(reef_mean_acrop)   ~ reef_mean_acrop,
      !is.na(sec_year_mean)     ~ sec_year_mean,
      TRUE                       ~ overall_mean
    ),
    imputation_method = case_when(
      !is.na(prop_acropora_obs) ~ "Observed",
      !is.na(idw_val)            ~ "Spatial-Temporal IDW",
      !is.na(reef_mean_acrop)   ~ "Local Reef Mean",
      !is.na(sec_year_mean)     ~ "Sector Annual Mean",
      TRUE                       ~ "Overall Dataset Mean"
    )
  )

cat("\n=== Step 4: Summary of Vectorized Interpolated Prop. Acropora ===\n")
print(table(target_grid$imputation_method))
cat("\nSummary stats of IDW prop_acropora:\n")
print(summary(target_grid$prop_acropora_idw))

# Save interpolated object
saveRDS(target_grid, "data/processed/acropora_interpolated_df.rds")
cat("\nSaved interpolated dataset to data/processed/acropora_interpolated_df.rds\n")
