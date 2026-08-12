# Script to test and compare:
# 1. Cheung Benchmark vs Recreated Zarr on Combined (N=262, 2016-2020)
# 2. Modern Expanded Dataset (1998-2024, including 2022 and 2024)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(gbm)
})

# Load datasets
zarr_df <- read.csv("data/processed/cheung_recreated_gbr_full.csv")
load("data/Cheungetal2025/01_sstvar_blchrf.RData")

pred_vars <- c("mcur_90", "cloudp_90", "secc3m", "histmDHW6", "yrsince6", "histmDHW4", "yrsince4", "winyear_sd", "winyear_mean")

# 1. Benchmark Predictors (2016, 2017, 2020)
bm_preds <- sstvar_blch2 %>%
  select(LABEL_id, year, lat = Y, lon = X, Sector, all_of(pred_vars)) %>%
  group_by(LABEL_id, year, Sector) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(year = as.numeric(year), LABEL_id = toupper(trimws(as.character(LABEL_id))))

# 2. Recreated Zarr Predictors (2016, 2017, 2020, 2022, 2024)
zarr_preds <- zarr_df %>%
  select(LABEL_id = LABEL_ID, year, lat, lon, all_of(pred_vars[pred_vars %in% names(zarr_df)])) %>%
  group_by(LABEL_id, year) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(year = as.numeric(year), LABEL_id = toupper(trimws(as.character(LABEL_id))))

cat("Benchmark unique reef-years:", nrow(bm_preds), "\n")
cat("Zarr unique reef-years:", nrow(zarr_preds), "\n")

# Save summary object
saveRDS(list(bm_preds = bm_preds, zarr_preds = zarr_preds), "data/processed/cheung_comparison_benchmark.rds")
cat("Saved comparison object.\n")
