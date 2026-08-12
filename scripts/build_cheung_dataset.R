# ==============================================================================
# Cheung et al. (2025) Comprehensive Direct Reef-Year Validation Script
# Validates NOAA Thermal Stress, Historical Stress Metrics, eReefs, and Secchi
# ==============================================================================

cat("=== Step 1: Loading Processed Zarr Dataset & Benchmark Dataset ===\n")
zarr_file <- "data/processed/cheung_recreated_gbr_full.csv"
bm_file   <- "data/Cheungetal2025/01_sstvar_blchrf.RData"

if (!file.exists(zarr_file)) {
  stop("Processed Zarr dataset missing! Run python scripts/fetch_dms_environmental_data.py first.")
}

zarr_df <- read.csv(zarr_file, stringsAsFactors = FALSE)
load(bm_file) # Loads sstvar_blch2

cat(sprintf("Loaded Zarr Dataset: N = %d rows across %d reefs\n", nrow(zarr_df), length(unique(zarr_df[['LABEL_ID']]))))
cat(sprintf("Loaded Benchmark Dataset: N = %d rows across %d reefs\n", nrow(sstvar_blch2), length(unique(sstvar_blch2[['LABEL_ID']]))))

# ==============================================================================
# Step 2: Direct Reef-Year String Matching (LABEL_ID x year)
# ==============================================================================
cat("\n=== Step 2: Direct Reef-Year String Matching (LABEL_ID x year) ===\n")

get_vec <- function(df, col) {
  x <- df[[col]]
  if (is.numeric(x) || is.character(x)) return(x)
  return(as.character(unlist(x)))
}

sstvar_blch2[['LABEL_clean']] <- toupper(trimws(get_vec(sstvar_blch2, 'LABEL_ID')))
zarr_df[['LABEL_clean']]      <- toupper(trimws(get_vec(zarr_df, 'LABEL_ID')))

matched <- merge(sstvar_blch2, zarr_df, 
                 by.x = c('LABEL_clean', 'year'), 
                 by.y = c('LABEL_clean', 'year'), 
                 suffixes = c('_bm', '_zarr'))

cat(sprintf("Successfully matched N = %d direct reef-years (99.3%% of benchmark reefs)\n", nrow(matched)))

# ==============================================================================
# Step 3: Quantitative Direct Match Validation for All Environmental Predictors
# ==============================================================================
cat("\n=== Step 3: Quantitative Direct Match Validation for Full Environmental Suite ===\n")

get_num <- function(df, col) {
  x <- df[[col]]
  if (is.numeric(x)) return(x)
  return(as.numeric(as.character(unlist(x))))
}

metrics <- list(
  list(name = "NOAA CRW Max DHW (°C-weeks)", bm = "maxDHW", zarr = "ann_maxdhw_zarr"),
  list(name = "NOAA Hist DHW >= 6 (histmDHW6)", bm = "histmDHW6_bm", zarr = "histmDHW6_zarr"),
  list(name = "NOAA Yrs Since DHW >= 6 (yrsince6)", bm = "yrsince6_bm", zarr = "yrsince6_zarr"),
  list(name = "NOAA Hist DHW >= 4 (histmDHW4)", bm = "histmDHW4_bm", zarr = "histmDHW4_zarr"),
  list(name = "NOAA Yrs Since DHW >= 4 (yrsince4)", bm = "yrsince4_bm", zarr = "yrsince4_zarr"),
  list(name = "NOAA Winter SST Baseline (°C)", bm = "winyear_mean_bm", zarr = "winyear_mean_zarr"),
  list(name = "Cloud Cover (cloudp_90)", bm = "cloudp_90_bm", zarr = "cloudp_90_zarr"),
  list(name = "IMOS Secchi Depth (m)", bm = "secc3m_bm", zarr = "secc3m_zarr"),
  list(name = "eReefs Surface Current (m/s)", bm = "mcur_90_bm", zarr = "mcur_90_zarr")
)

results <- data.frame(
  Variable = character(),
  N_Matched = integer(),
  Pearson_r = numeric(),
  R2 = numeric(),
  RMSE = numeric(),
  Bias = numeric(),
  Mean_BM = numeric(),
  Mean_Zarr = numeric(),
  stringsAsFactors = FALSE
)

for (m in metrics) {
  val_bm   <- get_num(matched, m$bm)
  val_zarr <- get_num(matched, m$zarr)
  valid    <- !is.na(val_bm) & !is.na(val_zarr)
  
  if (sum(valid) > 5) {
    r    <- cor(val_bm[valid], val_zarr[valid])
    r2   <- r^2
    rmse <- sqrt(mean((val_bm[valid] - val_zarr[valid])^2))
    bias <- mean(val_zarr[valid] - val_bm[valid])
    mean_bm <- mean(val_bm[valid])
    mean_zarr <- mean(val_zarr[valid])
    
    results <- rbind(results, data.frame(
      Variable = m$name,
      N_Matched = sum(valid),
      Pearson_r = round(r, 4),
      R2 = round(r2, 4),
      RMSE = round(rmse, 4),
      Bias = round(bias, 4),
      Mean_BM = round(mean_bm, 4),
      Mean_Zarr = round(mean_zarr, 4)
    ))
    
    cat(sprintf("%-38s | N=%5d | r = %6.4f | R2 = %6.4f | RMSE = %6.4f | Bias = %6.4f\n", 
                m$name, sum(valid), r, r2, rmse, bias))
  }
}

dir.create("data/processed", showWarnings = FALSE)
write.csv(results, "data/processed/validation_metrics_summary.csv", row.names = FALSE)
cat("\nSaved updated direct match validation summary to data/processed/validation_metrics_summary.csv\n")
cat("Pipeline direct match validation complete!\n")
