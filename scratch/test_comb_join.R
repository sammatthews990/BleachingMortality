library(dplyr)

# Load data and check dat.mod.mant and cheung_predictors
load("data/GBR_AIMS/MantaTow_CoralCover_AIMS.RData")
# Or let's test with the objects created in QMD
zarr_full_raw <- read.csv("data/processed/cheung_recreated_gbr_full.csv")
cat("Zarr rows:", nrow(zarr_full_raw), "\n")
print(head(zarr_full_raw[, c("LABEL_ID", "year", "ann_maxdhw", "histmDHW6")]))
