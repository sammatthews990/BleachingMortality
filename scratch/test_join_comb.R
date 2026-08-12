library(dplyr)

# Test creating dat.ml.comb.clean and matching with zarr
zarr_full_raw <- read.csv("data/processed/cheung_recreated_gbr_full.csv")

# Let's inspect how dat.ml.mant.clean has ReefID and DHWYear
# We can load the saved rdata if available or check names
cat("Zarr years:", unique(zarr_full_raw$year), "\n")
cat("Zarr LABEL_IDs count:", length(unique(zarr_full_raw$LABEL_ID)), "\n")
