library(dplyr)
zarr_full_raw <- read.csv("data/processed/cheung_recreated_gbr_full.csv")
cat("Zarr rows:", nrow(zarr_full_raw), "\n")
cat("Zarr LABEL_ID head:\n")
print(head(zarr_full_raw$LABEL_ID))

cheung <- read.csv("data/raw/cheung_etal_2025/13068_2025_2634_MOESM4_ESM.csv")
cat("Cheung rows:", nrow(cheung), "\n")
cat("Cheung head LABEL_ID and year:\n")
print(head(cheung[, c("LABEL_ID", "year", "maxDHW")]))

# Check match between cheung and zarr
zarr_sub <- zarr_full_raw %>%
  filter(year %in% c(2016, 2017, 2020)) %>%
  mutate(LABEL_ID = toupper(trimws(as.character(LABEL_ID))))

cheung_match <- cheung %>%
  mutate(LABEL_ID = toupper(trimws(as.character(LABEL_ID)))) %>%
  left_join(zarr_sub, by = c("LABEL_ID", "year"))

cat("Matched rows non-NA MaxDHW:", sum(!is.na(cheung_match$ann_maxdhw)), "out of", nrow(cheung_match), "\n")
