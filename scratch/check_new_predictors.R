library(dplyr)
source("scratch/check_disturbances_cheung.R") # produces dat.ml.mant.raw

cat("Available columns in dat.ml.mant.raw:\n")
print(colnames(dat.ml.mant.raw))

# Check for predictor names
years_cols <- grep("year|since|hist|win", colnames(dat.ml.mant.raw), value = TRUE, ignore.case = TRUE)
cat("\nMatching predictor candidate columns:\n")
print(years_cols)
