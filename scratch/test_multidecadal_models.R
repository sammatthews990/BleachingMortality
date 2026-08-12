# test_multidecadal_models.R
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gbm)

cat("Testing full model suite on recreated Zarr N=262 and Multi-Decadal N=465...\n")

# Load Zarr
zarr_full_raw <- read.csv("data/processed/cheung_recreated_gbr_full.csv")

# Let's test data prep and BRT + SINDy + GLM
cat("Zarr dimensions:", dim(zarr_full_raw), "\n")
