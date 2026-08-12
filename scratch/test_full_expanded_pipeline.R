# scratch/test_full_expanded_pipeline.R
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gbm)

cat("1. Testing full model suite logic in test script...\n")

# Load Zarr
zarr_full_raw <- read.csv("data/processed/cheung_recreated_gbr_full.csv")

# Create mock/simulated test data matching dat.ml.comb.clean structure
set.seed(42)
N_test <- 262
dat_sim <- data.frame(
    ReefName = paste0("Reef_", sample(1:50, N_test, replace=TRUE)),
    ReefID = paste0("09-", sample(100:400, N_test, replace=TRUE)),
    DHWYear = sample(c(2016, 2017, 2020), N_test, replace=TRUE),
    Sector = sample(c("Northern", "Central", "Southern"), N_test, replace=TRUE),
    MaxDHW.mean = pmax(0, rnorm(N_test, mean=7, sd=3.5)),
    secc3m = pmax(2, rnorm(N_test, mean=15, sd=5)),
    winyear_mean = rnorm(N_test, mean=24, sd=1.5),
    histmDHW6 = pmax(0, rnorm(N_test, mean=2, sd=2)),
    mcur_90 = pmax(0.05, rnorm(N_test, mean=0.3, sd=0.1)),
    winyear_sd = pmax(0.1, rnorm(N_test, mean=0.8, sd=0.2)),
    cloudp_90 = pmax(0.1, pmin(0.95, rnorm(N_test, mean=0.6, sd=0.15))),
    yrsince6 = sample(1:15, N_test, replace=TRUE),
    prop_acropora = pmax(0.01, pmin(0.99, rnorm(N_test, mean=0.4, sd=0.2))),
    Mort.prop = pmax(0, pmin(1, plogis(-2 + 0.35 * rnorm(N_test, mean=7, sd=3.5) - 0.05 * rnorm(N_test, mean=15, sd=5)))),
    Dataset = factor(sample(c("Manta Tow", "LTMP Benthic", "MMP Inshore"), N_test, replace=TRUE), levels = c("Manta Tow", "LTMP Benthic", "MMP Inshore"))
)

dat_sim$Mort.prop.nudge <- (dat_sim$Mort.prop * (N_test - 1) + 0.5) / N_test

# SINDy solver
stlsq <- function(Theta_mat, Y_vec, lambda = 0.08, max_iter = 10) {
    Xi <- solve(t(Theta_mat) %*% Theta_mat + 1e-4 * diag(ncol(Theta_mat))) %*% t(Theta_mat) %*% Y_vec
    for (k in 1:max_iter) {
        small_inds <- abs(Xi) < lambda
        Xi[small_inds] <- 0
        big_inds <- !small_inds
        if (sum(big_inds) == 0) break
        Xi[big_inds] <- solve(t(Theta_mat[, big_inds]) %*% Theta_mat[, big_inds] + 1e-4 * diag(sum(big_inds))) %*% t(Theta_mat[, big_inds]) %*% Y_vec
    }
    return(Xi)
}

# Test BRT fitting
f_brt <- Mort.prop ~ MaxDHW.mean + secc3m + winyear_mean + histmDHW6 + mcur_90 + winyear_sd + cloudp_90 + yrsince6 + Dataset
brt_fit <- gbm(
    f_brt, data = dat_sim, distribution = "gaussian",
    n.trees = 1000, interaction.depth = 2, shrinkage = 0.005,
    n.minobsinnode = 5, bag.fraction = 0.70, cv.folds = 5
)
best_t <- gbm.perf(brt_fit, method = "cv", plot.it = FALSE)
cat("BRT best tree:", best_t, "\n")

# Test 1D PDP bootstrap
grid_sz <- 30
B_test <- 20
p_vars <- c("MaxDHW.mean", "secc3m", "winyear_mean", "histmDHW6", "mcur_90", "winyear_sd", "cloudp_90", "yrsince6")
boot_list <- list()
for (pv in p_vars) {
    x_seq <- seq(min(dat_sim[[pv]]), max(dat_sim[[pv]]), length.out = grid_sz)
    mat_y <- matrix(NA, nrow = B_test, ncol = grid_sz)
    grid_base <- data.frame(lapply(dat_sim[p_vars], median))
    grid_base$Dataset <- factor("Manta Tow", levels = levels(factor(dat_sim$Dataset)))
    grid_full <- grid_base[rep(1, grid_sz), ]
    grid_full[[pv]] <- x_seq
    for (b in 1:B_test) {
        idx <- sample(1:nrow(dat_sim), replace=TRUE)
        fb <- gbm(f_brt, data=dat_sim[idx, ], distribution="gaussian", n.trees=best_t, interaction.depth=2, shrinkage=0.005, n.minobsinnode=5, bag.fraction=0.7, verbose=FALSE)
        mat_y[b, ] <- predict(fb, newdata=grid_full, n.trees=best_t, type="response")
    }
    boot_list[[pv]] <- data.frame(
        predictor = pv, x = x_seq,
        median_yhat = apply(mat_y, 2, median),
        lo_80 = apply(mat_y, 2, quantile, 0.10),
        hi_80 = apply(mat_y, 2, quantile, 0.90),
        lo_95 = apply(mat_y, 2, quantile, 0.025),
        hi_95 = apply(mat_y, 2, quantile, 0.975)
    )
}
pdp_df <- do.call(rbind, boot_list)
cat("PDP bootstrap generated successfully! Rows:", nrow(pdp_df), "\n")
