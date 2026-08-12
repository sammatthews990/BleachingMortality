# scratch/test_full_qmd_replacement.R
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(gbm)

cat("Testing full chunks for Section 1 (Zarr N=262) and Section 2 (Multi-Decadal N=465)...\n")

# Load Zarr
zarr_full_raw <- read.csv("data/processed/cheung_recreated_gbr_full.csv")

# Create mock data
set.seed(42)
N_comb <- 262
dat.ml.comb.clean <- data.frame(
    ReefName = paste0("Reef_", sample(1:50, N_comb, replace=TRUE)),
    ReefID = paste0("09-", sample(100:400, N_comb, replace=TRUE)),
    DHWYear = sample(c(2016, 2017, 2020), N_comb, replace=TRUE),
    Sector = sample(c("Northern", "Central", "Southern"), N_comb, replace=TRUE),
    MaxDHW.mean = pmax(0, rnorm(N_comb, mean=7, sd=3.5)),
    secc3m = pmax(2, rnorm(N_comb, mean=15, sd=5)),
    winyear_mean = rnorm(N_comb, mean=24, sd=1.5),
    histmDHW6 = pmax(0, rnorm(N_comb, mean=2, sd=2)),
    mcur_90 = pmax(0.05, rnorm(N_comb, mean=0.3, sd=0.1)),
    winyear_sd = pmax(0.1, rnorm(N_comb, mean=0.8, sd=0.2)),
    cloudp_90 = pmax(0.1, pmin(0.95, rnorm(N_comb, mean=0.6, sd=0.15))),
    yrsince6 = sample(1:15, N_comb, replace=TRUE),
    prop_acropora = pmax(0.01, pmin(0.99, rnorm(N_comb, mean=0.4, sd=0.2))),
    Mort.prop = pmax(0, pmin(1, plogis(-2 + 0.35 * rnorm(N_comb, mean=7, sd=3.5) - 0.05 * rnorm(N_comb, mean=15, sd=5)))),
    Dataset = factor(sample(c("Manta Tow", "LTMP Benthic", "MMP Inshore"), N_comb, replace=TRUE), levels = c("Manta Tow", "LTMP Benthic", "MMP Inshore"))
)
dat.ml.comb.clean$Mort.prop.nudge <- (dat.ml.comb.clean$Mort.prop * (N_comb - 1) + 0.5) / N_comb

# Benchmark metrics mock
met1_c <- c(Dev=35.1, Train_R2=0.351, CV_R2=0.320, Train_RMSE=0.142, CV_RMSE=0.145)
met2_c <- c(Dev=48.2, Train_R2=0.482, CV_R2=0.421, Train_RMSE=0.126, CV_RMSE=0.131)
met4_c <- c(Dev=44.0, Train_R2=0.440, CV_R2=0.395, Train_RMSE=0.130, CV_RMSE=0.135)
met5_rand_c <- c(Dev=54.5, Train_R2=0.545, CV_R2=0.472, Train_RMSE=0.118, CV_RMSE=0.125)
met5_sp_c <- c(Dev=54.5, Train_R2=0.545, CV_R2=0.412, Train_RMSE=0.118, CV_RMSE=0.134)

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

calc_comb_metrics <- function(y_true, y_pred, y_pred_cv = NULL, SS_tot = NULL) {
    if (is.null(SS_tot)) SS_tot <- sum((y_true - mean(y_true))^2)
    y_pred_b <- pmin(pmax(y_pred, 0), 1)
    ss_res <- sum((y_true - y_pred_b)^2)
    dev_exp <- (1 - ss_res / SS_tot) * 100
    tr_r2 <- max(0, cor(y_true, y_pred_b)^2)
    tr_rmse <- sqrt(mean((y_true - y_pred_b)^2))

    if (!is.null(y_pred_cv)) {
        y_cv_b <- pmin(pmax(y_pred_cv, 0), 1)
        cv_r2 <- max(0, cor(y_true, y_cv_b)^2)
        cv_rmse <- sqrt(mean((y_true - y_cv_b)^2))
    } else {
        cv_r2 <- NA
        cv_rmse <- NA
    }
    return(c(Dev = dev_exp, Train_R2 = tr_r2, CV_R2 = cv_r2, Train_RMSE = tr_rmse, CV_RMSE = cv_rmse))
}

# --- Test 1. Zarr N=262 Preparation & Fitting ---
cat("Running Stage 1 (Zarr N=262)...\n")
zarr_preds_sub <- zarr_full_raw %>%
    filter(year %in% c(2016, 2017, 2020)) %>%
    select(LABEL_ID, year, MaxDHW_zarr = ann_maxdhw, secc3m_zarr = secc3m,
           winyear_mean_zarr = winyear_mean, histmDHW6_zarr = histmDHW6,
           mcur_90_zarr = mcur_90, winyear_sd_zarr = winyear_sd,
           cloudp_90_zarr = cloudp_90, yrsince6_zarr = yrsince6) %>%
    mutate(LABEL_ID = toupper(trimws(as.character(LABEL_ID))))

dat.ml.val <- dat.ml.comb.clean %>%
    mutate(
        ReefID_clean  = toupper(trimws(as.character(ReefID))),
        DHWYear_num   = as.numeric(DHWYear)
    ) %>%
    left_join(
        zarr_preds_sub,
        by = c("ReefID_clean" = "LABEL_ID", "DHWYear_num" = "year")
    )

pred_z_vars <- c("MaxDHW_zarr", "secc3m_zarr", "winyear_mean_zarr", "histmDHW6_zarr",
                 "mcur_90_zarr", "winyear_sd_zarr", "cloudp_90_zarr", "yrsince6_zarr")
for (v in pred_z_vars) {
    if (v %in% names(dat.ml.val)) {
        med_v <- median(dat.ml.val[[v]], na.rm = TRUE)
        if (is.na(med_v)) med_v <- median(dat.ml.comb.clean[[gsub("_zarr", "", v)]], na.rm = TRUE)
        dat.ml.val[[v]][is.na(dat.ml.val[[v]])] <- med_v
    }
}
dat.ml.val$Dataset <- factor(dat.ml.val$Dataset, levels = c("Manta Tow", "LTMP Benthic", "MMP Inshore"))

y_val <- dat.ml.val$Mort.prop
SS_tot_val <- sum((y_val - mean(y_val))^2)
N_val <- length(y_val)
set.seed(42)
cv_folds_val <- sample(rep(1:5, length.out = N_val))

# 1. GLM Univariate
m1_zarr <- glm(Mort.prop ~ MaxDHW_zarr, data = dat.ml.val, family = quasibinomial)
p1_z_in <- predict(m1_zarr, type = "response")
p1_z_cv <- numeric(N_val)
for (f in 1:5) {
    tr <- dat.ml.val[cv_folds_val != f, ]
    va <- dat.ml.val[cv_folds_val == f, ]
    fit_f <- glm(Mort.prop ~ MaxDHW_zarr, data = tr, family = quasibinomial)
    p1_z_cv[cv_folds_val == f] <- predict(fit_f, newdata = va, type = "response")
}
met1_z <- calc_comb_metrics(y_val, p1_z_in, p1_z_cv, SS_tot_val)

# 2. GLM Full 8
f_full_z <- Mort.prop ~ MaxDHW_zarr + secc3m_zarr + winyear_mean_zarr + histmDHW6_zarr + mcur_90_zarr + winyear_sd_zarr + cloudp_90_zarr + yrsince6_zarr + prop_acropora + Dataset
m2_zarr <- glm(f_full_z, data = dat.ml.val, family = quasibinomial)
p2_z_in <- predict(m2_zarr, type = "response")
p2_z_cv <- numeric(N_val)
for (f in 1:5) {
    tr <- dat.ml.val[cv_folds_val != f, ]
    va <- dat.ml.val[cv_folds_val == f, ]
    fit_f <- glm(f_full_z, data = tr, family = quasibinomial)
    p2_z_cv[cv_folds_val == f] <- predict(fit_f, newdata = va, type = "response")
}
met2_z <- calc_comb_metrics(y_val, p2_z_in, p2_z_cv, SS_tot_val)

# 3. SINDy on Zarr
X_z <- dat.ml.val %>% select(MaxDHW_zarr, secc3m_zarr, winyear_mean_zarr, histmDHW6_zarr, mcur_90_zarr, winyear_sd_zarr, cloudp_90_zarr, yrsince6_zarr)
Theta_z <- cbind(
    1, X_z$MaxDHW_zarr, X_z$secc3m_zarr, X_z$winyear_mean_zarr, X_z$histmDHW6_zarr,
    X_z$mcur_90_zarr, X_z$winyear_sd_zarr, X_z$cloudp_90_zarr, X_z$yrsince6_zarr,
    X_z$MaxDHW_zarr^2, X_z$secc3m_zarr^2, X_z$winyear_mean_zarr^2, X_z$histmDHW6_zarr^2,
    X_z$MaxDHW_zarr * X_z$secc3m_zarr, X_z$MaxDHW_zarr * X_z$winyear_mean_zarr,
    X_z$MaxDHW_zarr * X_z$histmDHW6_zarr, X_z$MaxDHW_zarr * X_z$cloudp_90_zarr
)
Xi_zarr <- stlsq(Theta_z, y_val, lambda = 0.08)
p4_z_in <- Theta_z %*% Xi_zarr
p4_z_cv <- numeric(N_val)
for (f in 1:5) {
    tr_idx <- cv_folds_val != f
    va_idx <- cv_folds_val == f
    Xi_f <- stlsq(Theta_z[tr_idx, ], y_val[tr_idx], lambda = 0.08)
    p4_z_cv[va_idx] <- Theta_z[va_idx, ] %*% Xi_f
}
met4_z <- calc_comb_metrics(y_val, p4_z_in, p4_z_cv, SS_tot_val)

# 4. BRT on Zarr
f_brt_z <- Mort.prop ~ MaxDHW_zarr + secc3m_zarr + winyear_mean_zarr + histmDHW6_zarr + mcur_90_zarr + winyear_sd_zarr + cloudp_90_zarr + yrsince6_zarr + prop_acropora + Dataset
brt_zarr_fit <- gbm(
    f_brt_z, data = dat.ml.val, distribution = "gaussian",
    n.trees = 1500, interaction.depth = 2, shrinkage = 0.005,
    n.minobsinnode = 5, bag.fraction = 0.70, cv.folds = 5
)
best_zarr_t <- gbm.perf(brt_zarr_fit, method = "cv", plot.it = FALSE)
p5_z_in <- predict(brt_zarr_fit, newdata = dat.ml.val, n.trees = best_zarr_t, type = "response")

p5_z_cv <- numeric(N_val)
for (f in 1:5) {
    tr_f <- dat.ml.val[cv_folds_val != f, ]
    va_f <- dat.ml.val[cv_folds_val == f, ]
    fit_rf <- gbm(
        f_brt_z, data = tr_f, distribution = "gaussian", n.trees = 1000,
        interaction.depth = 2, shrinkage = 0.005, n.minobsinnode = 5, bag.fraction = 0.70
    )
    p5_z_cv[cv_folds_val == f] <- predict(fit_rf, newdata = va_f, n.trees = best_zarr_t, type = "response")
}
met5_z_rand <- calc_comb_metrics(y_val, p5_z_in, p5_z_cv, SS_tot_val)

sec_val <- dat.ml.val$Sector
if (is.null(sec_val) || length(unique(sec_val)) < 3) sec_val <- dat.ml.val$ReefName
unique_sec_val <- unique(sec_val)
p5_z_sp_cv <- numeric(N_val)
for (s in unique_sec_val) {
    val_i <- which(sec_val == s)
    tr_i  <- which(sec_val != s)
    if (length(tr_i) > 10) {
        fit_s <- gbm(
            f_brt_z, data = dat.ml.val[tr_i, ], distribution = "gaussian",
            n.trees = 1000, interaction.depth = 2, shrinkage = 0.005,
            n.minobsinnode = 5, bag.fraction = 0.70, verbose = FALSE
        )
        p5_z_sp_cv[val_i] <- predict(fit_s, newdata = dat.ml.val[val_i, ], n.trees = best_zarr_t, type = "response")
    } else {
        p5_z_sp_cv[val_i] <- p5_z_in[val_i]
    }
}
met5_z_sp <- calc_comb_metrics(y_val, p5_z_in, p5_z_sp_cv, SS_tot_val)

cat("Stage 1 completed successfully!\n")

# --- Test 2. Multi-Decadal N=465 Preparation & Fitting ---
cat("Running Stage 2 (Multi-Decadal N=465)...\n")
N_exp <- 465
dat.ml.expanded <- data.frame(
    ReefName = paste0("Reef_", sample(1:80, N_exp, replace=TRUE)),
    ReefID = paste0("09-", sample(100:600, N_exp, replace=TRUE)),
    DHWYear = sample(c(1998, 2002, 2016, 2017, 2020, 2022, 2024), N_exp, replace=TRUE),
    Sector = sample(c("Northern", "Central", "Southern"), N_exp, replace=TRUE),
    MaxDHW_val = pmax(0, rnorm(N_exp, mean=7.5, sd=4)),
    secc3m = pmax(2, rnorm(N_exp, mean=15, sd=5)),
    winyear_mean = rnorm(N_exp, mean=24, sd=1.5),
    histmDHW6 = pmax(0, rnorm(N_exp, mean=3, sd=2.5)),
    mcur_90 = pmax(0.05, rnorm(N_exp, mean=0.3, sd=0.1)),
    winyear_sd = pmax(0.1, rnorm(N_exp, mean=0.8, sd=0.2)),
    cloudp_90 = pmax(0.1, pmin(0.95, rnorm(N_exp, mean=0.6, sd=0.15))),
    yrsince6 = sample(1:15, N_exp, replace=TRUE),
    prop_acropora = pmax(0.01, pmin(0.99, rnorm(N_exp, mean=0.4, sd=0.2))),
    Mort.prop = pmax(0, pmin(1, plogis(-2.2 + 0.38 * rnorm(N_exp, mean=7.5, sd=4) - 0.06 * rnorm(N_exp, mean=15, sd=5)))),
    Dataset = factor(sample(c("Manta Tow", "LTMP Benthic", "MMP Inshore"), N_exp, replace=TRUE), levels = c("Manta Tow", "LTMP Benthic", "MMP Inshore"))
)
dat.ml.expanded$Mort.prop.nudge <- (dat.ml.expanded$Mort.prop * (N_exp - 1) + 0.5) / N_exp

y_exp <- dat.ml.expanded$Mort.prop
SS_tot_exp <- sum((y_exp - mean(y_exp))^2)
set.seed(42)
cv_folds_exp <- sample(rep(1:5, length.out = N_exp))

# 1. GLM Univariate
m1_exp <- glm(Mort.prop ~ MaxDHW_val, data = dat.ml.expanded, family = quasibinomial)
p1_exp_in <- predict(m1_exp, type = "response")
p1_exp_cv <- numeric(N_exp)
for (f in 1:5) {
    tr <- dat.ml.expanded[cv_folds_exp != f, ]
    va <- dat.ml.expanded[cv_folds_exp == f, ]
    fit_f <- glm(Mort.prop ~ MaxDHW_val, data = tr, family = quasibinomial)
    p1_exp_cv[cv_folds_exp == f] <- predict(fit_f, newdata = va, type = "response")
}
met1_exp <- calc_comb_metrics(y_exp, p1_exp_in, p1_exp_cv, SS_tot_exp)

# 2. GLM Full 8
f_full_exp <- Mort.prop ~ MaxDHW_val + secc3m + winyear_mean + histmDHW6 + mcur_90 + winyear_sd + cloudp_90 + yrsince6 + prop_acropora + Dataset
m2_exp <- glm(f_full_exp, data = dat.ml.expanded, family = quasibinomial)
p2_exp_in <- predict(m2_exp, type = "response")
p2_exp_cv <- numeric(N_exp)
for (f in 1:5) {
    tr <- dat.ml.expanded[cv_folds_exp != f, ]
    va <- dat.ml.expanded[cv_folds_exp == f, ]
    fit_f <- glm(f_full_exp, data = tr, family = quasibinomial)
    p2_exp_cv[cv_folds_exp == f] <- predict(fit_f, newdata = va, type = "response")
}
met2_exp <- calc_comb_metrics(y_exp, p2_exp_in, p2_exp_cv, SS_tot_exp)

# 3. SINDy on Expanded
X_exp <- dat.ml.expanded %>% select(MaxDHW_val, secc3m, winyear_mean, histmDHW6, mcur_90, winyear_sd, cloudp_90, yrsince6)
Theta_exp <- cbind(
    1, X_exp$MaxDHW_val, X_exp$secc3m, X_exp$winyear_mean, X_exp$histmDHW6,
    X_exp$mcur_90, X_exp$winyear_sd, X_exp$cloudp_90, X_exp$yrsince6,
    X_exp$MaxDHW_val^2, X_exp$secc3m^2, X_exp$winyear_mean^2, X_exp$histmDHW6^2,
    X_exp$MaxDHW_val * X_exp$secc3m, X_exp$MaxDHW_val * X_exp$winyear_mean,
    X_exp$MaxDHW_val * X_exp$histmDHW6, X_exp$MaxDHW_val * X_exp$cloudp_90
)
Xi_exp <- stlsq(Theta_exp, y_exp, lambda = 0.08)
p4_exp_in <- Theta_exp %*% Xi_exp
p4_exp_cv <- numeric(N_exp)
for (f in 1:5) {
    tr_idx <- cv_folds_exp != f
    va_idx <- cv_folds_exp == f
    Xi_f <- stlsq(Theta_exp[tr_idx, ], y_exp[tr_idx], lambda = 0.08)
    p4_exp_cv[va_idx] <- Theta_exp[va_idx, ] %*% Xi_f
}
met4_exp <- calc_comb_metrics(y_exp, p4_exp_in, p4_exp_cv, SS_tot_exp)

# 4. BRT on Expanded
f_brt_exp <- Mort.prop ~ MaxDHW_val + secc3m + winyear_mean + histmDHW6 + mcur_90 + winyear_sd + cloudp_90 + yrsince6 + prop_acropora + Dataset
brt_exp_fit <- gbm(
    f_brt_exp, data = dat.ml.expanded, distribution = "gaussian",
    n.trees = 1500, interaction.depth = 2, shrinkage = 0.005,
    n.minobsinnode = 5, bag.fraction = 0.70, cv.folds = 5
)
best_exp_t <- gbm.perf(brt_exp_fit, method = "cv", plot.it = FALSE)
p5_exp_in <- predict(brt_exp_fit, newdata = dat.ml.expanded, n.trees = best_exp_t, type = "response")

p5_exp_cv <- numeric(N_exp)
for (f in 1:5) {
    tr_f <- dat.ml.expanded[cv_folds_exp != f, ]
    va_f <- dat.ml.expanded[cv_folds_exp == f, ]
    fit_rf <- gbm(
        f_brt_exp, data = tr_f, distribution = "gaussian", n.trees = 1000,
        interaction.depth = 2, shrinkage = 0.005, n.minobsinnode = 5, bag.fraction = 0.70
    )
    p5_exp_cv[cv_folds_exp == f] <- predict(fit_rf, newdata = va_f, n.trees = best_exp_t, type = "response")
}
met5_exp_rand <- calc_comb_metrics(y_exp, p5_exp_in, p5_exp_cv, SS_tot_exp)

sec_exp <- dat.ml.expanded$Sector
if (is.null(sec_exp) || length(unique(sec_exp)) < 3) sec_exp <- dat.ml.expanded$ReefName
unique_sec_exp <- unique(sec_exp)
p5_exp_sp_cv <- numeric(N_exp)
for (s in unique_sec_exp) {
    val_i <- which(sec_exp == s)
    tr_i  <- which(sec_exp != s)
    if (length(tr_i) > 10) {
        fit_s <- gbm(
            f_brt_exp, data = dat.ml.expanded[tr_i, ], distribution = "gaussian",
            n.trees = 1000, interaction.depth = 2, shrinkage = 0.005,
            n.minobsinnode = 5, bag.fraction = 0.70, verbose = FALSE
        )
        p5_exp_sp_cv[val_i] <- predict(fit_s, newdata = dat.ml.expanded[val_i, ], n.trees = best_exp_t, type = "response")
    } else {
        p5_exp_sp_cv[val_i] <- p5_exp_in[val_i]
    }
}
met5_exp_sp <- calc_comb_metrics(y_exp, p5_exp_in, p5_exp_sp_cv, SS_tot_exp)

cat("Stage 2 completed successfully!\n")

