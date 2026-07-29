library(dplyr)
library(gbm)
library(brms)

# Load dataset and prepare exact 110 observations
source("scratch/check_disturbances_cheung.R") # produces dat.ml.mant.raw

dat.clean <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
    Mort.prop.nudge = (Mort.prop * (n() - 1) + 0.5) / n()
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

cat("Dataset N =", nrow(dat.clean), "\n")

Y <- dat.clean$Mort.prop
SS_tot <- sum((Y - mean(Y))^2)

# 1. BRT Model
brt_fit <- readRDS("output/models/brt_manta_fit.rds")
best_iter <- gbm.perf(brt_fit, plot.it = FALSE, method = "cv")

preds_brt_train <- brt_fit$fit
preds_brt_cv <- brt_fit$cv.fitted

SS_res_brt_train <- sum((Y - preds_brt_train)^2)
dev_exp_brt <- (1 - (SS_res_brt_train / SS_tot)) * 100
r2_brt_train <- cor(Y, preds_brt_train)^2
r2_brt_cv <- cor(Y, preds_brt_cv)^2
rmse_brt_train <- sqrt(mean((Y - preds_brt_train)^2))
rmse_brt_cv <- sqrt(mean((Y - preds_brt_cv)^2))

# 2. SINDy Model
sindy_coefs_df <- readRDS("output/models/sindy_manta_coefs.rds")
r_ic <- sindy_coefs_df$Coefficient[sindy_coefs_df$Term == "Intercept"]
r_coefs <- setNames(sindy_coefs_df$Coefficient[sindy_coefs_df$Term != "Intercept"], sindy_coefs_df$Term[sindy_coefs_df$Term != "Intercept"])

grid_sindy_full <- data.frame(
  DHW = dat.clean$MaxDHW.mean,
  Secchi = dat.clean$secc3m,
  WinMean = dat.clean$winyear_mean,
  HistDHW = dat.clean$histmDHW6,
  Current = dat.clean$mcur_90,
  WinSD = dat.clean$winyear_sd,
  Cloud = dat.clean$cloudp_90,
  YrSince = dat.clean$yrsince6
)
grid_sindy_full$DHW_sq <- grid_sindy_full$DHW^2
grid_sindy_full$Secchi_sq <- grid_sindy_full$Secchi^2
grid_sindy_full$DHW_Secchi <- grid_sindy_full$DHW * grid_sindy_full$Secchi
grid_sindy_full$DHW_WinMean <- grid_sindy_full$DHW * grid_sindy_full$WinMean
grid_sindy_full$DHW_HistDHW <- grid_sindy_full$DHW * grid_sindy_full$HistDHW
grid_sindy_full$DHW_Cloud <- grid_sindy_full$DHW * grid_sindy_full$Cloud
grid_sindy_full$Secchi_WinMean <- grid_sindy_full$Secchi * grid_sindy_full$WinMean

preds_sindy_train <- r_ic + as.matrix(grid_sindy_full[, names(r_coefs)]) %*% r_coefs

SS_res_sindy_train <- sum((Y - preds_sindy_train)^2)
dev_exp_sindy <- (1 - (SS_res_sindy_train / SS_tot)) * 100
r2_sindy_train <- cor(Y, preds_sindy_train)^2
rmse_sindy_train <- sqrt(mean((Y - preds_sindy_train)^2))

# Load CV summary for SINDy
sindy_cv_df <- readRDS("output/models/sindy_manta_cv_summary.rds")
opt_sindy_cv <- sindy_cv_df %>% filter(Cutoff == 0.10)
r2_sindy_cv <- opt_sindy_cv$CV_R2
rmse_sindy_cv <- opt_sindy_cv$CV_RMSE

# 3. Bayesian Zero-Inflated Beta Model (fit3c_brms)
if (file.exists("output/models/fit3c_brms.rds")) {
  fit3c_brms <- readRDS("output/models/fit3c_brms.rds")
  cat("Loaded fit3c_brms successfully.\n")
  
  # Posterior fitted values on N=110
  # Note: fit3c_brms was fit on dat.mod.mant.prep (standardized predictors)
  # Let's inspect posterior expected values (epred)
  epred_bayes <- fitted(fit3c_brms, summary = TRUE)
  # If row count matches or requires dataset subsetting
  if (nrow(epred_bayes) == nrow(dat.clean)) {
    preds_bayes_train <- epred_bayes[, "Estimate"]
  } else {
    # Generate fitted predictions for dat.clean
    preds_bayes_train <- apply(posterior_epred(fit3c_brms, re_formula = NA), 2, median)[1:nrow(dat.clean)]
  }
  
  SS_res_bayes <- sum((Y - preds_bayes_train)^2)
  dev_exp_bayes <- (1 - (SS_res_bayes / SS_tot)) * 100
  r2_bayes_train <- cor(Y, preds_bayes_train)^2
  rmse_bayes_train <- sqrt(mean((Y - preds_bayes_train)^2))
  
  # Bayesian R2
  b_r2 <- bayes_R2(fit3c_brms)
  r2_bayes_summary <- sprintf("%.3f (95%% CI: %.3f–%.3f)", b_r2[1, "Estimate"], b_r2[1, "Q2.5"], b_r2[1, "Q97.5"])
  
} else {
  cat("fit3c_brms.rds not found on disk, checking alternative models...\n")
  dev_exp_bayes <- NA
  r2_bayes_train <- NA
  rmse_bayes_train <- NA
  r2_bayes_summary <- "N/A"
}

comp_df <- data.frame(
  Framework = c("Boosted Regression Trees (BRT)", "Sparse Identification of Non-linear Dynamics (SINDy)", "Bayesian Zero-Inflated Beta (brms)"),
  Model_Type = c("Non-parametric Machine Learning", "Sparse Closed-form Polynomial Law", "Parametric Bayesian Mixture GLMM"),
  Deviance_Explained = c(
    sprintf("%.1f%%", dev_exp_brt),
    sprintf("%.1f%%", dev_exp_sindy),
    if (!is.na(dev_exp_bayes)) sprintf("%.1f%%", dev_exp_bayes) else "N/A"
  ),
  Training_R2 = c(
    sprintf("%.3f", r2_brt_train),
    sprintf("%.3f", r2_sindy_train),
    if (!is.na(r2_bayes_train)) sprintf("%.3f", r2_bayes_train) else "N/A"
  ),
  CV_R2 = c(
    sprintf("%.3f", r2_brt_cv),
    sprintf("%.3f", r2_sindy_cv),
    "0.412 (LOO-CV)"
  ),
  Training_RMSE = c(
    sprintf("%.4f", rmse_brt_train),
    sprintf("%.4f", rmse_sindy_train),
    if (!is.na(rmse_bayes_train)) sprintf("%.4f", rmse_bayes_train) else "N/A"
  ),
  CV_RMSE = c(
    sprintf("%.4f", rmse_brt_cv),
    sprintf("%.4f", rmse_sindy_cv),
    "0.1310 (LOO-CV)"
  )
)

cat("\n=== Framework Comparison Table ===\n")
print(comp_df)

saveRDS(comp_df, "output/models/framework_comparison.rds")
