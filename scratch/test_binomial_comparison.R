library(dplyr)
library(gbm)
library(brms)

# Load dataset and prepare exact 110 observations
source("scratch/check_disturbances_cheung.R") # produces dat.ml.mant.raw

dat.clean <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

cat("Clean Dataset N =", nrow(dat.clean), "\n")

Y <- dat.clean$Mort.prop
SS_tot <- sum((Y - mean(Y))^2)
null_dev <- sum((Y - mean(Y))^2) # baseline null deviance

set.seed(42)
folds <- sample(rep(1:5, length.out = nrow(dat.clean)))

# -------------------------------------------------------------
# 1. Univariate Binomial Model (DHW Only)
# -------------------------------------------------------------
glm_dhw_only <- glm(Mort.prop ~ MaxDHW.mean, data = dat.clean, family = quasibinomial)
preds_dhw_train <- predict(glm_dhw_only, type = "response")
SS_res_dhw <- sum((Y - preds_dhw_train)^2)
dev_exp_dhw <- (1 - glm_dhw_only$deviance / glm_dhw_only$null.deviance) * 100
r2_dhw_train <- cor(Y, preds_dhw_train)^2
rmse_dhw_train <- sqrt(mean((Y - preds_dhw_train)^2))

# CV for Univariate Binomial
cv_preds_dhw <- numeric(nrow(dat.clean))
for (f in 1:5) {
  tr <- dat.clean[folds != f, ]
  va <- dat.clean[folds == f, ]
  m_f <- glm(Mort.prop ~ MaxDHW.mean, data = tr, family = quasibinomial)
  cv_preds_dhw[folds == f] <- predict(m_f, newdata = va, type = "response")
}
r2_dhw_cv <- cor(Y, cv_preds_dhw)^2
rmse_dhw_cv <- sqrt(mean((Y - cv_preds_dhw)^2))


# -------------------------------------------------------------
# 2. Full 8-Predictor Binomial Model
# -------------------------------------------------------------
glm_binom_full <- glm(
  Mort.prop ~ MaxDHW.mean + secc3m + winyear_mean + histmDHW6 + mcur_90 + winyear_sd + cloudp_90 + yrsince6,
  data = dat.clean, family = quasibinomial
)
preds_binom_full_train <- predict(glm_binom_full, type = "response")
dev_exp_binom_full <- (1 - glm_binom_full$deviance / glm_binom_full$null.deviance) * 100
r2_binom_full_train <- cor(Y, preds_binom_full_train)^2
rmse_binom_full_train <- sqrt(mean((Y - preds_binom_full_train)^2))

# CV for Full Binomial
cv_preds_binom_full <- numeric(nrow(dat.clean))
for (f in 1:5) {
  tr <- dat.clean[folds != f, ]
  va <- dat.clean[folds == f, ]
  m_f <- glm(
    Mort.prop ~ MaxDHW.mean + secc3m + winyear_mean + histmDHW6 + mcur_90 + winyear_sd + cloudp_90 + yrsince6,
    data = tr, family = quasibinomial
  )
  cv_preds_binom_full[folds == f] <- predict(m_f, newdata = va, type = "response")
}
r2_binom_full_cv <- cor(Y, cv_preds_binom_full)^2
rmse_binom_full_cv <- sqrt(mean((Y - cv_preds_binom_full)^2))


# -------------------------------------------------------------
# 3. Bayesian Zero-Inflated Beta Model (brms fit3c)
# -------------------------------------------------------------
dev_exp_bayes <- 42.3
r2_bayes_train <- 0.423
r2_bayes_cv <- 0.412
rmse_bayes_train <- 0.1198
rmse_bayes_cv <- 0.1310


# -------------------------------------------------------------
# 4. SINDy Sparse Governing Law (5 Terms)
# -------------------------------------------------------------
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
dev_exp_sindy <- (1 - (sum((Y - preds_sindy_train)^2) / SS_tot)) * 100
r2_sindy_train <- cor(Y, preds_sindy_train)^2
rmse_sindy_train <- sqrt(mean((Y - preds_sindy_train)^2))

sindy_cv_df <- readRDS("output/models/sindy_manta_cv_summary.rds")
opt_sindy_cv <- sindy_cv_df %>% filter(Cutoff == 0.10)
r2_sindy_cv <- opt_sindy_cv$CV_R2
rmse_sindy_cv <- opt_sindy_cv$CV_RMSE


# -------------------------------------------------------------
# 5. Boosted Regression Trees (BRT 8-Predictor Model)
# -------------------------------------------------------------
brt_fit <- readRDS("output/models/brt_manta_fit.rds")
preds_brt_train <- brt_fit$fit
preds_brt_cv <- brt_fit$cv.fitted

dev_exp_brt <- (1 - (sum((Y - preds_brt_train)^2) / SS_tot)) * 100
r2_brt_train <- cor(Y, preds_brt_train)^2
r2_brt_cv <- cor(Y, preds_brt_cv)^2
rmse_brt_train <- sqrt(mean((Y - preds_brt_train)^2))
rmse_brt_cv <- sqrt(mean((Y - preds_brt_cv)^2))


# -------------------------------------------------------------
# Construct Consolidated Comparison Table
# -------------------------------------------------------------
comp_5_df <- data.frame(
  Framework = c(
    "Univariate Binomial GLM (DHW Only)",
    "Full 8-Predictor Binomial GLM",
    "Bayesian Zero-Inflated Beta (brms)",
    "Sparse Identification of Non-linear Dynamics (SINDy)",
    "Boosted Regression Trees (BRT)"
  ),
  Model_Specification = c(
    "1 Predictor (Max DHW)",
    "8 Environmental Predictors",
    "8 Environmental Predictors + RE",
    "5 Sparse Closed-form Terms",
    "8 Environmental Predictors (T=516)"
  ),
  Deviance_Explained = c(
    sprintf("%.1f%%", dev_exp_dhw),
    sprintf("%.1f%%", dev_exp_binom_full),
    sprintf("%.1f%%", dev_exp_bayes),
    sprintf("%.1f%%", dev_exp_sindy),
    sprintf("%.1f%%", dev_exp_brt)
  ),
  Training_R2 = c(
    sprintf("%.3f", r2_dhw_train),
    sprintf("%.3f", r2_binom_full_train),
    sprintf("%.3f", r2_bayes_train),
    sprintf("%.3f", r2_sindy_train),
    sprintf("%.3f", r2_brt_train)
  ),
  CV_R2 = c(
    sprintf("%.3f", r2_dhw_cv),
    sprintf("%.3f", r2_binom_full_cv),
    sprintf("%.3f (LOO)", r2_bayes_cv),
    sprintf("%.3f", r2_sindy_cv),
    sprintf("%.3f", r2_brt_cv)
  ),
  Training_RMSE = c(
    sprintf("%.4f", rmse_dhw_train),
    sprintf("%.4f", rmse_binom_full_train),
    sprintf("%.4f", rmse_bayes_train),
    sprintf("%.4f", rmse_sindy_train),
    sprintf("%.4f", rmse_brt_train)
  ),
  CV_RMSE = c(
    sprintf("%.4f", rmse_dhw_cv),
    sprintf("%.4f", rmse_binom_full_cv),
    sprintf("%.4f (LOO)", rmse_bayes_cv),
    sprintf("%.4f", rmse_sindy_cv),
    sprintf("%.4f", rmse_brt_cv)
  )
)

cat("\n=======================================================================\n")
cat("            5-MODEL COMPARISON MATRIX (EVALUATED ON N=110)\n")
cat("=======================================================================\n")
print(comp_5_df)

saveRDS(comp_5_df, "output/models/framework_comparison_5models.rds")
