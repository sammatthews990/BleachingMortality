library(dplyr)
library(gbm)
library(brms)

source("scratch/check_disturbances_cheung.R")

# Recreate dat.ml.ltmp exactly as done in DHW_Mortality_Attribution.qmd
cheung_predictors <- dat.cheung.sum %>%
  ungroup() %>%
  mutate(LABEL_id = toupper(trimws(LABEL_id)))

dat.ml.ltmp <- dat.mod.bent.hc %>%
  ungroup() %>%
  filter(project_code == "LTMP") %>%
  filter(!is.na(ReefID), !is.na(DHWYear)) %>%
  mutate(DHWYear = as.numeric(DHWYear)) %>%
  left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>%
  filter(!is.na(mcur_90)) %>%
  mutate(survey_type = "ltmp_benthic")

# Clean LTMP dataset (N=107)
acrop_comp <- df.AIMS.full %>%
  filter(
    data_type == "photo-transect", domain_category == "reef",
    purpose == "COMPOSITION", variable == "HARD CORAL"
  ) %>%
  mutate(Reef_Name_Clean = gsub(" REEF(S)?| ISLAND| IS", "", toupper(domain_name))) %>%
  left_join(dat.Ref.Clean %>% select(AIMS_REEF_NAME_Clean, ReefName),
    by = c("Reef_Name_Clean" = "AIMS_REEF_NAME_Clean")
  ) %>%
  filter(!is.na(ReefName)) %>%
  group_by(ReefName, report_year, depth) %>%
  summarise(
    total_hc = sum(mean, na.rm = TRUE),
    acrop_cover = sum(mean[grepl("Acropora", reefpage_category, ignore.case = TRUE)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(prop_acropora = ifelse(total_hc > 0, acrop_cover / total_hc, 0))

dat.ml.ltmp.clean <- dat.ml.ltmp %>%
  filter(!is.na(mcur_90), !is.na(Rel.Change), !is.na(MaxDHW.mean)) %>%
  left_join(acrop_comp %>% select(ReefName, report_year, depth, prop_acropora),
    by = c("ReefName", "report_year", "depth")
  ) %>%
  mutate(prop_acropora = replace_na(prop_acropora, median(prop_acropora, na.rm = TRUE))) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  ) %>%
  filter(is.finite(Mort.prop))

cat("Independent LTMP Benthic Test Dataset N =", nrow(dat.ml.ltmp.clean), "\n")

# Training Manta Tow Dataset (N=110)
dat.train <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

cat("Training Manta Tow Dataset N =", nrow(dat.train), "\n")

dat.test <- dat.ml.ltmp.clean
Y_test <- dat.test$Mort.prop
SS_tot_test <- sum((Y_test - mean(Y_test))^2)

eval_metrics <- function(y_true, y_pred, model_name) {
  y_pred_b <- pmin(pmax(y_pred, 0), 1)
  SS_res <- sum((y_true - y_pred_b)^2)
  r2 <- max(0, cor(y_true, y_pred_b)^2)
  dev_exp <- (1 - (SS_res / SS_tot_test)) * 100
  rmse <- sqrt(mean((y_true - y_pred_b)^2))
  
  data.frame(
    Framework = model_name,
    Test_Deviance_Explained = sprintf("%.1f%%", dev_exp),
    Test_R2 = sprintf("%.3f", r2),
    Test_RMSE = sprintf("%.4f", rmse)
  )
}

# 1. Univariate Binomial (DHW Only)
m1_dhw <- glm(Mort.prop ~ MaxDHW.mean, data = dat.train, family = quasibinomial)
preds_m1 <- predict(m1_dhw, newdata = dat.test, type = "response")
res1 <- eval_metrics(Y_test, preds_m1, "Univariate Binomial GLM (DHW Only)")

# 2. Full 8-Predictor Binomial
m2_binom <- glm(
  Mort.prop ~ MaxDHW.mean + secc3m + winyear_mean + histmDHW6 + mcur_90 + winyear_sd + cloudp_90 + yrsince6,
  data = dat.train, family = quasibinomial
)
preds_m2 <- predict(m2_binom, newdata = dat.test, type = "response")
res2 <- eval_metrics(Y_test, preds_m2, "Full 8-Predictor Binomial GLM")

# 3. Bayesian Zero-Inflated Beta Model (brms fit3c)
if (file.exists("output/models/fit3c_brms.rds")) {
  fit3c <- readRDS("output/models/fit3c_brms.rds")
  # Standardize test predictors using training scaling params
  scale_params <- list(
    MaxDHW.mean = c(mean = mean(dat.train$MaxDHW.mean), sd2 = 2 * sd(dat.train$MaxDHW.mean)),
    secc3m = c(mean = mean(dat.train$secc3m), sd2 = 2 * sd(dat.train$secc3m)),
    winyear_mean = c(mean = mean(dat.train$winyear_mean), sd2 = 2 * sd(dat.train$winyear_mean)),
    histmDHW6 = c(mean = mean(dat.train$histmDHW6), sd2 = 2 * sd(dat.train$histmDHW6)),
    mcur_90 = c(mean = mean(dat.train$mcur_90), sd2 = 2 * sd(dat.train$mcur_90)),
    winyear_sd = c(mean = mean(dat.train$winyear_sd), sd2 = 2 * sd(dat.train$winyear_sd)),
    cloudp_90 = c(mean = mean(dat.train$cloudp_90), sd2 = 2 * sd(dat.train$cloudp_90)),
    yrsince6 = c(mean = mean(dat.train$yrsince6), sd2 = 2 * sd(dat.train$yrsince6))
  )
  
  grid_test_brms <- data.frame(
    MaxDHW.mean_s = (dat.test$MaxDHW.mean - scale_params$MaxDHW.mean["mean"]) / scale_params$MaxDHW.mean["sd2"],
    secc3m_s = (dat.test$secc3m - scale_params$secc3m["mean"]) / scale_params$secc3m["sd2"],
    winyear_mean_s = (dat.test$winyear_mean - scale_params$winyear_mean["mean"]) / scale_params$winyear_mean["sd2"],
    histmDHW6_s = (dat.test$histmDHW6 - scale_params$histmDHW6["mean"]) / scale_params$histmDHW6["sd2"],
    mcur_90_s = (dat.test$mcur_90 - scale_params$mcur_90["mean"]) / scale_params$mcur_90["sd2"],
    winyear_sd_s = (dat.test$winyear_sd - scale_params$winyear_sd["mean"]) / scale_params$winyear_sd["sd2"],
    cloudp_90_s = (dat.test$cloudp_90 - scale_params$cloudp_90["mean"]) / scale_params$cloudp_90["sd2"],
    yrsince6_s = (dat.test$yrsince6 - scale_params$yrsince6["mean"]) / scale_params$yrsince6["sd2"]
  )
  
  preds_m3 <- apply(posterior_epred(fit3c, newdata = grid_test_brms, re_formula = NA), 2, median)
} else {
  m3_beta <- glm(
    (Mort.prop * (nrow(dat.train) - 1) + 0.5) / nrow(dat.train) ~ MaxDHW.mean + secc3m + winyear_mean + histmDHW6 + mcur_90 + winyear_sd + cloudp_90 + yrsince6,
    data = dat.train, family = quasibinomial
  )
  preds_m3 <- predict(m3_beta, newdata = dat.test, type = "response")
}
res3 <- eval_metrics(Y_test, preds_m3, "Bayesian Zero-Inflated Beta (brms)")

# 4. SINDy Sparse Governing Law
sindy_coefs_df <- readRDS("output/models/sindy_manta_coefs.rds")
r_ic <- sindy_coefs_df$Coefficient[sindy_coefs_df$Term == "Intercept"]
r_coefs <- setNames(sindy_coefs_df$Coefficient[sindy_coefs_df$Term != "Intercept"], sindy_coefs_df$Term[sindy_coefs_df$Term != "Intercept"])

grid_sindy_test <- data.frame(
  DHW = dat.test$MaxDHW.mean,
  Secchi = dat.test$secc3m,
  WinMean = dat.test$winyear_mean,
  HistDHW = dat.test$histmDHW6,
  Current = dat.test$mcur_90,
  WinSD = dat.test$winyear_sd,
  Cloud = dat.test$cloudp_90,
  YrSince = dat.test$yrsince6
)
grid_sindy_test$DHW_sq <- grid_sindy_test$DHW^2
grid_sindy_test$Secchi_sq <- grid_sindy_test$Secchi^2
grid_sindy_test$DHW_Secchi <- grid_sindy_test$DHW * grid_sindy_test$Secchi
grid_sindy_test$DHW_WinMean <- grid_sindy_test$DHW * grid_sindy_test$WinMean
grid_sindy_test$DHW_HistDHW <- grid_sindy_test$DHW * grid_sindy_test$HistDHW
grid_sindy_test$DHW_Cloud <- grid_sindy_test$DHW * grid_sindy_test$Cloud
grid_sindy_test$Secchi_WinMean <- grid_sindy_test$Secchi * grid_sindy_test$WinMean

preds_m4 <- r_ic + as.matrix(grid_sindy_test[, names(r_coefs)]) %*% r_coefs
res4 <- eval_metrics(Y_test, preds_m4, "Sparse Identification of Non-linear Dynamics (SINDy)")

# 5. Boosted Regression Trees (BRT)
brt_fit <- readRDS("output/models/brt_manta_fit.rds")
best_iter <- gbm.perf(brt_fit, plot.it = FALSE, method = "cv")

test_df_brt <- data.frame(
  MaxDHW.mean = dat.test$MaxDHW.mean,
  secc3m = dat.test$secc3m,
  winyear_mean = dat.test$winyear_mean,
  histmDHW6 = dat.test$histmDHW6,
  mcur_90 = dat.test$mcur_90,
  winyear_sd = dat.test$winyear_sd,
  cloudp_90 = dat.test$cloudp_90,
  yrsince6 = dat.test$yrsince6
)

preds_m5 <- predict(brt_fit, newdata = test_df_brt, n.trees = best_iter, type = "response")
res5 <- eval_metrics(Y_test, preds_m5, "Boosted Regression Trees (BRT)")

# Combine into summary table
ltmp_eval_df <- bind_rows(res1, res2, res3, res4, res5)

cat("\n=======================================================================\n")
cat("      INDEPENDENT TEST SET EVALUATION ON LTMP BENTHIC DATASET (N=107)\n")
cat("=======================================================================\n")
print(ltmp_eval_df)

saveRDS(ltmp_eval_df, "output/models/ltmp_independent_test_metrics.rds")
