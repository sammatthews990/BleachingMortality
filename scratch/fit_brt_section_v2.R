library(dplyr)
library(gbm)
library(ggplot2)
library(patchwork)
library(purrr)
library(pROC)

# Load cleaned dataset prepared with strict disturbance filter
source("scratch/check_disturbances_cheung.R") # produces dat.ml.mant.raw

dat.ml.mant.clean <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
    Mort.prop.nudge = (Mort.prop * (n() - 1) + 0.5) / n()
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

cat("dat.ml.mant.clean N =", nrow(dat.ml.mant.clean), "\n")

predictors <- c("MaxDHW.mean", "secc3m", "winyear_mean", "histmDHW6", "mcur_90", "winyear_sd", "cloudp_90", "yrsince6")
pred_labels <- c(
  "MaxDHW.mean" = "Max DHW (°C-weeks)",
  "secc3m" = "Secchi Depth (m)",
  "winyear_mean" = "Winter SST Mean (°C)",
  "histmDHW6" = "Prior Exposure (>6 DHW)",
  "mcur_90" = "Current Speed (90d)",
  "winyear_sd" = "Winter SST SD",
  "cloudp_90" = "Cloud Cover (90d)",
  "yrsince6" = "Years Since >6 DHW"
)

# 1. Fit main 8-predictor BRT model
set.seed(42)
brt_fit <- gbm(
  formula = Mort.prop ~ MaxDHW.mean + secc3m + winyear_mean + histmDHW6 + yrsince6 + mcur_90 + winyear_sd + cloudp_90,
  data = dat.ml.mant.clean,
  distribution = "gaussian",
  n.trees = 2500,
  interaction.depth = 3,
  shrinkage = 0.005,
  bag.fraction = 0.75,
  cv.folds = 5
)

best_iter <- gbm.perf(brt_fit, plot.it = FALSE, method = "cv")
cat("Optimal boosting trees (5-fold CV):", best_iter, "\n")

# 2. Compute Performance & Summary Statistics
null_dev <- sum((dat.ml.mant.clean$Mort.prop - mean(dat.ml.mant.clean$Mort.prop))^2)
res_dev <- sum((dat.ml.mant.clean$Mort.prop - brt_fit$fit)^2)
dev_explained <- (1 - (res_dev / null_dev)) * 100

train_preds <- brt_fit$fit
cv_preds <- brt_fit$cv.fitted

rmse_train <- sqrt(mean((dat.ml.mant.clean$Mort.prop - train_preds)^2))
rmse_cv <- sqrt(mean((dat.ml.mant.clean$Mort.prop - cv_preds)^2))
mae_train <- mean(abs(dat.ml.mant.clean$Mort.prop - train_preds))
mae_cv <- mean(abs(dat.ml.mant.clean$Mort.prop - cv_preds))

r2_train <- cor(dat.ml.mant.clean$Mort.prop, train_preds)^2
r2_cv <- cor(dat.ml.mant.clean$Mort.prop, cv_preds)^2

thresh <- 0.15
actual_binary <- factor(ifelse(dat.ml.mant.clean$Mort.prop >= thresh, "High", "Low"), levels = c("Low", "High"))
pred_binary_train <- factor(ifelse(train_preds >= thresh, "High", "Low"), levels = c("Low", "High"))
pred_binary_cv <- factor(ifelse(cv_preds >= thresh, "High", "Low"), levels = c("Low", "High"))

acc_train <- mean(actual_binary == pred_binary_train)
acc_cv <- mean(actual_binary == pred_binary_cv)

roc_train <- roc(actual_binary, train_preds, quiet = TRUE)
roc_cv <- roc(actual_binary, cv_preds, quiet = TRUE)
auc_train <- as.numeric(auc(roc_train))
auc_cv <- as.numeric(auc(roc_cv))

metrics_df <- data.frame(
  Metric = c(
    "Predictor Count",
    "Deviance Explained (%)",
    "CV R² (Correlation Squared)",
    "Training R²",
    "CV RMSE (Root Mean Sq. Error)",
    "Training RMSE",
    "CV MAE (Mean Abs. Error)",
    "Training MAE",
    "Classification Threshold (Mortality)",
    "CV Classification Accuracy (%)",
    "Training Classification Accuracy (%)",
    "CV AUC (Area Under ROC)",
    "Training AUC"
  ),
  Value = c(
    "8 Environmental Predictors",
    sprintf("%.1f%%", dev_explained),
    sprintf("%.3f", r2_cv),
    sprintf("%.3f", r2_train),
    sprintf("%.4f", rmse_cv),
    sprintf("%.4f", rmse_train),
    sprintf("%.4f", mae_cv),
    sprintf("%.4f", mae_train),
    sprintf("%.0f%%", thresh * 100),
    sprintf("%.1f%%", acc_cv * 100),
    sprintf("%.1f%%", acc_train * 100),
    sprintf("%.3f", auc_cv),
    sprintf("%.3f", auc_train)
  )
)

cat("\n=== Summary Statistics Table (8 Predictors) ===\n")
print(metrics_df)

var_imp <- summary(brt_fit, n.trees = best_iter, plotit = FALSE)
var_imp$Variable_Label <- pred_labels[as.character(var_imp$var)]

# 3. Bootstrap 1D Partial Dependence Plots (80% and 95% CIs, B = 100)
B <- 100
cat("\nRunning", B, "bootstrap iterations for 8 Predictors (80% & 95% CIs)...\n")

grids <- map(predictors, function(p) {
  seq(min(dat.ml.mant.clean[[p]], na.rm = TRUE),
      max(dat.ml.mant.clean[[p]], na.rm = TRUE),
      length.out = 50)
})
names(grids) <- predictors

set.seed(123)
boot_pdp_list <- list()

for (b in 1:B) {
  boot_idx <- sample(seq_len(nrow(dat.ml.mant.clean)), replace = TRUE)
  boot_data <- dat.ml.mant.clean[boot_idx, ]
  
  boot_fit <- gbm(
    formula = Mort.prop ~ MaxDHW.mean + secc3m + winyear_mean + histmDHW6 + yrsince6 + mcur_90 + winyear_sd + cloudp_90,
    data = boot_data,
    distribution = "gaussian",
    n.trees = best_iter,
    interaction.depth = 3,
    shrinkage = 0.005,
    bag.fraction = 0.75,
    verbose = FALSE
  )
  
  for (p in predictors) {
    grid_vals <- grids[[p]]
    pdp_vals <- map_dbl(grid_vals, function(x_val) {
      synth_data <- boot_data
      synth_data[[p]] <- x_val
      mean(predict(boot_fit, newdata = synth_data, n.trees = best_iter))
    })
    
    boot_pdp_list[[length(boot_pdp_list) + 1]] <- data.frame(
      boot = b,
      predictor = p,
      x = grid_vals,
      yhat = pdp_vals
    )
  }
}

boot_pdp_df <- bind_rows(boot_pdp_list)

pdp_summary <- boot_pdp_df %>%
  group_by(predictor, x) %>%
  summarise(
    median_yhat = median(yhat),
    lo_95 = quantile(yhat, 0.025),
    hi_95 = quantile(yhat, 0.975),
    lo_80 = quantile(yhat, 0.100),
    hi_80 = quantile(yhat, 0.900),
    .groups = "drop"
  ) %>%
  mutate(pred_label = pred_labels[predictor])

# 4. Generate 2-Way Partial Dependence Interactions
cat("\nComputing 2-Way Interaction PDPs...\n")

# Interaction 1: MaxDHW.mean x secc3m (DHW x Water Clarity)
grid_dhw_sec <- expand.grid(
  MaxDHW.mean = seq(min(dat.ml.mant.clean$MaxDHW.mean), max(dat.ml.mant.clean$MaxDHW.mean), length.out = 40),
  secc3m = seq(min(dat.ml.mant.clean$secc3m), max(dat.ml.mant.clean$secc3m), length.out = 40)
) %>% mutate(
  winyear_mean = median(dat.ml.mant.clean$winyear_mean),
  cloudp_90 = median(dat.ml.mant.clean$cloudp_90),
  histmDHW6 = median(dat.ml.mant.clean$histmDHW6),
  yrsince6 = median(dat.ml.mant.clean$yrsince6),
  mcur_90 = median(dat.ml.mant.clean$mcur_90),
  winyear_sd = median(dat.ml.mant.clean$winyear_sd)
)
grid_dhw_sec$pred <- predict(brt_fit, newdata = grid_dhw_sec, n.trees = best_iter)

# Interaction 2: MaxDHW.mean x winyear_mean (DHW x Winter SST Mean)
grid_dhw_winm <- expand.grid(
  MaxDHW.mean = seq(min(dat.ml.mant.clean$MaxDHW.mean), max(dat.ml.mant.clean$MaxDHW.mean), length.out = 40),
  winyear_mean = seq(min(dat.ml.mant.clean$winyear_mean), max(dat.ml.mant.clean$winyear_mean), length.out = 40)
) %>% mutate(
  secc3m = median(dat.ml.mant.clean$secc3m),
  cloudp_90 = median(dat.ml.mant.clean$cloudp_90),
  histmDHW6 = median(dat.ml.mant.clean$histmDHW6),
  yrsince6 = median(dat.ml.mant.clean$yrsince6),
  mcur_90 = median(dat.ml.mant.clean$mcur_90),
  winyear_sd = median(dat.ml.mant.clean$winyear_sd)
)
grid_dhw_winm$pred <- predict(brt_fit, newdata = grid_dhw_winm, n.trees = best_iter)

# Save cached objects
saveRDS(brt_fit, "output/models/brt_manta_fit.rds")
saveRDS(var_imp, "output/models/brt_manta_var_imp.rds")
saveRDS(pdp_summary, "output/models/brt_manta_pdp_boot.rds")
saveRDS(metrics_df, "output/models/brt_manta_metrics.rds")
saveRDS(grid_dhw_sec, "output/models/brt_manta_pdp_2way_sec.rds")
saveRDS(grid_dhw_winm, "output/models/brt_manta_pdp_2way_winm.rds")

cat("=== BRT fitting with 8 predictors complete! ===\n")
