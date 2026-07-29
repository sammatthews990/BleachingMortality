library(dplyr)
library(gbm)
library(pROC)

source("scratch/check_disturbances_cheung.R") # produces dat.ml.mant.raw

dat.clean <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

cat("dat.clean N =", nrow(dat.clean), "\n")

preds_6 <- c("MaxDHW.mean", "secc3m", "cloudp_90", "histmDHW6", "mcur_90", "winyear_sd")
preds_8 <- c("MaxDHW.mean", "secc3m", "cloudp_90", "histmDHW6", "yrsince6", "mcur_90", "winyear_sd", "winyear_mean")

# Check NA counts for 8 predictors
cat("\nNA counts in 8 predictors:\n")
sapply(dat.clean[preds_8], function(x) sum(is.na(x)))

# Filter dataset to complete cases on 8 predictors if any
dat.clean.8 <- dat.clean %>% filter(if_all(all_of(preds_8), is.finite))
cat("dat.clean.8 N =", nrow(dat.clean.8), "\n")

# Fit 6-predictor model
set.seed(42)
brt_6 <- gbm(
  formula = Mort.prop ~ MaxDHW.mean + secc3m + cloudp_90 + histmDHW6 + mcur_90 + winyear_sd,
  data = dat.clean.8,
  distribution = "gaussian",
  n.trees = 2500,
  interaction.depth = 3,
  shrinkage = 0.005,
  bag.fraction = 0.75,
  cv.folds = 5
)
iter_6 <- gbm.perf(brt_6, plot.it = FALSE, method = "cv")

# Fit 8-predictor model
set.seed(42)
brt_8 <- gbm(
  formula = Mort.prop ~ MaxDHW.mean + secc3m + cloudp_90 + histmDHW6 + yrsince6 + mcur_90 + winyear_sd + winyear_mean,
  data = dat.clean.8,
  distribution = "gaussian",
  n.trees = 2500,
  interaction.depth = 3,
  shrinkage = 0.005,
  bag.fraction = 0.75,
  cv.folds = 5
)
iter_8 <- gbm.perf(brt_8, plot.it = FALSE, method = "cv")

# Diagnostics helper
calc_metrics <- function(mod, d, iter) {
  null_dev <- sum((d$Mort.prop - mean(d$Mort.prop))^2)
  res_dev <- sum((d$Mort.prop - mod$fit)^2)
  dev_exp <- (1 - (res_dev / null_dev)) * 100
  
  train_preds <- mod$fit
  cv_preds <- mod$cv.fitted
  
  rmse_train <- sqrt(mean((d$Mort.prop - train_preds)^2))
  rmse_cv <- sqrt(mean((d$Mort.prop - cv_preds)^2))
  
  r2_train <- cor(d$Mort.prop, train_preds)^2
  r2_cv <- cor(d$Mort.prop, cv_preds)^2
  
  actual_bin <- factor(ifelse(d$Mort.prop >= 0.15, "High", "Low"), levels = c("Low", "High"))
  auc_cv <- as.numeric(auc(roc(actual_bin, cv_preds, quiet = TRUE)))
  acc_cv <- mean(actual_bin == factor(ifelse(cv_preds >= 0.15, "High", "Low"), levels = c("Low", "High")))
  
  data.frame(
    Optimal_Trees = iter,
    Deviance_Explained = sprintf("%.1f%%", dev_exp),
    CV_R2 = sprintf("%.3f", r2_cv),
    CV_RMSE = sprintf("%.4f", rmse_cv),
    Train_RMSE = sprintf("%.4f", rmse_train),
    CV_Accuracy = sprintf("%.1f%%", acc_cv * 100),
    CV_AUC = sprintf("%.3f", auc_cv)
  )
}

res_cmp <- bind_rows(
  "6 Predictors" = calc_metrics(brt_6, dat.clean.8, iter_6),
  "8 Predictors (+yrsince6, +winyear_mean)" = calc_metrics(brt_8, dat.clean.8, iter_8),
  .id = "Model"
)

cat("\n=== Comparison Results ===\n")
print(res_cmp)

cat("\n=== Variable Importance for 8 Predictors ===\n")
print(summary(brt_8, n.trees = iter_8, plotit = FALSE))
