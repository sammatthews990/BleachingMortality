library(dplyr)
library(gbm)

source("scratch/check_disturbances_cheung.R")

# Clean Manta Tow Training Dataset (N=110)
dat.train <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

cat("Manta Training Dataset N =", nrow(dat.train), "\n")

predictor_vars <- c("MaxDHW.mean", "secc3m", "winyear_mean", "histmDHW6", "mcur_90", "winyear_sd", "cloudp_90", "yrsince6")
f_str <- as.formula(paste("Mort.prop ~", paste(predictor_vars, collapse = " + ")))

# -------------------------------------------------------------
# 1. Standard (Unregularized) BRT Baseline (Depth 3, lr 0.01, minobs 5)
# -------------------------------------------------------------
set.seed(42)
brt_orig <- gbm(
  f_str,
  data = dat.train,
  distribution = "gaussian",
  n.trees = 1000,
  interaction.depth = 3,
  shrinkage = 0.01,
  n.minobsinnode = 5,
  bag.fraction = 0.5,
  cv.folds = 5
)
best_iter_orig <- gbm.perf(brt_orig, method = "cv", plot.it = FALSE)

# -------------------------------------------------------------
# 2. Regularized Pairwise BRT (Restricted Depth = 2, Lower Shrinkage, Higher MinObs, Bagging 0.7)
# -------------------------------------------------------------
set.seed(42)
brt_reg2 <- gbm(
  f_str,
  data = dat.train,
  distribution = "gaussian",
  n.trees = 1500,
  interaction.depth = 2,      # Restricted depth = 2 (pairwise interactions only)
  shrinkage = 0.005,          # Lower learning rate
  n.minobsinnode = 10,        # Higher min node size (smoother leaf values)
  bag.fraction = 0.70,        # Higher subsample fraction
  cv.folds = 5
)
best_iter_reg2 <- gbm.perf(brt_reg2, method = "cv", plot.it = FALSE)

# -------------------------------------------------------------
# 3. Additive Regularized BRT (Restricted Depth = 1, Additive Stumps)
# -------------------------------------------------------------
set.seed(42)
brt_reg1 <- gbm(
  f_str,
  data = dat.train,
  distribution = "gaussian",
  n.trees = 1500,
  interaction.depth = 1,      # Pure decision stumps (strictly main effects)
  shrinkage = 0.005,
  n.minobsinnode = 10,
  bag.fraction = 0.70,
  cv.folds = 5
)
best_iter_reg1 <- gbm.perf(brt_reg1, method = "cv", plot.it = FALSE)

# -------------------------------------------------------------
# 4. Spatial / Sector-Blocked Cross-Validation Evaluation
# -------------------------------------------------------------
unique_sec <- unique(dat.train$SECTOR)
set.seed(42)
sec_fold_map <- data.frame(
  SECTOR = unique_sec,
  spatial_fold = sample(rep(1:5, length.out = length(unique_sec)))
)
dat.train <- dat.train %>% left_join(sec_fold_map, by = "SECTOR")

spatial_cv_evaluate <- function(depth, shrinkage, minobs, bag, n_trees_eval = 400) {
  preds_sp <- numeric(nrow(dat.train))
  for (f in 1:5) {
    tr_f <- dat.train %>% filter(spatial_fold != f)
    va_f <- dat.train %>% filter(spatial_fold == f)
    
    fit_f <- gbm(
      f_str,
      data = tr_f,
      distribution = "gaussian",
      n.trees = 1000,
      interaction.depth = depth,
      shrinkage = shrinkage,
      n.minobsinnode = minobs,
      bag.fraction = bag
    )
    preds_sp[dat.train$spatial_fold == f] <- predict(fit_f, newdata = va_f, n.trees = n_trees_eval, type = "response")
  }
  
  y_tr <- dat.train$Mort.prop
  ss_res <- sum((y_tr - preds_sp)^2)
  ss_tot <- sum((y_tr - mean(y_tr))^2)
  sp_r2 <- max(0, cor(y_tr, preds_sp)^2)
  sp_dev <- (1 - ss_res / ss_tot) * 100
  sp_rmse <- sqrt(mean((y_tr - preds_sp)^2))
  
  return(c(Dev = sprintf("%.1f%%", sp_dev), R2 = sprintf("%.3f", sp_r2), RMSE = sprintf("%.4f", sp_rmse)))
}

sp_orig <- spatial_cv_evaluate(depth = 3, shrinkage = 0.01, minobs = 5, bag = 0.5, n_trees_eval = best_iter_orig)
sp_reg2 <- spatial_cv_evaluate(depth = 2, shrinkage = 0.005, minobs = 10, bag = 0.7, n_trees_eval = best_iter_reg2)
sp_reg1 <- spatial_cv_evaluate(depth = 1, shrinkage = 0.005, minobs = 10, bag = 0.7, n_trees_eval = best_iter_reg1)

# Evaluation summary table
eval_summary <- data.frame(
  Model = c(
    "Original BRT (Depth=3, lr=0.01, minobs=5)",
    "Regularized Pairwise BRT (Depth=2, lr=0.005, minobs=10)",
    "Regularized Additive BRT (Depth=1, lr=0.005, minobs=10)"
  ),
  Optimal_Trees = c(best_iter_orig, best_iter_reg2, best_iter_reg1),
  Train_Deviance = c(
    sprintf("%.1f%%", (1 - sum((dat.train$Mort.prop - predict(brt_orig, dat.train, n.trees=best_iter_orig))^2)/sum((dat.train$Mort.prop - mean(dat.train$Mort.prop))^2))*100),
    sprintf("%.1f%%", (1 - sum((dat.train$Mort.prop - predict(brt_reg2, dat.train, n.trees=best_iter_reg2))^2)/sum((dat.train$Mort.prop - mean(dat.train$Mort.prop))^2))*100),
    sprintf("%.1f%%", (1 - sum((dat.train$Mort.prop - predict(brt_reg1, dat.train, n.trees=best_iter_reg1))^2)/sum((dat.train$Mort.prop - mean(dat.train$Mort.prop))^2))*100)
  ),
  Standard_CV_R2 = c(
    sprintf("%.3f", cor(dat.train$Mort.prop, predict(brt_orig, dat.train, n.trees=best_iter_orig))^2),
    sprintf("%.3f", cor(dat.train$Mort.prop, predict(brt_reg2, dat.train, n.trees=best_iter_reg2))^2),
    sprintf("%.3f", cor(dat.train$Mort.prop, predict(brt_reg1, dat.train, n.trees=best_iter_reg1))^2)
  ),
  Spatial_Blocked_CV_Deviance = c(sp_orig["Dev"], sp_reg2["Dev"], sp_reg1["Dev"]),
  Spatial_Blocked_CV_R2 = c(sp_orig["R2"], sp_reg2["R2"], sp_reg1["R2"]),
  Spatial_Blocked_CV_RMSE = c(sp_orig["RMSE"], sp_reg2["RMSE"], sp_reg1["RMSE"])
)

cat("\n=======================================================================\n")
cat("      BRT OVERFITTING MITIGATION & SPATIAL BLOCKED CV COMPARISON\n")
cat("=======================================================================\n")
print(eval_summary)

saveRDS(brt_reg2, "output/models/brt_manta_reg_fit.rds")
saveRDS(eval_summary, "output/models/brt_regularization_comparison.rds")
