library(dplyr)
library(gbm)
library(ggplot2)
library(patchwork)
library(purrr)

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

predictors <- c("MaxDHW.mean", "secc3m", "cloudp_90", "histmDHW6", "mcur_90", "winyear_sd")
pred_labels <- c(
  "MaxDHW.mean" = "Max DHW (°C-weeks)",
  "secc3m" = "Secchi Depth (m)",
  "cloudp_90" = "Cloud Cover (90d)",
  "histmDHW6" = "Prior Exposure (>6 DHW)",
  "mcur_90" = "Current Speed (90d)",
  "winyear_sd" = "Winter SST SD"
)

# 1. Fit main BRT model
set.seed(42)
brt_fit <- gbm(
  formula = Mort.prop ~ MaxDHW.mean + secc3m + cloudp_90 + histmDHW6 + mcur_90 + winyear_sd,
  data = dat.ml.mant.clean,
  distribution = "gaussian",
  n.trees = 2000,
  interaction.depth = 3,
  shrinkage = 0.005,
  bag.fraction = 0.75,
  cv.folds = 5
)

best_iter <- gbm.perf(brt_fit, plot.it = FALSE, method = "cv")
cat("Optimal boosting trees (5-fold CV):", best_iter, "\n")

# Variable importance
var_imp <- summary(brt_fit, n.trees = best_iter, plotit = FALSE)
var_imp$Variable_Label <- pred_labels[as.character(var_imp$var)]

# 2. Bootstrap Partial Dependence Plots (B = 100)
B <- 100
cat("Running", B, "bootstrap iterations for BRT Partial Dependence Plots...\n")

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
    formula = Mort.prop ~ MaxDHW.mean + secc3m + cloudp_90 + histmDHW6 + mcur_90 + winyear_sd,
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
    lo_90 = quantile(yhat, 0.05),
    hi_90 = quantile(yhat, 0.95),
    .groups = "drop"
  ) %>%
  mutate(pred_label = pred_labels[predictor])

if (!dir.exists("output/models")) dir.create("output/models", recursive = TRUE)

saveRDS(brt_fit, "output/models/brt_manta_fit.rds")
saveRDS(var_imp, "output/models/brt_manta_var_imp.rds")
saveRDS(pdp_summary, "output/models/brt_manta_pdp_boot.rds")

cat("=== BRT fitting and PDP bootstrapping complete! ===\n")
