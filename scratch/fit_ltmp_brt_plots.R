library(dplyr)
library(gbm)
library(ggplot2)

load("Data/aims_ltmp/aims_ltmp.RData")
load("data/Cheungetal2025/01_sstvar_blchrf.RData")

# Rebuild dat.ml.ltmp.clean matching logic
dat.AIMSRef <- reef_info %>% distinct(Reef_Name, ReefID, ReefName, SECT_NAME, SECTOR)
dat.Ref.Clean <- dat.AIMSRef %>%
  mutate(AIMS_REEF_NAME_Clean = gsub(" REEF(S)?| ISLAND| IS", "", toupper(Reef_Name))) %>%
  group_by(AIMS_REEF_NAME_Clean) %>%
  slice(1) %>%
  ungroup()

acrop_comp <- df.AIMS.full %>%
  filter(data_type == "photo-transect", domain_category == "reef",
         purpose == "COMPOSITION", variable == "HARD CORAL") %>%
  mutate(Reef_Name_Clean = gsub(" REEF(S)?| ISLAND| IS", "", toupper(domain_name))) %>%
  left_join(dat.Ref.Clean %>% select(AIMS_REEF_NAME_Clean, ReefName),
            by = c("Reef_Name_Clean" = "AIMS_REEF_NAME_Clean")) %>%
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
            by = c("ReefName", "report_year", "depth")) %>%
  mutate(
    prop_acropora = replace_na(prop_acropora, median(prop_acropora, na.rm = TRUE)),
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  )

N_ltmp <- nrow(dat.ml.ltmp.clean)
cat("LTMP Dataset loaded:", N_ltmp, "observations\n")

predictor_vars <- c("MaxDHW.mean", "secc3m", "winyear_mean", "histmDHW6", "mcur_90", "winyear_sd", "cloudp_90", "yrsince6")
f_str <- as.formula(paste("Mort.prop ~", paste(predictor_vars, collapse = " + ")))

# 1. Fit Regularized BRT (Depth=2)
set.seed(42)
brt_fit_ltmp <- gbm(
  f_str,
  data = dat.ml.ltmp.clean,
  distribution = "gaussian",
  n.trees = 1500,
  interaction.depth = 2,
  shrinkage = 0.005,
  n.minobsinnode = 10,
  bag.fraction = 0.70,
  cv.folds = 5
)
best_iter_ltmp <- gbm.perf(brt_fit_ltmp, method = "cv", plot.it = FALSE)
cat("Best tree iteration (5-fold CV):", best_iter_ltmp, "\n")
saveRDS(brt_fit_ltmp, "output/models/brt_ltmp_fit.rds")

# 2. Variable Importance
var_imp_ltmp <- summary(brt_fit_ltmp, n.trees = best_iter_ltmp, plotit = FALSE)
var_imp_ltmp$Variable_Label <- c(
  "Max DHW (°C-weeks)", "Secchi Depth (m)", "Winter SST Mean (°C)",
  "Prior Severe Bleaching (>6 DHW)", "Current Speed (90th)", "Winter SST SD (°C)",
  "Cloud Cover (90th)", "Years Since >6 DHW"
)[match(var_imp_ltmp$var, predictor_vars)]
saveRDS(var_imp_ltmp, "output/models/brt_ltmp_var_imp.rds")

# 3. Bootstrapped Partial Dependence (B = 100)
set.seed(42)
B <- 100
grid_size <- 50
boot_pdp_list <- list()

labels_map <- c(
  MaxDHW.mean = "Max DHW (°C-weeks)",
  secc3m = "Secchi Depth (m)",
  winyear_mean = "Winter SST Mean (°C)",
  histmDHW6 = "Prior Exposure (>6 DHW)",
  mcur_90 = "Current Speed (90th)",
  winyear_sd = "Winter SST SD (°C)",
  cloudp_90 = "Cloud Cover (90th)",
  yrsince6 = "Years Since >6 DHW"
)

cat("Computing Bootstrapped PDPs...\n")
for (p_var in predictor_vars) {
  x_seq <- seq(min(dat.ml.ltmp.clean[[p_var]], na.rm = TRUE), max(dat.ml.ltmp.clean[[p_var]], na.rm = TRUE), length.out = grid_size)
  mat_yhat <- matrix(NA, nrow = B, ncol = grid_size)

  for (b in 1:B) {
    boot_idx <- sample(1:N_ltmp, replace = TRUE)
    d_boot <- dat.ml.ltmp.clean[boot_idx, ]

    fit_b <- gbm(
      f_str, data = d_boot, distribution = "gaussian",
      n.trees = best_iter_ltmp, interaction.depth = 2, shrinkage = 0.005,
      n.minobsinnode = 10, bag.fraction = 0.70, verbose = FALSE
    )

    pdp_b <- plot(fit_b, i.var = p_var, grid.levels = grid_size, return.grid = TRUE, continuous.resolution = grid_size)
    mat_yhat[b, ] <- pdp_b$yhat
  }

  df_pdp <- data.frame(
    predictor = p_var,
    pred_label = labels_map[p_var],
    x = x_seq,
    median_yhat = apply(mat_yhat, 2, median),
    lo_80 = apply(mat_yhat, 2, quantile, 0.10),
    hi_80 = apply(mat_yhat, 2, quantile, 0.90),
    lo_95 = apply(mat_yhat, 2, quantile, 0.025),
    hi_95 = apply(mat_yhat, 2, quantile, 0.975)
  )
  boot_pdp_list[[p_var]] <- df_pdp
}

pdp_boot_df_ltmp <- do.call(rbind, boot_pdp_list)
saveRDS(pdp_boot_df_ltmp, "output/models/brt_ltmp_pdp_boot.rds")

# 4. 2-Way Environmental Interactions (DHW x Secchi, DHW x WinterSST, DHW x Acropora)
cat("Computing 2-Way Interaction Grids...\n")
dhw_seq <- seq(min(dat.ml.ltmp.clean$MaxDHW.mean, na.rm = TRUE), max(dat.ml.ltmp.clean$MaxDHW.mean, na.rm = TRUE), length.out = 30)
sec_seq <- seq(min(dat.ml.ltmp.clean$secc3m, na.rm = TRUE), max(dat.ml.ltmp.clean$secc3m, na.rm = TRUE), length.out = 30)
winm_seq <- seq(min(dat.ml.ltmp.clean$winyear_mean, na.rm = TRUE), max(dat.ml.ltmp.clean$winyear_mean, na.rm = TRUE), length.out = 30)
acrop_seq <- seq(min(dat.ml.ltmp.clean$prop_acropora, na.rm = TRUE), max(dat.ml.ltmp.clean$prop_acropora, na.rm = TRUE), length.out = 30)

grid_sec <- expand.grid(MaxDHW.mean = dhw_seq, secc3m = sec_seq)
for (v in setdiff(predictor_vars, c("MaxDHW.mean", "secc3m"))) {
  grid_sec[[v]] <- median(dat.ml.ltmp.clean[[v]], na.rm = TRUE)
}
grid_sec$pred <- predict(brt_fit_ltmp, newdata = grid_sec, n.trees = best_iter_ltmp, type = "response")
saveRDS(grid_sec, "output/models/brt_ltmp_pdp_2way_sec.rds")

grid_winm <- expand.grid(MaxDHW.mean = dhw_seq, winyear_mean = winm_seq)
for (v in setdiff(predictor_vars, c("MaxDHW.mean", "winyear_mean"))) {
  grid_winm[[v]] <- median(dat.ml.ltmp.clean[[v]], na.rm = TRUE)
}
grid_winm$pred <- predict(brt_fit_ltmp, newdata = grid_winm, n.trees = best_iter_ltmp, type = "response")
saveRDS(grid_winm, "output/models/brt_ltmp_pdp_2way_winm.rds")

grid_acrop <- expand.grid(MaxDHW.mean = dhw_seq, prop_acropora = acrop_seq)
for (v in setdiff(predictor_vars, c("MaxDHW.mean", "prop_acropora"))) {
  grid_acrop[[v]] <- median(dat.ml.ltmp.clean[[v]], na.rm = TRUE)
}
grid_acrop$pred <- predict(brt_fit_ltmp, newdata = grid_acrop, n.trees = best_iter_ltmp, type = "response")
saveRDS(grid_acrop, "output/models/brt_ltmp_pdp_2way_acrop.rds")

cat("LTMP BRT plots data saved successfully!\n")
