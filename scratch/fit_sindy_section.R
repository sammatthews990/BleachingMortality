library(dplyr)
library(ggplot2)
library(patchwork)

source("scratch/check_disturbances_cheung.R") # produces dat.ml.mant.raw

dat.clean <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

Y <- dat.clean$Mort.prop

raw_df <- data.frame(
  DHW = dat.clean$MaxDHW.mean,
  Secchi = dat.clean$secc3m,
  WinMean = dat.clean$winyear_mean,
  HistDHW = dat.clean$histmDHW6,
  Current = dat.clean$mcur_90,
  WinSD = dat.clean$winyear_sd,
  Cloud = dat.clean$cloudp_90,
  YrSince = dat.clean$yrsince6
)

# Build SINDy Candidate Feature Library Theta(X)
lib <- list()
for (nm in colnames(raw_df)) lib[[nm]] <- raw_df[[nm]]

# Key interactions & non-linear features
lib[["DHW_sq"]] <- raw_df$DHW^2
lib[["Secchi_sq"]] <- raw_df$Secchi^2
lib[["DHW_Secchi"]] <- raw_df$DHW * raw_df$Secchi
lib[["DHW_WinMean"]] <- raw_df$DHW * raw_df$WinMean
lib[["DHW_HistDHW"]] <- raw_df$DHW * raw_df$HistDHW
lib[["DHW_Cloud"]] <- raw_df$DHW * raw_df$Cloud
lib[["Secchi_WinMean"]] <- raw_df$Secchi * raw_df$WinMean
lib[["max_DHW_4"]] <- pmax(0, raw_df$DHW - 4)
lib[["max_DHW_6"]] <- pmax(0, raw_df$DHW - 6)

Theta_raw <- as.data.frame(lib)
Theta_mean <- colMeans(Theta_raw)
Theta_sd <- apply(Theta_raw, 2, sd)
Theta_sd[Theta_sd == 0] <- 1

Theta_scaled <- as.matrix(scale(Theta_raw, center = Theta_mean, scale = Theta_sd))

stlsq <- function(X, y, cutoff = 0.1, max_iter = 50) {
  fit <- lm.fit(X, y)
  coefs <- fit$coefficients
  coefs[is.na(coefs)] <- 0
  
  for (k in 1:max_iter) {
    small_inds <- abs(coefs) < cutoff
    coefs[small_inds] <- 0
    big_inds <- !small_inds
    
    if (sum(big_inds) == 0) break
    
    fit <- lm.fit(X[, big_inds, drop = FALSE], y)
    coefs_new <- numeric(ncol(X))
    names(coefs_new) <- colnames(X)
    coefs_new[big_inds] <- fit$coefficients
    
    if (max(abs(coefs_new - coefs)) < 1e-6) break
    coefs <- coefs_new
  }
  
  coefs
}

# Cross-Validation Evaluation of SINDy
set.seed(42)
folds <- sample(rep(1:5, length.out = nrow(dat.clean)))
cv_results <- list()

cutoffs <- c(0.05, 0.08, 0.10, 0.15, 0.18, 0.20)

for (cut in cutoffs) {
  cv_preds <- numeric(nrow(dat.clean))
  
  for (f in 1:5) {
    train_idx <- folds != f
    val_idx <- folds == f
    
    X_tr <- Theta_scaled[train_idx, ]
    y_tr <- Y[train_idx]
    
    xi_tr <- stlsq(X_tr, y_tr, cutoff = cut)
    active <- xi_tr[xi_tr != 0]
    
    if (length(active) > 0) {
      r_coefs <- active / Theta_sd[names(active)]
      r_ic <- mean(y_tr) - sum(active * Theta_mean[names(active)] / Theta_sd[names(active)])
      
      X_val_raw <- as.matrix(Theta_raw[val_idx, names(active), drop = FALSE])
      cv_preds[val_idx] <- r_ic + X_val_raw %*% r_coefs
    } else {
      cv_preds[val_idx] <- mean(y_tr)
    }
  }
  
  cv_r2 <- cor(Y, cv_preds)^2
  cv_rmse <- sqrt(mean((Y - cv_preds)^2))
  
  xi_full <- stlsq(Theta_scaled, Y, cutoff = cut)
  act_full <- xi_full[xi_full != 0]
  
  cv_results[[length(cv_results) + 1]] <- data.frame(
    Cutoff = cut,
    Active_Terms = length(act_full),
    CV_R2 = cv_r2,
    CV_RMSE = cv_rmse
  )
}

cv_summary_df <- bind_rows(cv_results)

xi_opt <- stlsq(Theta_scaled, Y, cutoff = 0.10)
act_opt <- xi_opt[xi_opt != 0]

r_coefs_opt <- act_opt / Theta_sd[names(act_opt)]
r_ic_opt <- mean(Y) - sum(act_opt * Theta_mean[names(act_opt)] / Theta_sd[names(act_opt)])

feature_labels <- c(
  "Intercept" = "Intercept",
  "DHW" = "DHW (°C-weeks)",
  "Secchi" = "Secchi Depth (m)",
  "WinMean" = "Winter SST Mean (°C)",
  "HistDHW" = "Prior Exposure (>6 DHW)",
  "DHW_sq" = "DHW² (Thermal Acceleration)",
  "Secchi_sq" = "Secchi² (Clarity Non-linearity)",
  "DHW_Secchi" = "DHW × Secchi Depth",
  "DHW_WinMean" = "DHW × Winter SST Mean",
  "DHW_HistDHW" = "DHW × Prior Exposure",
  "Secchi_WinMean" = "Secchi × Winter SST Mean",
  "max_DHW_4" = "max(0, DHW - 4)",
  "max_DHW_6" = "max(0, DHW - 6)"
)

eq_parts <- c(sprintf("Mortality = %.4f", r_ic_opt))
for (nm in names(r_coefs_opt)) {
  val <- r_coefs_opt[nm]
  sign_str <- if (val >= 0) " + " else " - "
  eq_parts <- c(eq_parts, sprintf("%s%.5f·[%s]", sign_str, abs(val), nm))
}
sindy_eq_str <- paste(eq_parts, collapse = " ")

sindy_coefs_df <- data.frame(
  Term = c("Intercept", names(r_coefs_opt)),
  Feature_Name = c("Intercept", feature_labels[names(r_coefs_opt)]),
  Coefficient = c(r_ic_opt, as.numeric(r_coefs_opt)),
  Interpretation = c(
    "Baseline mortality intercept",
    "Linear thermal stress coefficient",
    "Quadratic thermal stress acceleration (non-linear threshold scaling)",
    "Non-linear water clarity depth effect",
    "Compound thermal amplification (DHW × Winter SST interaction)",
    "Light penetration × thermal interaction (Secchi × Winter SST)"
  )[1:(length(r_coefs_opt) + 1)]
)

if (!dir.exists("output/models")) dir.create("output/models", recursive = TRUE)

saveRDS(sindy_coefs_df, "output/models/sindy_manta_coefs.rds")
saveRDS(sindy_eq_str, "output/models/sindy_manta_eq_str.rds")
saveRDS(cv_summary_df, "output/models/sindy_manta_cv_summary.rds")

grid_dhw_sindy <- data.frame(
  DHW = seq(min(raw_df$DHW), max(raw_df$DHW), length.out = 100),
  Secchi = median(raw_df$Secchi),
  WinMean = median(raw_df$WinMean),
  HistDHW = median(raw_df$HistDHW),
  Current = median(raw_df$Current),
  WinSD = median(raw_df$WinSD),
  Cloud = median(raw_df$Cloud),
  YrSince = median(raw_df$YrSince)
)
grid_dhw_sindy$DHW_sq <- grid_dhw_sindy$DHW^2
grid_dhw_sindy$Secchi_sq <- grid_dhw_sindy$Secchi^2
grid_dhw_sindy$DHW_Secchi <- grid_dhw_sindy$DHW * grid_dhw_sindy$Secchi
grid_dhw_sindy$DHW_WinMean <- grid_dhw_sindy$DHW * grid_dhw_sindy$WinMean
grid_dhw_sindy$DHW_HistDHW <- grid_dhw_sindy$DHW * grid_dhw_sindy$HistDHW
grid_dhw_sindy$DHW_Cloud <- grid_dhw_sindy$DHW * grid_dhw_sindy$Cloud
grid_dhw_sindy$Secchi_WinMean <- grid_dhw_sindy$Secchi * grid_dhw_sindy$WinMean
grid_dhw_sindy$max_DHW_4 <- pmax(0, grid_dhw_sindy$DHW - 4)
grid_dhw_sindy$max_DHW_6 <- pmax(0, grid_dhw_sindy$DHW - 6)

grid_dhw_sindy$pred <- r_ic_opt + as.matrix(grid_dhw_sindy[, names(r_coefs_opt), drop = FALSE]) %*% r_coefs_opt

saveRDS(grid_dhw_sindy, "output/models/sindy_manta_dhw_grid.rds")
cat("=== SINDy pipeline completed successfully! ===\n")
