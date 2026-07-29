library(dplyr)
library(glmnet)
library(ggplot2)

source("scratch/check_disturbances_cheung.R") # produces dat.ml.mant.raw

dat.clean <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

cat("Clean dataset size N =", nrow(dat.clean), "\n")

# Target vector
Y <- dat.clean$Mort.prop

# Predictors
DHW <- dat.clean$MaxDHW.mean
Secchi <- dat.clean$secc3m
WinMean <- dat.clean$winyear_mean
HistDHW <- dat.clean$histmDHW6
Current <- dat.clean$mcur_90
WinSD <- dat.clean$winyear_sd
Cloud <- dat.clean$cloudp_90
YrSince <- dat.clean$yrsince6

# Construct SINDy Candidate Feature Library Theta(X)
# Include linear, quadratic, cubic, cross-interactions, and threshold basis functions
raw_df <- data.frame(
  DHW = DHW,
  Secchi = Secchi,
  WinMean = WinMean,
  HistDHW = HistDHW,
  Current = Current,
  WinSD = WinSD,
  Cloud = Cloud,
  YrSince = YrSince
)

# Build polynomial and interaction library up to degree 2 (and select degree 3 DHW terms)
build_sindy_library <- function(df) {
  lib <- list()
  
  # Degree 1 (Linear)
  for (nm in colnames(df)) {
    lib[[nm]] <- df[[nm]]
  }
  
  # Degree 2 (Quadratic & 2-Way Interactions)
  cols <- colnames(df)
  for (i in 1:length(cols)) {
    for (j in i:length(cols)) {
      nm1 <- cols[i]
      nm2 <- cols[j]
      if (i == j) {
        lib[[paste0(nm1, "^2")]] <- df[[nm1]]^2
      } else {
        lib[[paste0(nm1, ":", nm2)]] <- df[[nm1]] * df[[nm2]]
      }
    }
  }
  
  # Non-linear ecological threshold functions (e.g. max(0, DHW - 4), DHW^3)
  lib[["DHW^3"]] <- df$DHW^3
  lib[["max(0, DHW - 4)"]] <- pmax(0, df$DHW - 4)
  lib[["max(0, DHW - 6)"]] <- pmax(0, df$DHW - 6)
  lib[["DHW * max(0, Secchi - 25)"]] <- df$DHW * pmax(0, df$Secchi - 25)
  lib[["DHW / Secchi"]] <- df$DHW / pmax(df$Secchi, 1)
  
  as.data.frame(lib)
}

Theta_raw <- build_sindy_library(raw_df)
cat("Candidate Library Theta shape:", nrow(Theta_raw), "x", ncol(Theta_raw), "\n")

# Standardize Library for STLSQ (Sequential Thresholded Least Squares)
Theta_mean <- colMeans(Theta_raw)
Theta_sd <- apply(Theta_raw, 2, sd)
# avoid division by zero if sd == 0
Theta_sd[Theta_sd == 0] <- 1

Theta_scaled <- scale(Theta_raw, center = Theta_mean, scale = Theta_sd)

# Implement STLSQ (Sequential Thresholded Least Squares - SINDy core algorithm)
stlsq <- function(X, y, cutoff = 0.05, max_iter = 25) {
  # Initial OLS fit
  fit <- lm.fit(X, y)
  coefs <- fit$coefficients
  coefs[is.na(coefs)] <- 0
  
  for (k in 1:max_iter) {
    small_inds <- abs(coefs) < cutoff
    coefs[small_inds] <- 0
    big_inds <- !small_inds
    
    if (sum(big_inds) == 0) break
    
    # Re-fit OLS on active features
    fit <- lm.fit(X[, big_inds, drop = FALSE], y)
    coefs_new <- numeric(ncol(X))
    names(coefs_new) <- colnames(X)
    coefs_new[big_inds] <- fit$coefficients
    
    if (max(abs(coefs_new - coefs)) < 1e-6) break
    coefs <- coefs_new
  }
  
  coefs
}

# Test a grid of cutoff thresholds for STLSQ
cutoffs <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.08, 0.10)
results_list <- list()

for (cut in cutoffs) {
  xi_scaled <- stlsq(Theta_scaled, Y, cutoff = cut)
  
  # Unscale coefficients back to raw feature units
  # Y_hat = Intercept + sum_j (xi_scaled_j * (X_j - mean_j) / sd_j)
  # Intercept_raw = mean(Y) - sum_j (xi_scaled_j * mean_j / sd_j)
  # Coef_raw_j = xi_scaled_j / sd_j
  active <- xi_scaled[xi_scaled != 0]
  
  raw_coefs <- active / Theta_sd[names(active)]
  raw_intercept <- mean(Y) - sum(active * Theta_mean[names(active)] / Theta_sd[names(active)])
  
  # Predictions
  preds <- raw_intercept + as.matrix(Theta_raw[, names(active), drop = FALSE]) %*% raw_coefs
  
  r2 <- cor(Y, preds)^2
  rmse <- sqrt(mean((Y - preds)^2))
  n_terms <- length(active)
  
  cat(sprintf("\nCutoff: %.2f | Active Terms: %d | R2: %.3f | RMSE: %.4f\n", cut, n_terms, r2, rmse))
  cat("Equation:\n  Mortality = ", sprintf("%.4f", raw_intercept))
  for (nm in names(raw_coefs)) {
    val <- raw_coefs[nm]
    sign_str <- if (val >= 0) " + " else " - "
    cat(sprintf("%s%.5f * [%s]", sign_str, abs(val), nm))
  }
  cat("\n")
}

# Also run Sparse Lasso / ElasticNet for comparison
cv_lasso <- cv.glmnet(as.matrix(Theta_scaled), Y, alpha = 1, nfolds = 5)
best_lambda <- cv_lasso$lambda.1se
lasso_coefs <- coef(cv_lasso, s = best_lambda)
cat("\n=== Lasso (SINDy variant via L1 Regularization) ===\n")
print(lasso_coefs[lasso_coefs[, 1] != 0, , drop = FALSE])
