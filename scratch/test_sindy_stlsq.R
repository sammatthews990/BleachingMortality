library(dplyr)
library(glmnet)

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

# Add key ecological interactions and non-linear threshold terms
lib[["DHW^2"]] <- raw_df$DHW^2
lib[["Secchi^2"]] <- raw_df$Secchi^2
lib[["DHW * Secchi"]] <- raw_df$DHW * raw_df$Secchi
lib[["DHW * WinMean"]] <- raw_df$DHW * raw_df$WinMean
lib[["DHW * HistDHW"]] <- raw_df$DHW * raw_df$HistDHW
lib[["DHW * Cloud"]] <- raw_df$DHW * raw_df$Cloud
lib[["Secchi * WinMean"]] <- raw_df$Secchi * raw_df$WinMean

# Key non-linear ecological threshold functions
lib[["max(0, DHW - 4)"]] <- pmax(0, raw_df$DHW - 4)
lib[["max(0, DHW - 6)"]] <- pmax(0, raw_df$DHW - 6)
lib[["DHW * max(0, Secchi - 25)"]] <- raw_df$DHW * pmax(0, raw_df$Secchi - 25)

Theta_raw <- as.data.frame(lib)

Theta_mean <- colMeans(Theta_raw)
Theta_sd <- apply(Theta_raw, 2, sd)
Theta_sd[Theta_sd == 0] <- 1

Theta_scaled <- scale(Theta_raw, center = Theta_mean, scale = Theta_sd)

stlsq <- function(X, y, cutoff = 0.1, max_iter = 50) {
  # OLS fit on standardized features
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

cutoffs <- c(0.05, 0.08, 0.10, 0.12, 0.15, 0.18, 0.20, 0.25)

for (cut in cutoffs) {
  xi_scaled <- stlsq(as.matrix(Theta_scaled), Y, cutoff = cut)
  active <- xi_scaled[xi_scaled != 0]
  
  if (length(active) > 0) {
    raw_coefs <- active / Theta_sd[names(active)]
    raw_intercept <- mean(Y) - sum(active * Theta_mean[names(active)] / Theta_sd[names(active)])
    preds <- raw_intercept + as.matrix(Theta_raw[, names(active), drop = FALSE]) %*% raw_coefs
    
    r2 <- cor(Y, preds)^2
    rmse <- sqrt(mean((Y - preds)^2))
    
    cat(sprintf("\n=======================================================\n"))
    cat(sprintf("SINDy STLSQ Threshold (lambda) = %.2f | Terms: %d | R2: %.3f | RMSE: %.4f\n", cut, length(active), r2, rmse))
    cat(sprintf("Equ: Mortality = %.4f", raw_intercept))
    for (nm in names(raw_coefs)) {
      val <- raw_coefs[nm]
      sign_str <- if (val >= 0) " + " else " - "
      cat(sprintf("%s%.5f * [%s]", sign_str, abs(val), nm))
    }
    cat("\n")
  }
}
