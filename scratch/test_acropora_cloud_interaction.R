# scratch/test_acropora_cloud_interaction.R
# Sourcing the main data pipeline
cat("Sourcing data pipeline...\n")
source("scratch/DHW_Mortality_Attribution.R")
cat("Finished sourcing data pipeline.\n\n")

library(dplyr)
library(ggplot2)
library(brms)

# ----------------- LTMP Benthic Comparison -----------------
cat("=========================================\n")
cat("       LTMP BENTHIC MODEL COMPARISON     \n")
cat("=========================================\n")

# Frequentist Original
fit_qbin_ltmp_orig <- fit_qbin_ltmp
# Frequentist Alternative: % acropora * cloudp_90 instead of MaxDHW.mean * prop_acropora
fit_qbin_ltmp_alt <- glm(
  Mort.prop.nudge ~ MaxDHW.mean * secc3m + MaxDHW.mean * cloudp_90 +
                    MaxDHW.mean * histmDHW6 + MaxDHW.mean * mcur_90 +
                    prop_acropora * cloudp_90 + MaxDHW.mean + winyear_sd,
  data = dat.ml.ltmp.clean,
  family = quasibinomial(link = "logit")
)

cat("=== Frequentist LTMP Original (Acropora * DHW) ===\n")
print(summary(fit_qbin_ltmp_orig)$coefficients["MaxDHW.mean:prop_acropora", , drop = FALSE])
cat("\n=== Frequentist LTMP Alternative (Acropora * Cloud Cover) ===\n")
print(summary(fit_qbin_ltmp_alt)$coefficients["cloudp_90:prop_acropora", , drop = FALSE])

# Let's compare R^2 (quasibinomial pseudo-R^2 based on deviance)
r2_pseudo <- function(mod) {
  1 - (mod$deviance / mod$null.deviance)
}
cat("\nLTMP Original deviance R²:", round(r2_pseudo(fit_qbin_ltmp_orig), 4), "\n")
cat("LTMP Alternative deviance R²:", round(r2_pseudo(fit_qbin_ltmp_alt), 4), "\n")


# ----------------- MMP Inshore Comparison -----------------
cat("\n=========================================\n")
cat("        MMP INSHORE MODEL COMPARISON     \n")
cat("=========================================\n")

# Frequentist Original
fit_qbin_mmp_orig <- fit_qbin_mmp
# Frequentist Alternative: % acropora * cloudp_90 instead of MaxDHW.mean * prop_acropora
fit_qbin_mmp_alt <- glm(
  Mort.prop.nudge ~ MaxDHW.mean * secc3m + MaxDHW.mean * cloudp_90 +
                    MaxDHW.mean * histmDHW6 + MaxDHW.mean * winyear_sd +
                    MaxDHW.mean * mcur_90 + prop_acropora * cloudp_90 +
                    MaxDHW.mean + depth,
  data = dat.ml.mmp.clean,
  family = quasibinomial(link = "logit")
)

cat("=== Frequentist MMP Original (Acropora * DHW) ===\n")
print(summary(fit_qbin_mmp_orig)$coefficients["MaxDHW.mean:prop_acropora", , drop = FALSE])
cat("\n=== Frequentist MMP Alternative (Acropora * Cloud Cover) ===\n")
print(summary(fit_qbin_mmp_alt)$coefficients["cloudp_90:prop_acropora", , drop = FALSE])

cat("\nMMP Original deviance R²:", round(r2_pseudo(fit_qbin_mmp_orig), 4), "\n")
cat("MMP Alternative deviance R²:", round(r2_pseudo(fit_qbin_mmp_alt), 4), "\n")


# ----------------- Bayesian brms models comparison -----------------
cat("\n=========================================\n")
cat("        BAYESIAN MODEL COMPARISON        \n")
cat("=========================================\n")

# We can fit alternative Bayesian Beta models in brms and run LOO-CV comparison!
cache_ltmp_beta_alt <- "output/models/brms_ltmp_beta_alt.rds"
if (file.exists(cache_ltmp_beta_alt)) {
  fit_brms_ltmp_alt <- readRDS(cache_ltmp_beta_alt)
  cat("Loaded cached LTMP Beta Alternative model.\n")
} else {
  cat("Fitting LTMP Beta Alternative model (Acropora * Cloud Cover)...\n")
  fit_brms_ltmp_alt <- brm(
    Mort.prop.nudge ~ MaxDHW_s * cloudp_90_s + MaxDHW_s * histmDHW6_s +
                      MaxDHW_s * mcur_90_s + prop_acropora_s * cloudp_90_s +
                      secc3m_s + winyear_sd_s + MaxDHW_s + (1 | ReefName),
    data = dat.ml.ltmp.clean, family = Beta(link = "logit"),
    prior = c(prior(normal(-1, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.97)
  )
  saveRDS(fit_brms_ltmp_alt, cache_ltmp_beta_alt)
}

cache_mmp_beta_alt <- "output/models/brms_mmp_beta_alt.rds"
if (file.exists(cache_mmp_beta_alt)) {
  fit_brms_mmp_alt <- readRDS(cache_mmp_beta_alt)
  cat("Loaded cached MMP Beta Alternative model.\n")
} else {
  cat("Fitting MMP Beta Alternative model (Acropora * Cloud Cover)...\n")
  fit_brms_mmp_alt <- brm(
    Mort.prop.nudge ~ MaxDHW_s * cloudp_90_s + prop_acropora_s * cloudp_90_s +
                      secc3m_s + depth + MaxDHW_s + (1 | ReefName),
    data = dat.ml.mmp.clean, family = Beta(link = "logit"),
    prior = c(prior(normal(-1, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(2), class = "sd")),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.98)
  )
  saveRDS(fit_brms_mmp_alt, cache_mmp_beta_alt)
}

# Compare Bayesian Models using LOO-CV and Bayes R^2
cat("\n=== LTMP Bayesian Comparison ===\n")
cat("Original Beta R²:", round(median(bayes_R2(fit_brms_ltmp)), 4), "\n")
cat("Alternative Beta R²:", round(median(bayes_R2(fit_brms_ltmp_alt)), 4), "\n")
cat("\nLOO Comparison (LTMP):\n")
loo_ltmp_orig <- loo(fit_brms_ltmp)
loo_ltmp_alt <- loo(fit_brms_ltmp_alt)
print(loo_compare(loo_ltmp_orig, loo_ltmp_alt))

cat("\n=== MMP Bayesian Comparison ===\n")
cat("Original Beta R²:", round(median(bayes_R2(fit_brms_mmp)), 4), "\n")
cat("Alternative Beta R²:", round(median(bayes_R2(fit_brms_mmp_alt)), 4), "\n")
cat("\nLOO Comparison (MMP):\n")
loo_mmp_orig <- loo(fit_brms_mmp)
loo_mmp_alt <- loo(fit_brms_mmp_alt)
print(loo_compare(loo_mmp_orig, loo_mmp_alt))

cat("\n=== Bayesian Parameter Estimates for the Interaction term ===\n")
cat("LTMP cloudp_90_s:prop_acropora_s posterior:\n")
print(summary(fit_brms_ltmp_alt)$fixed["cloudp_90_s:prop_acropora_s", , drop = FALSE])
cat("MMP cloudp_90_s:prop_acropora_s posterior:\n")
print(summary(fit_brms_mmp_alt)$fixed["cloudp_90_s:prop_acropora_s", , drop = FALSE])
