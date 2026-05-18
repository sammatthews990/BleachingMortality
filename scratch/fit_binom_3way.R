# scratch/fit_binom_3way.R
source("scratch/DHW_Mortality_Attribution.R")
library(brms)

cat("Fitting Bayesian LTMP Binomial-OLRE 3-way model...\n")
cache_ltmp_binom_3way <- "output/models/brms_ltmp_binom_3way.rds"
if (!file.exists(cache_ltmp_binom_3way)) {
  fit_brmsb_ltmp_3way <- brm(
    y_binom | trials(n_trials) ~ MaxDHW_s * cloudp_90_s * prop_acropora_s + 
                      MaxDHW_s * histmDHW6_s + MaxDHW_s * mcur_90_s + 
                      secc3m_s + winyear_sd_s + (1 | ReefName) + (1 | obs_id),
    data = dat.ml.ltmp.clean, family = binomial(link = "logit"),
    prior = c(prior(normal(-1, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.97)
  )
  saveRDS(fit_brmsb_ltmp_3way, cache_ltmp_binom_3way)
  cat("Saved LTMP Binomial-OLRE 3-way model.\n")
} else {
  cat("LTMP Binomial-OLRE 3-way model already cached.\n")
}

cat("Fitting Bayesian MMP Binomial-OLRE 3-way model...\n")
cache_mmp_binom_3way <- "output/models/brms_mmp_binom_3way.rds"
if (!file.exists(cache_mmp_binom_3way)) {
  fit_brmsb_mmp_3way <- brm(
    y_binom | trials(n_trials) ~ MaxDHW_s * cloudp_90_s * prop_acropora_s + 
                      secc3m_s + depth + (1 | ReefName) + (1 | obs_id),
    data = dat.ml.mmp.clean, family = binomial(link = "logit"),
    prior = c(prior(normal(-1, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(2), class = "sd")),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.98)
  )
  saveRDS(fit_brmsb_mmp_3way, cache_mmp_binom_3way)
  cat("Saved MMP Binomial-OLRE 3-way model.\n")
} else {
  cat("MMP Binomial-OLRE 3-way model already cached.\n")
}
