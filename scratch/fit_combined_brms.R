# Fit Bayesian models for combined 3-dataset model suite
suppressPackageStartupMessages({
    library(dplyr)
    library(brms)
    library(readr)
})

# Load environment / run minimal qmd code to obtain individual clean datasets
lines <- readLines("scratch/prep.R")

# Execute up to qbin-mmp-model
qbin_mmp_line <- grep("qbin-mmp-model", lines)[1]
eval(parse(text = lines[1:(qbin_mmp_line + 50)]))

cat("Manta Tow rows:", nrow(dat.ml.mant.clean), "\n")
cat("LTMP Benthic rows:", nrow(dat.ml.ltmp.clean), "\n")
cat("MMP Inshore rows:", nrow(dat.ml.mmp.clean), "\n")

# Combine datasets into unified dataframe
cols_needed <- c(
    "ReefName", "Sector", "MaxDHW.mean", "secc3m", "winyear_mean",
    "histmDHW6", "mcur_90", "winyear_sd", "cloudp_90", "yrsince6",
    "prop_acropora", "Mort.prop", "Rel.Change"
)

d_mant <- dat.ml.mant.clean %>% select(any_of(cols_needed)) %>% mutate(Dataset = "Manta Tow")
d_ltmp <- dat.ml.ltmp.clean %>% select(any_of(cols_needed)) %>% mutate(Dataset = "LTMP Benthic")
d_mmp  <- dat.ml.mmp.clean  %>% select(any_of(cols_needed)) %>% mutate(Dataset = "MMP Inshore")

dat.ml.comb.clean <- bind_rows(d_mant, d_ltmp, d_mmp)

# Handle missing Sector if any
if (!"Sector" %in% names(dat.ml.comb.clean) || any(is.na(dat.ml.comb.clean$Sector))) {
    dat.ml.comb.clean$Sector <- ifelse(is.na(dat.ml.comb.clean$Sector), dat.ml.comb.clean$ReefName, dat.ml.comb.clean$Sector)
}

# Impute median for any missing predictors across combined set
pred_8 <- c("MaxDHW.mean", "secc3m", "winyear_mean", "histmDHW6", "mcur_90", "winyear_sd", "cloudp_90", "yrsince6", "prop_acropora")
for (v in pred_8) {
    if (v %in% names(dat.ml.comb.clean)) {
        dat.ml.comb.clean[[v]][is.na(dat.ml.comb.clean[[v]])] <- median(dat.ml.comb.clean[[v]], na.rm = TRUE)
    }
}

N_comb <- nrow(dat.ml.comb.clean)
dat.ml.comb.clean <- dat.ml.comb.clean %>%
    mutate(
        Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
        Mort.prop.nudge = (Mort.prop * (N_comb - 1) + 0.5) / N_comb,
        Dataset = factor(Dataset, levels = c("Manta Tow", "LTMP Benthic", "MMP Inshore")),
        obs_id = 1:n(),
        n_trials = 100L,
        y_binom = as.integer(round(Mort.prop.nudge * 100))
    )

s2 <- function(x) (x - mean(x, na.rm = TRUE)) / (2 * sd(x, na.rm = TRUE))
dat.ml.comb.clean <- dat.ml.comb.clean %>%
    mutate(
        MaxDHW_s = s2(MaxDHW.mean),
        secc3m_s = s2(secc3m),
        winyear_mean_s = s2(winyear_mean),
        histmDHW6_s = s2(histmDHW6),
        mcur_90_s = s2(mcur_90),
        winyear_sd_s = s2(winyear_sd),
        cloudp_90_s = s2(cloudp_90),
        yrsince6_s = s2(yrsince6),
        prop_acropora_s = s2(prop_acropora)
    )

cat("Combined dataset successfully created with N =", nrow(dat.ml.comb.clean), "observations across", length(unique(dat.ml.comb.clean$ReefName)), "reefs.\n")

# 1. Fit Bayesian Beta GLMM
cache_comb_beta <- "output/models/brms_combined_beta.rds"
if (!file.exists(cache_comb_beta)) {
    cat("Fitting Combined Bayesian Beta GLMM...\n")
    fit_brms_comb <- brm(
        Mort.prop.nudge ~ MaxDHW_s * secc3m_s + MaxDHW_s * cloudp_90_s + MaxDHW_s * prop_acropora_s +
            winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + Dataset +
            (1 | Dataset/ReefName),
        data = dat.ml.comb.clean, family = Beta(link = "logit"),
        prior = c(
            prior(normal(-1, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(2), class = "sd")
        ),
        chains = 4, iter = 3000, warmup = 1500, cores = 4, seed = 42,
        control = list(adapt_delta = 0.96)
    )
    saveRDS(fit_brms_comb, cache_comb_beta)
    cat("Saved Combined Beta model to", cache_comb_beta, "\n")
} else {
    cat("Combined Beta model already exists.\n")
}

# 2. Fit Bayesian Binomial-OLRE Model
cache_comb_binom <- "output/models/brms_combined_binomial_olre.rds"
if (!file.exists(cache_comb_binom)) {
    cat("Fitting Combined Bayesian Binomial-OLRE Model...\n")
    fit_brmsb_comb <- brm(
        y_binom | trials(n_trials) ~ MaxDHW_s * secc3m_s + MaxDHW_s * cloudp_90_s + MaxDHW_s * prop_acropora_s +
            winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + Dataset +
            (1 | Dataset/ReefName) + (1 | obs_id),
        data = dat.ml.comb.clean, family = binomial(link = "logit"),
        prior = c(
            prior(normal(-1, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(2), class = "sd")
        ),
        chains = 4, iter = 3000, warmup = 1500, cores = 4, seed = 42,
        control = list(adapt_delta = 0.96)
    )
    saveRDS(fit_brmsb_comb, cache_comb_binom)
    cat("Saved Combined Binomial-OLRE model to", cache_comb_binom, "\n")
} else {
    cat("Combined Binomial-OLRE model already exists.\n")
}

cat("All Bayesian combined models successfully prepared.\n")
