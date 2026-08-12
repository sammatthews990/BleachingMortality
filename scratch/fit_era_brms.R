suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(brms)
  library(patchwork)
})

# Load the cached environments
e_m <- new.env()
lazyLoad('DHW_Mortality_Attribution_cache/html/manta-model-data_d3d5fd24b8f548ff1210864d2d490674', envir = e_m)
dat.mod.mant <- e_m[['dat.mod.mant']] %>% st_drop_geometry()

e_b <- new.env()
lazyLoad('DHW_Mortality_Attribution_cache/html/benthic-data_3436ddf577b9428741a2f0fe11c82b03', envir = e_b)
dat.mod.bent.hc <- e_b[['dat.mod.bent.hc']] %>% st_drop_geometry()

zarr_full_raw <- read.csv('data/processed/cheung_recreated_gbr_full.csv')

d_m_all <- dat.mod.mant %>%
    ungroup() %>%
    filter(!is.na(ReefID), !is.na(DHWYear)) %>%
    mutate(DHWYear = as.numeric(DHWYear),
           ReefID_clean = toupper(trimws(as.character(ReefID)))) %>%
    left_join(
        zarr_full_raw %>% mutate(LABEL_ID_clean = toupper(trimws(as.character(LABEL_ID)))),
        by = c('ReefID_clean' = 'LABEL_ID_clean', 'DHWYear' = 'year')
    ) %>%
    mutate(Dataset = 'Manta Tow', survey_depth = 7.5)

d_l_all <- dat.mod.bent.hc %>%
    ungroup() %>%
    filter(project_code == 'LTMP', !is.na(ReefID), !is.na(DHWYear)) %>%
    mutate(DHWYear = as.numeric(DHWYear),
           ReefID_clean = toupper(trimws(as.character(ReefID)))) %>%
    left_join(
        zarr_full_raw %>% mutate(LABEL_ID_clean = toupper(trimws(as.character(LABEL_ID)))),
        by = c('ReefID_clean' = 'LABEL_ID_clean', 'DHWYear' = 'year')
    ) %>%
    mutate(Dataset = 'LTMP Benthic')

d_p_all <- dat.mod.bent.hc %>%
    ungroup() %>%
    filter(project_code == 'MMP', !is.na(ReefID), !is.na(DHWYear)) %>%
    mutate(DHWYear = as.numeric(DHWYear),
           ReefID_clean = toupper(trimws(as.character(ReefID)))) %>%
    left_join(
        zarr_full_raw %>% mutate(LABEL_ID_clean = toupper(trimws(as.character(LABEL_ID)))),
        by = c('ReefID_clean' = 'LABEL_ID_clean', 'DHWYear' = 'year')
    ) %>%
    mutate(Dataset = 'MMP Inshore')

dat.ml.expanded <- bind_rows(d_m_all, d_l_all, d_p_all) %>%
    mutate(MaxDHW_val = coalesce(ann_maxdhw, MaxDHW.mean)) %>%
    filter(!is.na(Rel.Change), !is.na(MaxDHW_val))

# Define Era: 3 eras
dat.ml.expanded <- dat.ml.expanded %>%
    mutate(
        Era = case_when(
            DHWYear <= 2002 ~ "Historical (1998-2002)",
            DHWYear >= 2016 & DHWYear <= 2020 ~ "Modern Baseline (2016-2020)",
            DHWYear >= 2022 ~ "Contemporary (2022-2024)",
            TRUE ~ "Intermediate"
        ),
        Era = factor(Era, levels = c("Historical (1998-2002)", "Modern Baseline (2016-2020)", "Contemporary (2022-2024)"))
    ) %>%
    filter(!is.na(Era))

# Also define discrete Major Event Year (1998, 2002, 2016, 2017, 2020, 2022, 2024)
dat.ml.expanded <- dat.ml.expanded %>%
    mutate(
        EventYear = as.character(DHWYear),
        EventYear = factor(EventYear, levels = c("1998", "2002", "2016", "2017", "2020", "2022", "2024"))
    )

# Prepare scaled variables
N_all <- nrow(dat.ml.expanded)
dat.ml.expanded <- as.data.frame(dat.ml.expanded)
for (v in c("secc3m", "winyear_mean", "histmDHW6", "mcur_90", "winyear_sd", "cloudp_90", "yrsince6", "prop_acropora")) {
    if (v %in% names(dat.ml.expanded)) {
        vals <- dat.ml.expanded[[v]]
        med_val <- median(vals, na.rm = TRUE)
        if (is.na(med_val)) med_val <- 0
        vals[is.na(vals)] <- med_val
        dat.ml.expanded[[v]] <- vals
    }
}

mean_dhw <- mean(dat.ml.expanded$MaxDHW_val, na.rm=TRUE)
sd2_dhw <- 2 * sd(dat.ml.expanded$MaxDHW_val, na.rm=TRUE)

dat.ml.expanded <- dat.ml.expanded %>%
    mutate(
        Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
        Mort.prop.nudge = (Mort.prop * (N_all - 1) + 0.5) / N_all,
        MaxDHW_s = (MaxDHW_val - mean_dhw) / sd2_dhw,
        Year_s = (DHWYear - 2016) / 10,
        secc3m_s = (secc3m - 17.7) / 11.0,
        cloudp_90_s = (cloudp_90 - 0.72) / 0.2,
        prop_acropora_s = (prop_acropora - 0.35) / 0.3,
        winyear_sd_s = (winyear_sd - mean(winyear_sd)) / (2 * sd(winyear_sd)),
        histmDHW6_s = (histmDHW6 - mean(histmDHW6)) / (2 * sd(histmDHW6)),
        mcur_90_s = (mcur_90 - mean(mcur_90)) / (2 * sd(mcur_90)),
        winyear_mean_s = (winyear_mean - mean(winyear_mean)) / (2 * sd(winyear_mean)),
        yrsince6_s = (yrsince6 - mean(yrsince6)) / (2 * sd(yrsince6)),
        ReefName_clean = coalesce(as.character(Reef_Name), as.character(ReefName), as.character(ReefID))
    )

# Fit Era Bayesian Model
cache_era_brms <- "output/models/brms_era_beta.rds"
if (file.exists(cache_era_brms)) {
    fit_brms_era <- readRDS(cache_era_brms)
} else {
    fit_brms_era <- brm(
        Mort.prop.nudge ~ MaxDHW_s * Era + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          winyear_sd_s + histmDHW6_s + mcur_90_s + Dataset + (1 | ReefName_clean),
        data = dat.ml.expanded, family = Beta(link = "logit"),
        prior = c(
            prior(normal(-1, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")
        ),
        chains = 4, iter = 2500, warmup = 1000, cores = 4, seed = 42,
        control = list(adapt_delta = 0.95)
    )
    saveRDS(fit_brms_era, cache_era_brms)
}

cat("BRMS Era model loaded!\n")
print(summary(fit_brms_era))
