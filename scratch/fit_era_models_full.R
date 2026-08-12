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

if (!"prop_acropora" %in% names(dat.ml.expanded)) {
    dat.ml.expanded$prop_acropora <- 0.35
}

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
        Dataset = factor(Dataset, levels = c("Manta Tow", "LTMP Benthic", "MMP Inshore")),
        ReefName_clean = coalesce(as.character(Reef_Name), as.character(ReefName), as.character(ReefID))
    )

# Fit Era Bayesian Model
cache_era_brms <- "output/models/brms_era_beta.rds"
if (file.exists(cache_era_brms)) {
    cat("Loading cached brms_era_beta.rds...\n")
    fit_brms_era <- readRDS(cache_era_brms)
} else {
    cat("Fitting brms_era_beta...\n")
    fit_brms_era <- brm(
        Mort.prop.nudge ~ MaxDHW_s * Era + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          winyear_sd_s + histmDHW6_s + mcur_90_s + Dataset + (1 | ReefName_clean),
        data = dat.ml.expanded, family = Beta(link = "logit"),
        prior = c(
            prior(normal(-1, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")
        ),
        chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 42,
        control = list(adapt_delta = 0.95)
    )
    saveRDS(fit_brms_era, cache_era_brms)
}

# Fit Event Year Bayesian Model
cache_event_brms <- "output/models/brms_eventyear_beta.rds"
dat_events <- dat.ml.expanded %>% filter(!is.na(EventYear))
if (file.exists(cache_event_brms)) {
    cat("Loading cached brms_eventyear_beta.rds...\n")
    fit_brms_event <- readRDS(cache_event_brms)
} else {
    cat("Fitting brms_eventyear_beta...\n")
    fit_brms_event <- brm(
        Mort.prop.nudge ~ MaxDHW_s * EventYear + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          Dataset + (1 | ReefName_clean),
        data = dat_events, family = Beta(link = "logit"),
        prior = c(
            prior(normal(-1, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")
        ),
        chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 42,
        control = list(adapt_delta = 0.95)
    )
    saveRDS(fit_brms_event, cache_event_brms)
}

# Fit Secular Continuous Time Trend Model
cache_secular_brms <- "output/models/brms_secular_beta.rds"
if (file.exists(cache_secular_brms)) {
    cat("Loading cached brms_secular_beta.rds...\n")
    fit_brms_secular <- readRDS(cache_secular_brms)
} else {
    cat("Fitting brms_secular_beta...\n")
    fit_brms_secular <- brm(
        Mort.prop.nudge ~ MaxDHW_s * Year_s + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          winyear_sd_s + histmDHW6_s + mcur_90_s + Dataset + (1 | ReefName_clean),
        data = dat.ml.expanded, family = Beta(link = "logit"),
        prior = c(
            prior(normal(-1, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")
        ),
        chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 42,
        control = list(adapt_delta = 0.95)
    )
    saveRDS(fit_brms_secular, cache_secular_brms)
}

cat("All models fitted/loaded successfully!\n")

# Now generate a 4-panel temporal comparison figure:
# Panel 1: Dose-Response across the 3 Eras
dhw_seq <- seq(0, 15, length.out = 150)
grid_era <- expand.grid(
    MaxDHW_val = dhw_seq,
    Era = levels(dat.ml.expanded$Era)
) %>%
    mutate(
        MaxDHW_s = (MaxDHW_val - mean_dhw) / sd2_dhw,
        secc3m_s = 0, cloudp_90_s = 0, prop_acropora_s = 0,
        winyear_sd_s = 0, histmDHW6_s = 0, mcur_90_s = 0,
        Dataset = factor("Manta Tow", levels = levels(dat.ml.expanded$Dataset))
    )

pp_era <- posterior_epred(fit_brms_era, newdata = grid_era, re_formula = NA)
grid_era$pred <- apply(pp_era, 2, median)
grid_era$pred_lo <- apply(pp_era, 2, quantile, 0.05)
grid_era$pred_hi <- apply(pp_era, 2, quantile, 0.95)

p1_era <- ggplot() +
    geom_point(data = dat.ml.expanded, aes(x = MaxDHW_val, y = Mort.prop, colour = Era), alpha = 0.25, size = 1) +
    geom_ribbon(data = grid_era, aes(x = MaxDHW_val, ymin = pred_lo, ymax = pred_hi, fill = Era), alpha = 0.18, colour = NA) +
    geom_line(data = grid_era, aes(x = MaxDHW_val, y = pred, colour = Era), linewidth = 1.1) +
    scale_colour_manual(values = c("Historical (1998-2002)" = "#2b5c8f", "Modern Baseline (2016-2020)" = "#e66101", "Contemporary (2022-2024)" = "#5e3c99")) +
    scale_fill_manual(values = c("Historical (1998-2002)" = "#2b5c8f", "Modern Baseline (2016-2020)" = "#e66101", "Contemporary (2022-2024)" = "#5e3c99")) +
    labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", title = "(A) Multi-Decadal Dose-Response Across Eras", colour = "Era", fill = "Era") +
    coord_cartesian(ylim = c(0, 0.8), xlim = c(0, 15)) +
    theme_bw(base_size = 11) +
    theme(legend.position = c(0.35, 0.82), legend.background = element_rect(fill = alpha('white', 0.8)))

# Panel 2: Year-by-Year Dose-Response Trajectories for Major Bleaching Events
grid_event <- expand.grid(
    MaxDHW_val = dhw_seq,
    EventYear = levels(dat_events$EventYear)
) %>%
    mutate(
        MaxDHW_s = (MaxDHW_val - mean_dhw) / sd2_dhw,
        secc3m_s = 0, cloudp_90_s = 0, prop_acropora_s = 0,
        Dataset = factor("Manta Tow", levels = levels(dat.ml.expanded$Dataset))
    )

pp_event <- posterior_epred(fit_brms_event, newdata = grid_event, re_formula = NA)
grid_event$pred <- apply(pp_event, 2, median)
grid_event$pred_lo <- apply(pp_event, 2, quantile, 0.05)
grid_event$pred_hi <- apply(pp_event, 2, quantile, 0.95)

p2_event <- ggplot(grid_event, aes(x = MaxDHW_val, y = pred, colour = EventYear, fill = EventYear)) +
    geom_line(linewidth = 1) +
    scale_colour_viridis_d(option = "turbo", name = "Bleaching Year") +
    scale_fill_viridis_d(option = "turbo", name = "Bleaching Year") +
    labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", title = "(B) Event-Specific Dose-Response Trajectories") +
    coord_cartesian(ylim = c(0, 0.8), xlim = c(0, 15)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "right")

# Panel 3: Temporal Trend in Predicted Mortality at Fixed Thermal Stress Levels (4, 8, 12 DHW)
year_seq <- seq(1998, 2024, length.out = 100)
grid_secular <- expand.grid(
    year = year_seq,
    dhw_level = c(4, 8, 12)
) %>%
    mutate(
        Year_s = (year - 2016) / 10,
        MaxDHW_s = (dhw_level - mean_dhw) / sd2_dhw,
        secc3m_s = 0, cloudp_90_s = 0, prop_acropora_s = 0,
        winyear_sd_s = 0, histmDHW6_s = 0, mcur_90_s = 0,
        Dataset = factor("Manta Tow", levels = levels(dat.ml.expanded$Dataset)),
        Stress = factor(paste0(dhw_level, " °C-weeks"), levels = c("4 °C-weeks", "8 °C-weeks", "12 °C-weeks"))
    )

pp_secular <- posterior_epred(fit_brms_secular, newdata = grid_secular, re_formula = NA)
grid_secular$pred <- apply(pp_secular, 2, median)
grid_secular$pred_lo <- apply(pp_secular, 2, quantile, 0.05)
grid_secular$pred_hi <- apply(pp_secular, 2, quantile, 0.95)

p3_secular <- ggplot(grid_secular, aes(x = year, y = pred, colour = Stress, fill = Stress)) +
    geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.15, colour = NA) +
    geom_line(linewidth = 1.1) +
    geom_vline(xintercept = c(1998, 2002, 2016, 2017, 2020, 2022, 2024), linetype = "dotted", colour = "grey60") +
    scale_colour_manual(values = c("4 °C-weeks" = "#2c7bb6", "8 °C-weeks" = "#fdae61", "12 °C-weeks" = "#d7191c")) +
    scale_fill_manual(values = c("4 °C-weeks" = "#2c7bb6", "8 °C-weeks" = "#fdae61", "12 °C-weeks" = "#d7191c")) +
    labs(x = "Survey Year", y = "Predicted Mortality Proportion", title = "(C) Secular Vulnerability Trend at Fixed DHW (1998–2024)", colour = "Thermal Stress", fill = "Thermal Stress") +
    coord_cartesian(ylim = c(0, 0.7), xlim = c(1998, 2024)) +
    theme_bw(base_size = 11) +
    theme(legend.position = c(0.25, 0.82), legend.background = element_rect(fill = alpha('white', 0.8)))

# Panel 4: Distribution of Observed Heat Stress vs Mortality by Era
p4_dist <- ggplot(dat.ml.expanded, aes(x = MaxDHW_val, y = Mort.prop, colour = Era)) +
    geom_point(alpha = 0.45, size = 1.5) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 4), se = TRUE, alpha = 0.15, linewidth = 1) +
    scale_colour_manual(values = c("Historical (1998-2002)" = "#2b5c8f", "Modern Baseline (2016-2020)" = "#e66101", "Contemporary (2022-2024)" = "#5e3c99")) +
    facet_wrap(~Era, ncol = 3) +
    labs(x = "Max DHW (°C-weeks)", y = "Observed Mortality Proportion", title = "(D) Empirical GAM Trajectories Stratified by Era") +
    coord_cartesian(ylim = c(0, 0.8), xlim = c(0, 15)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey95"))

composite_temporal_plot <- (p1_era + p2_event) / (p3_secular + p4_dist)

dir.create("output/plots", showWarnings = FALSE, recursive = TRUE)
ggsave("output/plots/multidecadal_era_temporal_comparison.png", composite_temporal_plot, width = 13, height = 10, dpi = 200)
cat("Saved composite plot to output/plots/multidecadal_era_temporal_comparison.png!\n")
