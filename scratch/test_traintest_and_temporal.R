suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(brms)
  library(gbm)
  library(ranger)
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

# Define Era & EventYear
dat.ml.expanded <- dat.ml.expanded %>%
    mutate(
        Era = case_when(
            DHWYear <= 2002 ~ "Historical (1998-2002)",
            DHWYear >= 2016 & DHWYear <= 2020 ~ "Modern Baseline (2016-2020)",
            DHWYear >= 2022 ~ "Contemporary (2022-2024)",
            TRUE ~ "Intermediate"
        ),
        Era = factor(Era, levels = c("Historical (1998-2002)", "Modern Baseline (2016-2020)", "Contemporary (2022-2024)")),
        EventYear = as.character(DHWYear),
        EventYear = factor(EventYear, levels = c("1998", "2002", "2016", "2017", "2020", "2022", "2024"))
    )

# Prepare scaled variables with all predictors separate
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
sd2_dhw  <- 2 * sd(dat.ml.expanded$MaxDHW_val, na.rm=TRUE)

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

cat("=========================================================\n")
cat("1. FITTING FULL TEMPORAL & ERA BAYESIAN BETA GLMMS\n")
cat("=========================================================\n")

# A. Secular Model with All 9 Separate Covariates + DHW * Year
cache_secular_brms <- "output/models/brms_secular_beta.rds"
if (file.exists(cache_secular_brms)) {
    fit_brms_secular <- readRDS(cache_secular_brms)
} else {
    fit_brms_secular <- brm(
        Mort.prop.nudge ~ MaxDHW_s * Year_s + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + 
                          Dataset + (1 | ReefName_clean),
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
cat("Loaded brms_secular_beta.rds\n")

# B. Era Model with All Separate Covariates
cache_era_brms <- "output/models/brms_era_beta.rds"
if (file.exists(cache_era_brms)) {
    fit_brms_era <- readRDS(cache_era_brms)
} else {
    fit_brms_era <- brm(
        Mort.prop.nudge ~ MaxDHW_s * Era + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + 
                          Dataset + (1 | ReefName_clean),
        data = dat.ml.expanded %>% filter(!is.na(Era)), family = Beta(link = "logit"),
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
cat("Loaded brms_era_beta.rds\n")

# C. Event Year Model
cache_event_brms <- "output/models/brms_eventyear_beta.rds"
dat_events <- dat.ml.expanded %>% filter(!is.na(EventYear))
if (file.exists(cache_event_brms)) {
    fit_brms_event <- readRDS(cache_event_brms)
} else {
    fit_brms_event <- brm(
        Mort.prop.nudge ~ MaxDHW_s * EventYear + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + 
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
cat("Loaded brms_eventyear_beta.rds\n")

cat("\n=========================================================\n")
cat("2. TRAIN-TEST SPLIT: TRAIN <= 2022 (N=", sum(dat.ml.expanded$DHWYear <= 2022), ") -> PREDICT 2024 (N=", sum(dat.ml.expanded$DHWYear == 2024), ")\n")
cat("=========================================================\n")

dat_train <- dat.ml.expanded %>% filter(DHWYear <= 2022)
dat_test  <- dat.ml.expanded %>% filter(DHWYear == 2024)

# 1. Fit Bayesian Beta GLMM on Training Data
cache_brms_train <- "output/models/brms_train_beta_pre2024.rds"
if (file.exists(cache_brms_train)) {
    fit_brms_train <- readRDS(cache_brms_train)
} else {
    fit_brms_train <- brm(
        Mort.prop.nudge ~ MaxDHW_s + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + 
                          Dataset + (1 | ReefName_clean),
        data = dat_train, family = Beta(link = "logit"),
        prior = c(
            prior(normal(-1, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")
        ),
        chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 42,
        control = list(adapt_delta = 0.95)
    )
    saveRDS(fit_brms_train, cache_brms_train)
}
cat("Loaded fit_brms_train\n")

# 2. Fit Quasibinomial GLM on Training Data
fit_qbin_train <- glm(
    Mort.prop ~ MaxDHW_s + secc3m_s + cloudp_90_s + prop_acropora_s + 
                winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + Dataset,
    data = dat_train, family = quasibinomial(link = "logit")
)

# 3. Fit Random Forest on Training Data (via ranger)
set.seed(42)
rf_train_formula <- Mort.prop ~ MaxDHW_s + secc3m_s + cloudp_90_s + prop_acropora_s + 
                               winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + Dataset
fit_rf_train <- ranger(rf_train_formula, data = dat_train, num.trees = 500, importance = "impurity")

# 4. Fit Boosted Regression Trees (BRT/GBM) on Training Data
fit_brt_train <- gbm(
    rf_train_formula, data = dat_train, distribution = "gaussian",
    n.trees = 1500, interaction.depth = 3, shrinkage = 0.005, cv.folds = 5, n.cores = 2
)
best_trees <- gbm.perf(fit_brt_train, plot.it = FALSE, method = "cv")

# 5. Fit Bayesian Beta GLMM with Secular Trend (DHW * Year) on Training Data
cache_brms_train_temporal <- "output/models/brms_train_temporal_pre2024.rds"
if (file.exists(cache_brms_train_temporal)) {
    fit_brms_train_temporal <- readRDS(cache_brms_train_temporal)
} else {
    fit_brms_train_temporal <- brm(
        Mort.prop.nudge ~ MaxDHW_s * Year_s + secc3m_s + cloudp_90_s + prop_acropora_s + 
                          winyear_sd_s + histmDHW6_s + mcur_90_s + winyear_mean_s + yrsince6_s + 
                          Dataset + (1 | ReefName_clean),
        data = dat_train, family = Beta(link = "logit"),
        prior = c(
            prior(normal(-1, 1.5), class = "Intercept"),
            prior(normal(0, 1), class = "b"),
            prior(exponential(1), class = "sd")
        ),
        chains = 4, iter = 2000, warmup = 1000, cores = 4, seed = 42,
        control = list(adapt_delta = 0.95)
    )
    saveRDS(fit_brms_train_temporal, cache_brms_train_temporal)
}
cat("Loaded fit_brms_train_temporal\n")

cat("\n=========================================================\n")
cat("3. PREDICTING 2024 EVENT OUT-OF-SAMPLE AND EVALUATING\n")
cat("=========================================================\n")

# Out-of-sample predictions
pp_brms <- posterior_epred(fit_brms_train, newdata = dat_test, re_formula = NA, allow_new_levels = TRUE)
pred_brms_med <- apply(pp_brms, 2, median)
pred_brms_lo  <- apply(pp_brms, 2, quantile, 0.05)
pred_brms_hi  <- apply(pp_brms, 2, quantile, 0.95)

pp_brms_temp <- posterior_epred(fit_brms_train_temporal, newdata = dat_test, re_formula = NA, allow_new_levels = TRUE)
pred_brms_temp_med <- apply(pp_brms_temp, 2, median)

pred_qbin <- predict(fit_qbin_train, newdata = dat_test, type = "response")
pred_rf   <- predict(fit_rf_train, data = dat_test)$predictions
pred_brt  <- predict(fit_brt_train, newdata = dat_test, n.trees = best_trees)

obs_2024 <- dat_test$Mort.prop

calc_eval <- function(obs, pred, name) {
    res <- obs - pred
    rmse <- sqrt(mean(res^2))
    mae  <- mean(abs(res))
    r_val <- cor(obs, pred)
    rho_val <- cor(obs, pred, method = "spearman")
    # Calibration slope & intercept
    cal_m <- lm(obs ~ pred)
    alpha_cal <- coef(cal_m)[1]
    beta_cal  <- coef(cal_m)[2]
    data.frame(
        Model = name,
        RMSE = round(rmse, 4),
        MAE = round(mae, 4),
        Pearson_r = round(r_val, 4),
        Spearman_rho = round(rho_val, 4),
        Cal_Intercept = round(alpha_cal, 4),
        Cal_Slope = round(beta_cal, 4)
    )
}

eval_df <- bind_rows(
    calc_eval(obs_2024, pred_brms_med, "Bayesian Beta GLMM (Pre-2024)"),
    calc_eval(obs_2024, pred_brms_temp_med, "Bayesian Beta + Temporal Trend (Pre-2024)"),
    calc_eval(obs_2024, pred_qbin, "Quasibinomial GLM (Pre-2024)"),
    calc_eval(obs_2024, pred_rf, "Random Forest (Pre-2024)"),
    calc_eval(obs_2024, pred_brt, "Boosted Regression Trees (Pre-2024)")
)

# 90% CI Coverage for Bayesian model
coverage_90 <- mean(obs_2024 >= pred_brms_lo & obs_2024 <= pred_brms_hi)
eval_df$`Coverage_90_CI` <- c(paste0(round(coverage_90 * 100, 1), "%"), "-", "-", "-", "-")

print(as.data.frame(eval_df))

# Save metrics
saveRDS(eval_df, "output/models/traintest_2024_evaluation_metrics.rds")

cat("\n=========================================================\n")
cat("4. GENERATING 4-PANEL 2024 OUT-OF-SAMPLE EVALUATION FIGURE\n")
cat("=========================================================\n")

# Prepare dataframe for plotting
df_eval_plot <- dat_test %>%
    mutate(
        pred_brms = pred_brms_med,
        pred_lo = pred_brms_lo,
        pred_hi = pred_brms_hi,
        pred_qbin = pred_qbin,
        pred_rf = pred_rf,
        pred_brt = pred_brt,
        residual_brms = Mort.prop - pred_brms
    )

# Panel A: Observed vs Predicted 2024 Calibration (Multi-Model)
df_comp_long <- bind_rows(
    data.frame(Obs = obs_2024, Pred = pred_brms_med, Model = "Bayesian Beta GLMM", Dataset = dat_test$Dataset),
    data.frame(Obs = obs_2024, Pred = pred_qbin, Model = "Quasibinomial GLM", Dataset = dat_test$Dataset),
    data.frame(Obs = obs_2024, Pred = pred_rf, Model = "Random Forest", Dataset = dat_test$Dataset),
    data.frame(Obs = obs_2024, Pred = pred_brt, Model = "Boosted Regression Trees", Dataset = dat_test$Dataset)
)

p_eval_A <- ggplot(df_comp_long, aes(x = Pred, y = Obs, colour = Model)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey40", linewidth = 0.9) +
    geom_point(alpha = 0.45, size = 1.3) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1.1) +
    scale_colour_manual(values = c("Bayesian Beta GLMM" = "#2b5c8f", "Quasibinomial GLM" = "#e66101", "Random Forest" = "#33a02c", "Boosted Regression Trees" = "#e31a1c")) +
    facet_wrap(~Model, ncol = 2) +
    labs(
        x = "Out-of-Sample Predicted Mortality (Trained on 1998–2022)",
        y = "Observed 2024 Mortality Proportion",
        title = "(A) 2024 Out-of-Sample Calibration by Model Framework"
    ) +
    coord_cartesian(xlim = c(0, 0.7), ylim = c(0, 0.8)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey95"))

# Panel B: Pre-2024 Bayesian Dose-Response Curve vs 2024 Observed Points
dhw_seq <- seq(0, 15, length.out = 200)
grid_pred_2024 <- data.frame(
    MaxDHW_val = dhw_seq,
    MaxDHW_s = (dhw_seq - mean_dhw) / sd2_dhw,
    secc3m_s = 0, cloudp_90_s = 0, prop_acropora_s = 0,
    winyear_sd_s = 0, histmDHW6_s = 0, mcur_90_s = 0,
    winyear_mean_s = 0, yrsince6_s = 0,
    Dataset = factor("Manta Tow", levels = levels(dat.ml.expanded$Dataset))
)
pp_curve_train <- posterior_epred(fit_brms_train, newdata = grid_pred_2024, re_formula = NA)
grid_pred_2024$pred <- apply(pp_curve_train, 2, median)
grid_pred_2024$pred_lo <- apply(pp_curve_train, 2, quantile, 0.05)
grid_pred_2024$pred_hi <- apply(pp_curve_train, 2, quantile, 0.95)

max_dhw_train <- max(dat_train$MaxDHW_val, na.rm = TRUE)

p_eval_B <- ggplot() +
    annotate("rect", xmin = max_dhw_train, xmax = 15, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.4) +
    annotate("text", x = (max_dhw_train + 15)/2, y = 0.05, label = "2024 Extrapolation Zone", fontface = "italic", colour = "grey40", size = 3) +
    geom_ribbon(data = grid_pred_2024, aes(x = MaxDHW_val, ymin = pred_lo, ymax = pred_hi), fill = "#2b5c8f", alpha = 0.2) +
    geom_line(data = grid_pred_2024, aes(x = MaxDHW_val, y = pred), colour = "#2b5c8f", linewidth = 1.2) +
    geom_point(data = dat_test, aes(x = MaxDHW_val, y = Mort.prop, colour = Dataset), size = 1.8, alpha = 0.7) +
    scale_colour_brewer(palette = "Set2") +
    labs(
        x = "Max DHW (°C-weeks)", y = "Mortality Proportion",
        colour = "2024 Survey Dataset",
        title = "(B) Pre-2024 Model Forecast vs Observed 2024 Surveys"
    ) +
    coord_cartesian(xlim = c(0, 15), ylim = c(0, 0.8)) +
    theme_bw(base_size = 11) +
    theme(legend.position = c(0.3, 0.82), legend.background = element_rect(fill = alpha("white", 0.85)))

# Panel C: Regional Residual Error Distribution (Spatial Transferability)
df_eval_plot <- df_eval_plot %>%
    mutate(
        Sector_clean = coalesce(as.character(SECT_NAME), as.character(SECTOR), as.character(Region), "Other GBR")
    )


p_eval_C <- ggplot(df_eval_plot, aes(x = Sector_clean, y = residual_brms, fill = Sector_clean)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_boxplot(alpha = 0.7, outlier.size = 1.2) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 1.2) +
    scale_fill_brewer(palette = "Dark2") +
    labs(
        x = "GBR Sector", y = "Prediction Error (Observed - Predicted)",
        title = "(C) 2024 Regional Forecast Errors (Spatial Generalization)"
    ) +
    coord_cartesian(ylim = c(-0.4, 0.6)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))

# Panel D: Stress-Dependent Residuals (Extrapolation Bias)
p_eval_D <- ggplot(df_eval_plot, aes(x = MaxDHW_val, y = residual_brms)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(aes(colour = Dataset), size = 1.8, alpha = 0.7) +
    geom_smooth(method = "loess", colour = "firebrick", fill = "firebrick", alpha = 0.15, linewidth = 1) +
    scale_colour_brewer(palette = "Set2") +
    labs(
        x = "Max DHW (°C-weeks)", y = "Prediction Error (Observed - Predicted)",
        title = "(D) Forecast Residuals vs 2024 Thermal Stress"
    ) +
    coord_cartesian(xlim = c(0, 15), ylim = c(-0.4, 0.6)) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

composite_traintest_2024 <- (p_eval_A + p_eval_B) / (p_eval_C + p_eval_D)

ggsave("output/plots/traintest_2024_out_of_sample_evaluation.png", composite_traintest_2024, width = 13, height = 10.5, dpi = 200)
cat("Saved composite plot to output/plots/traintest_2024_out_of_sample_evaluation.png!\n")
