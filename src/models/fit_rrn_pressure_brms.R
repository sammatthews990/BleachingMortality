# Formal zero-one-inflated beta models for the selected RRN pressure set.
# The two variants differ only by modelled COTS pressure, allowing a direct
# LOO and Bayesian-R2 comparison after the blocked BRT prediction screen.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(posterior)
    library(readr)
})
source('src/lib/joint_compound_model_helpers.R')

set.seed(202411L)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores(logical = FALSE)))

output_dir <- 'output/rrn_pressure_assessment'
model_dir <- 'output/models/rrn_pressure'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

model_variant <- trimws(Sys.getenv('RRN_BRMS_VARIANT', 'weather_rrn'))
allowed_variants <- c('weather_rrn', 'weather_rrn_cots')
if (!model_variant %in% allowed_variants) {
    stop('RRN_BRMS_VARIANT must be one of: ', paste(allowed_variants, collapse = ', '))
}

selected_predictors <- c(
    joint_compound_predictors,
    'wqc_freqcc12', 'log1p_cyc_maxHrs4mw'
)
if (model_variant == 'weather_rrn_cots') {
    selected_predictors <- c(selected_predictors, 'log1p_cot_idwmeanpertow')
}

mean_terms <- c(
    '0 + programme_factor',
    paste0(selected_predictors, '_z'),
    'dhw_excess4_z:prop_acropora_pre_z',
    'dhw_excess8_z:prop_acropora_pre_z',
    'ann_maxdhw_z:log_coastal_rain30_z',
    'ann_maxdhw_z:wqc_freqcc12_z',
    '(1 | reef_event_effect)'
)
boundary_terms <- c(
    '0 + programme_factor',
    'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z',
    'prop_acropora_pre_z', 'observed_pre_cover_z',
    'log_coastal_rain30_z', 'wqc_freqcc12_z',
    'log1p_cyc_maxHrs4mw_z'
)
if (model_variant == 'weather_rrn_cots') {
    boundary_terms <- c(boundary_terms, 'log1p_cot_idwmeanpertow_z')
}

formula <- bf(
    as.formula(paste('mortality_prop ~', paste(mean_terms, collapse = ' + '))),
    phi = ~ 0 + programme_factor,
    zoi = as.formula(paste('~', paste(boundary_terms, collapse = ' + '))),
    coi = ~ 0 + programme_factor
)

programme_coefficients <- paste0('programme_factor', joint_programmes)
priors <- c(
    prior(normal(0, 0.5), class = 'b'),
    prior(exponential(1), class = 'sd'),
    prior(normal(0, 0.7), class = 'b', dpar = 'zoi'),
    prior(normal(2, 0.8), class = 'b', dpar = 'phi'),
    prior(normal(-4, 1.2), class = 'b', dpar = 'coi'),
    prior(normal(0.15, 0.2), class = 'b', coef = 'ann_maxdhw_z'),
    prior(normal(0.2, 0.25), class = 'b', coef = 'dhw_excess4_z'),
    prior(normal(0.35, 0.3), class = 'b', coef = 'dhw_excess8_z'),
    prior(
        normal(0.2, 0.25), class = 'b',
        coef = 'ann_maxdhw_z:log_coastal_rain30_z'
    ),
    prior(
        normal(0, 0.3), class = 'b',
        coef = 'ann_maxdhw_z:wqc_freqcc12_z'
    )
)
for (coefficient in programme_coefficients) {
    priors <- c(
        priors,
        set_prior('normal(-2, 1)', class = 'b', coef = coefficient)
    )
}

data <- load_joint_compound_rows()
prepared <- prepare_joint_predictors(data, data)
model_file <- file.path(
    model_dir, paste0('rrn_pressure_brms_', model_variant, '.rds')
)
model_version <- paste0('rrn_pressure_brms_', model_variant, '_v2')
cached <- if (file.exists(model_file)) readRDS(model_file) else NULL

if (identical(Sys.getenv('RRN_BRMS_DRY_RUN'), '1')) {
    make_stancode(
        formula, prepared$analysis,
        family = zero_one_inflated_beta(), prior = priors
    )
    cat('RRN pressure BRMS ', model_variant, ' dry run passed.\n', sep = '')
    quit(save = 'no', status = 0)
}

cache_current <- !is.null(cached) &&
    identical(cached$model_version, model_version) &&
    identical(cached$predictors, selected_predictors)
if (cache_current) {
    fit <- cached$fit
    cat('Loaded RRN pressure BRMS ', model_variant, '\n', sep = '')
} else {
    fit <- brm(
        formula = formula,
        data = prepared$analysis,
        family = zero_one_inflated_beta(),
        prior = priors,
        backend = 'rstan',
        chains = as.integer(Sys.getenv('RRN_BRMS_CHAINS', '4')),
        cores = min(
            as.integer(Sys.getenv('RRN_BRMS_CHAINS', '4')),
            getOption('mc.cores', 1L)
        ),
        iter = as.integer(Sys.getenv('RRN_BRMS_ITER', '2400')),
        warmup = as.integer(Sys.getenv('RRN_BRMS_WARMUP', '1200')),
        seed = 202411L,
        control = list(adapt_delta = 0.995, max_treedepth = 14),
        refresh = 50
    )
    cached <- list(
        fit = fit, preprocessing = prepared$preprocessing,
        predictors = selected_predictors,
        model_version = model_version, n = nrow(data)
    )
}

draws <- posterior::as_draws_df(fit)
effect_variables <- c(
    rainfall_mean = 'b_log_coastal_rain30_z',
    wq_mean = 'b_wqc_freqcc12_z',
    cyclone_mean = 'b_log1p_cyc_maxHrs4mw_z',
    dhw_by_rainfall_mean = 'b_ann_maxdhw_z:log_coastal_rain30_z',
    dhw_by_wq_mean = 'b_ann_maxdhw_z:wqc_freqcc12_z',
    wq_zero_one = 'b_zoi_wqc_freqcc12_z',
    cyclone_zero_one = 'b_zoi_log1p_cyc_maxHrs4mw_z'
)
if (model_variant == 'weather_rrn_cots') {
    effect_variables <- c(
        effect_variables,
        cots_mean = 'b_log1p_cot_idwmeanpertow_z',
        cots_zero_one = 'b_zoi_log1p_cot_idwmeanpertow_z'
    )
}
missing_draws <- setdiff(effect_variables, names(draws))
if (length(missing_draws) > 0L) {
    stop('Missing expected BRMS draws: ', paste(missing_draws, collapse = ', '))
}

effects <- bind_rows(lapply(names(effect_variables), function(effect_name) {
    values <- as.numeric(draws[[effect_variables[[effect_name]]]])
    tibble(
        model = model_variant,
        effect = effect_name,
        posterior_mean = mean(values),
        posterior_sd = sd(values),
        q025 = quantile(values, 0.025),
        q50 = quantile(values, 0.5),
        q975 = quantile(values, 0.975),
        probability_positive = mean(values > 0)
    )
}))
write_csv(
    effects,
    file.path(output_dir, paste0('rrn_brms_effects_', model_variant, '.csv'))
)

summary_draws <- posterior::summarise_draws(
    posterior::as_draws_array(fit),
    rhat = posterior::rhat,
    ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
)
nuts <- nuts_params(fit)
diagnostics <- tibble(
    model = model_variant,
    max_rhat = max(summary_draws$rhat, na.rm = TRUE),
    min_bulk_ess = min(summary_draws$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(summary_draws$ess_tail, na.rm = TRUE),
    divergences = sum(nuts$Value[nuts$Parameter == 'divergent__']),
    max_treedepth_hits = sum(
        nuts$Value[nuts$Parameter == 'treedepth__'] >= 14
    )
)
write_csv(
    diagnostics,
    file.path(output_dir, paste0('rrn_brms_diagnostics_', model_variant, '.csv'))
)

loo_result <- loo(fit)
r2_draws <- as.numeric(bayes_R2(fit, re_formula = NA, summary = FALSE))
fit_metrics <- tibble(
    model = model_variant,
    elpd_loo = loo_result$estimates['elpd_loo', 'Estimate'],
    elpd_loo_se = loo_result$estimates['elpd_loo', 'SE'],
    looic = loo_result$estimates['looic', 'Estimate'],
    looic_se = loo_result$estimates['looic', 'SE'],
    bayes_r2_mean = mean(r2_draws),
    bayes_r2_q025 = quantile(r2_draws, 0.025),
    bayes_r2_q975 = quantile(r2_draws, 0.975)
)
write_csv(
    fit_metrics,
    file.path(output_dir, paste0('rrn_brms_fit_', model_variant, '.csv'))
)

cached$loo <- loo_result
cached$fit_metrics <- fit_metrics
saveRDS(cached, model_file)

print(effects)
print(diagnostics)
print(fit_metrics)
