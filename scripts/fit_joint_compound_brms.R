# Fit the joint LTMP--manta--MMP compound-stress model. Event-held-out folds
# test transfer to new bleaching events; 2024 reef blocks test spatial transfer
# of the novel rainfall/calm mechanism while retaining other 2024 reefs.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(posterior)
    library(readr)
})
source('scripts/joint_compound_model_helpers.R')

set.seed(202409L)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores(logical = FALSE)))

output_dir <- 'output/joint_compound_models'
model_dir <- 'output/models/joint_compound'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

cv_chains <- as.integer(Sys.getenv('JOINT_BRMS_CHAINS', '2'))
cv_iter <- as.integer(Sys.getenv('JOINT_BRMS_ITER', '1400'))
cv_warmup <- as.integer(Sys.getenv('JOINT_BRMS_WARMUP', '700'))
production_chains <- as.integer(Sys.getenv('JOINT_BRMS_PRODUCTION_CHAINS', '4'))
production_iter <- as.integer(Sys.getenv('JOINT_BRMS_PRODUCTION_ITER', '2400'))
production_warmup <- as.integer(Sys.getenv('JOINT_BRMS_PRODUCTION_WARMUP', '1200'))
adapt_delta <- as.numeric(Sys.getenv('JOINT_BRMS_ADAPT_DELTA', '0.995'))
refit_failed <- identical(Sys.getenv('JOINT_BRMS_REFIT_DIAGNOSTICS'), '1')

fit_diagnostics <- function(fit) {
    draws <- posterior::summarise_draws(
        posterior::as_draws_array(fit),
        rhat = posterior::rhat,
        ess_bulk = posterior::ess_bulk,
        ess_tail = posterior::ess_tail
    )
    nuts <- nuts_params(fit)
    tibble(
        max_rhat = max(draws$rhat, na.rm = TRUE),
        min_bulk_ess = min(draws$ess_bulk, na.rm = TRUE),
        min_tail_ess = min(draws$ess_tail, na.rm = TRUE),
        divergences = sum(nuts$Value[nuts$Parameter == 'divergent__']),
        max_treedepth_hits = sum(
            nuts$Value[nuts$Parameter == 'treedepth__'] >= 14
        )
    )
}

diagnostics_pass <- function(x) {
    x$max_rhat <= 1.01 && x$min_bulk_ess >= 100 &&
        x$min_tail_ess >= 100 && x$divergences == 0 &&
        x$max_treedepth_hits == 0
}

fit_model <- function(template, formula, priors, data, seed,
                      chains, iter, warmup) {
    arguments <- list(
        formula = formula, data = data,
        family = zero_one_inflated_beta(), prior = priors,
        backend = 'rstan', chains = chains,
        cores = min(chains, getOption('mc.cores', 1L)),
        iter = iter, warmup = warmup, seed = seed,
        control = list(adapt_delta = adapt_delta, max_treedepth = 14),
        refresh = 0, silent = 2
    )
    if (is.null(template)) return(do.call(brm, arguments))
    update(
        template, newdata = data, recompile = FALSE,
        chains = chains, cores = min(chains, getOption('mc.cores', 1L)),
        iter = iter, warmup = warmup, seed = seed,
        control = list(adapt_delta = adapt_delta, max_treedepth = 14),
        refresh = 0, silent = 2
    )
}

data <- load_joint_compound_rows()
formula <- make_joint_compound_formula()
priors <- joint_compound_priors()
model_version <- 'joint_compound_weather_v2'

if (identical(Sys.getenv('JOINT_BRMS_DRY_RUN'), '1')) {
    prepared <- prepare_joint_predictors(data, data)
    make_stancode(
        formula, prepared$analysis,
        family = zero_one_inflated_beta(), prior = priors
    )
    cat('Joint compound BRMS dry run passed.\n')
    quit(save = 'no', status = 0)
}

schemes <- c('leave_one_event_out', 'reef_blocked_2024')
requested_scheme <- trimws(Sys.getenv('JOINT_BRMS_SCHEME', ''))
if (nzchar(requested_scheme)) schemes <- requested_scheme

all_predictions <- tibble()
all_diagnostics <- tibble()
template <- NULL

for (scheme in schemes) {
    folds <- if (scheme == 'leave_one_event_out') {
        sort(unique(data$event_year))
    } else {
        sort(unique(data$joint_reef_fold[data$event_year == 2024L]))
    }
    for (fold in folds) {
        assessment_index <- joint_assessment_index(data, scheme, fold)
        analysis_index <- joint_analysis_index(data, scheme, fold)
        analysis <- data[analysis_index, , drop = FALSE]
        assessment <- data[assessment_index, , drop = FALSE]
        prepared <- prepare_joint_predictors(analysis, assessment)
        cache_file <- file.path(
            model_dir,
            paste0('joint_brms_v2_', scheme, '_', fold, '.rds')
        )
        cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
        cache_current <- !is.null(cached) &&
            identical(cached$model_version, model_version) &&
            identical(cached$predictors, joint_compound_predictors)
        if (cache_current && refit_failed) {
            cache_current <- diagnostics_pass(cached$diagnostics)
        }

        if (cache_current) {
            fit <- cached$fit
            predictions <- cached$predictions
            diagnostics <- cached$diagnostics
            cat('Loaded joint BRMS:', scheme, fold, '\n')
        } else {
            cat('Fitting joint BRMS:', scheme, fold, '\n')
            flush.console()
            seed <- 202409L + as.integer(fold) +
                if_else(scheme == 'reef_blocked_2024', 1000L, 0L)
            fit <- fit_model(
                template, formula, priors, prepared$analysis, seed,
                cv_chains, cv_iter, cv_warmup
            )
            predictions <- joint_prediction_components(
                fit, prepared$assessment
            )
            diagnostics <- fit_diagnostics(fit)
            saveRDS(
                list(
                    fit = fit, predictions = predictions,
                    preprocessing = prepared$preprocessing,
                    predictors = joint_compound_predictors,
                    diagnostics = diagnostics, model_version = model_version,
                    scheme = scheme, fold = fold
                ),
                cache_file
            )
        }
        template <- fit

        all_predictions <- bind_rows(
            all_predictions,
            assessment |>
                transmute(
                    programme_key, source_observation_id, ReefID, ReefName,
                    event_year, region_block, joint_reef_fold,
                    ann_maxdhw, prop_acropora_pre, observed_pre_cover,
                    log_coastal_rain30, era5_wind_mean,
                    era5_wind_calm_fraction, era5_coastal_distance_km,
                    DISTURBANCE_TYPE, description,
                    observed_mortality = mortality_prop
                ) |>
                bind_cols(predictions) |>
                mutate(
                    learner = 'joint_compound_brms', scheme = scheme,
                    fold = as.character(fold),
                    residual = observed_mortality - predicted_mortality,
                    absolute_error = abs(residual)
                )
        )
        all_diagnostics <- bind_rows(
            all_diagnostics,
            diagnostics |>
                mutate(
                    scheme = scheme, fold = as.character(fold),
                    n_analysis = nrow(analysis),
                    n_assessment = nrow(assessment)
                )
        )
    }
}

prepared <- prepare_joint_predictors(data, data)
production_file <- file.path(model_dir, 'joint_brms_v2_production.rds')
production <- if (file.exists(production_file)) readRDS(production_file) else NULL
production_current <- !is.null(production) &&
    identical(production$model_version, model_version) &&
    identical(production$predictors, joint_compound_predictors)
if (production_current && refit_failed) {
    production_current <- diagnostics_pass(production$diagnostics)
}
if (!production_current) {
    cat('Fitting joint BRMS production model\n')
    production_fit <- fit_model(
        template, formula, priors, prepared$analysis, 202409L,
        production_chains, production_iter, production_warmup
    )
    production <- list(
        fit = production_fit,
        preprocessing = prepared$preprocessing,
        predictors = joint_compound_predictors,
        diagnostics = fit_diagnostics(production_fit),
        model_version = model_version, n = nrow(data)
    )
    saveRDS(production, production_file)
}

write_csv(data, file.path(output_dir, 'joint_compound_rows.csv'))
prediction_file <- file.path(output_dir, 'joint_brms_predictions.csv')
diagnostic_file <- file.path(output_dir, 'joint_brms_diagnostics.csv')
if (file.exists(prediction_file)) {
    retained <- read_csv(prediction_file, show_col_types = FALSE) |>
        mutate(
            source_observation_id = as.character(source_observation_id),
            fold = as.character(fold)
        ) |>
        filter(!scheme %in% unique(all_predictions$scheme))
    all_predictions <- all_predictions |>
        mutate(source_observation_id = as.character(source_observation_id))
    all_predictions <- bind_rows(retained, all_predictions)
}
if (file.exists(diagnostic_file)) {
    retained <- read_csv(diagnostic_file, show_col_types = FALSE) |>
        mutate(fold = as.character(fold)) |>
        filter(!scheme %in% unique(all_diagnostics$scheme))
    all_diagnostics <- bind_rows(retained, all_diagnostics)
}
write_csv(all_predictions, prediction_file)
write_csv(all_diagnostics, diagnostic_file)
print(all_predictions |>
    group_by(scheme, programme_key) |>
    summarise(
        n = n(), rmse = sqrt(mean(residual^2)),
        mae = mean(absolute_error),
        severe_mae = mean(absolute_error[observed_mortality >= 0.2]),
        bias = mean(residual), .groups = 'drop'
    ))
