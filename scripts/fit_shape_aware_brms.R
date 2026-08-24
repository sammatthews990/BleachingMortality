# Fit the shape-aware Bayesian candidate under leave-one-event-out validation,
# then refit on all events. Existing formal-model caches are never overwritten.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(posterior)
    library(readr)
})
source('scripts/shape_aware_model_helpers.R')

set.seed(202408L)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores(logical = FALSE)))

output_dir <- 'output/shape_aware_models'
model_dir <- 'output/models/shape_aware'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

cv_chains <- as.integer(Sys.getenv('SHAPE_BRMS_CHAINS', '2'))
cv_iter <- as.integer(Sys.getenv('SHAPE_BRMS_ITER', '2000'))
cv_warmup <- as.integer(Sys.getenv('SHAPE_BRMS_WARMUP', '1000'))
cv_adapt_delta <- as.numeric(Sys.getenv('SHAPE_BRMS_ADAPT_DELTA', '0.995'))
production_chains <- as.integer(Sys.getenv('SHAPE_BRMS_PRODUCTION_CHAINS', '4'))
production_iter <- as.integer(Sys.getenv('SHAPE_BRMS_PRODUCTION_ITER', '2400'))
production_warmup <- as.integer(Sys.getenv('SHAPE_BRMS_PRODUCTION_WARMUP', '1200'))
refit_failed_diagnostics <- identical(
    Sys.getenv('SHAPE_BRMS_REFIT_DIAGNOSTICS'), '1'
)

diagnostics_pass <- function(diagnostics) {
    diagnostics$divergences == 0 &&
        diagnostics$max_treedepth_hits == 0 &&
        diagnostics$max_rhat <= 1.01 &&
        diagnostics$min_bulk_ess >= 100 &&
        diagnostics$min_tail_ess >= 100
}

fit_shape_model <- function(template, formula, family, priors, data, seed,
                            chains, iter, warmup) {
    arguments <- list(
        formula = formula, data = data, family = family, prior = priors,
        backend = 'rstan', chains = chains,
        cores = min(chains, getOption('mc.cores', 1L)),
        iter = iter, warmup = warmup, seed = seed,
        control = list(adapt_delta = cv_adapt_delta, max_treedepth = 14),
        refresh = 0, silent = 2
    )
    if (is.null(template)) return(do.call(brm, arguments))
    update(
        template, newdata = data, recompile = FALSE,
        chains = chains, cores = min(chains, getOption('mc.cores', 1L)),
        iter = iter, warmup = warmup, seed = seed,
        control = list(adapt_delta = cv_adapt_delta, max_treedepth = 14),
        refresh = 0, silent = 2
    )
}

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

if (identical(Sys.getenv('SHAPE_BRMS_DRY_RUN'), '1')) {
    data <- load_shape_rows('ltmp')
    predictors <- shape_predictors('ltmp')
    prepared <- prepare_fold_predictors(data, data, predictors)
    formula <- make_shape_formula('ltmp', predictors)
    priors <- shape_priors('ltmp')
    print(get_prior(formula, prepared$analysis, family = shape_family('ltmp')))
    make_stancode(
        formula, prepared$analysis, family = shape_family('ltmp'),
        prior = priors
    )
    cat('Shape-aware BRMS dry run passed.\n')
    quit(save = 'no', status = 0)
}

programme_keys <- names(validation_files)
requested_programme <- trimws(Sys.getenv('SHAPE_PROGRAMME_KEY', ''))
if (nzchar(requested_programme)) programme_keys <- requested_programme

all_predictions <- tibble()
all_diagnostics <- tibble()

for (programme_key in programme_keys) {
    data <- load_shape_rows(programme_key)
    predictors <- shape_predictors(programme_key)
    formula <- make_shape_formula(programme_key, predictors)
    family <- shape_family(programme_key)
    priors <- shape_priors(programme_key)
    template <- NULL

    for (fold in sort(unique(data$event_year))) {
        assessment_index <- data$event_year == fold
        analysis <- data[!assessment_index, , drop = FALSE]
        assessment <- data[assessment_index, , drop = FALSE]
        prepared <- prepare_fold_predictors(analysis, assessment, predictors)
        cache_file <- file.path(
            model_dir,
            paste0('shape_brms_', programme_key, '_event_', fold, '.rds')
        )
        cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
        cache_current <- !is.null(cached) &&
            identical(cached$predictors, predictors) &&
            identical(cached$model_version, 'hinge_acropora_v1')
        if (cache_current && refit_failed_diagnostics) {
            cache_current <- diagnostics_pass(cached$diagnostics)
        }

        if (cache_current) {
            fit <- cached$fit
            predictions <- cached$predictions
            diagnostics <- cached$diagnostics
            cat('Loaded shape-aware BRMS:', programme_key, fold, '\n')
        } else {
            cat('Fitting shape-aware BRMS:', programme_key, fold, '\n')
            flush.console()
            fit <- fit_shape_model(
                template, formula, family, priors, prepared$analysis,
                seed = 202408L + fold, chains = cv_chains,
                iter = cv_iter, warmup = cv_warmup
            )
            predictions <- shape_prediction_components(
                fit, prepared$assessment, programme_key
            )
            diagnostics <- fit_diagnostics(fit)
            saveRDS(
                list(
                    fit = fit, preprocessing = prepared$preprocessing,
                    predictors = predictors, predictions = predictions,
                    diagnostics = diagnostics,
                    model_version = 'hinge_acropora_v1'
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
                    event_year, region_block,
                    ann_maxdhw, prop_acropora_pre, observed_pre_cover,
                    observed_mortality = mortality_prop
                ) |>
                bind_cols(predictions) |>
                mutate(
                    learner = 'shape_brms',
                    scheme = 'leave_one_event_out', fold = as.character(fold)
                )
        )
        all_diagnostics <- bind_rows(
            all_diagnostics,
            diagnostics |>
                mutate(
                    programme_key, fold = as.character(fold),
                    n_analysis = nrow(analysis),
                    n_assessment = nrow(assessment)
                )
        )
    }

    prepared <- prepare_fold_predictors(data, data, predictors)
    cache_file <- file.path(
        model_dir, paste0('shape_brms_', programme_key, '_production.rds')
    )
    cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
    cache_current <- !is.null(cached) &&
        identical(cached$predictors, predictors) &&
        identical(cached$model_version, 'hinge_acropora_v1')
    if (cache_current && refit_failed_diagnostics) {
        cache_current <- diagnostics_pass(cached$diagnostics)
    }
    if (!cache_current) {
        cat('Fitting production shape-aware BRMS:', programme_key, '\n')
        production_fit <- fit_shape_model(
            template, formula, family, priors, prepared$analysis,
            seed = 202408L, chains = production_chains,
            iter = production_iter, warmup = production_warmup
        )
        saveRDS(
            list(
                fit = production_fit,
                preprocessing = prepared$preprocessing,
                predictors = predictors,
                diagnostics = fit_diagnostics(production_fit),
                event_years = sort(unique(data$event_year)), n = nrow(data),
                model_version = 'hinge_acropora_v1'
            ),
            cache_file
        )
    }
}

if (any(!is.finite(all_predictions$predicted_mortality))) {
    stop('Shape-aware BRMS produced invalid predictions.')
}
write_csv(
    all_predictions,
    file.path(output_dir, 'shape_brms_event_predictions.csv')
)
write_csv(
    all_diagnostics,
    file.path(output_dir, 'shape_brms_diagnostics.csv')
)
write_csv(
    prediction_metrics(all_predictions),
    file.path(output_dir, 'shape_brms_metrics.csv')
)
print(prediction_metrics(all_predictions))
