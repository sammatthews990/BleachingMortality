# Fit and validate one shared reef-event ecological mortality response with
# separate manta, LTMP and MMP observation layers.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(posterior)
    library(readr)
})
source('src/lib/latent_mortality_model_helpers.R')

set.seed(202412L)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores(logical = FALSE)))

output_dir <- 'output/latent_mortality'
model_dir <- 'output/models/latent_mortality'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

formula <- make_latent_mortality_formula()
priors <- latent_mortality_priors()
data <- load_joint_compound_rows()
model_version <- 'latent_mortality_v3'

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

fit_one <- function(template, model_data, seed, chains, iter, warmup) {
    arguments <- list(
        formula = formula, data = model_data,
        family = zero_one_inflated_beta(), prior = priors,
        backend = 'rstan', chains = chains,
        cores = min(chains, getOption('mc.cores', 1L)),
        iter = iter, warmup = warmup, seed = seed,
        control = list(adapt_delta = 0.995, max_treedepth = 14),
        refresh = 0, silent = 2
    )
    if (is.null(template)) return(do.call(brm, arguments))
    update(
        template, newdata = model_data, recompile = FALSE,
        chains = chains, cores = min(chains, getOption('mc.cores', 1L)),
        iter = iter, warmup = warmup, seed = seed,
        control = list(adapt_delta = 0.995, max_treedepth = 14),
        refresh = 0, silent = 2
    )
}

prepared_all <- prepare_latent_fold(data, data)
if (identical(Sys.getenv('LATENT_BRMS_DRY_RUN'), '1')) {
    make_stancode(
        formula, prepared_all$analysis,
        family = zero_one_inflated_beta(), prior = priors
    )
    cat('Latent mortality BRMS dry run passed.\n')
    quit(save = 'no', status = 0)
}

cv_chains <- as.integer(Sys.getenv('LATENT_BRMS_CHAINS', '2'))
cv_iter <- as.integer(Sys.getenv('LATENT_BRMS_ITER', '1400'))
cv_warmup <- as.integer(Sys.getenv('LATENT_BRMS_WARMUP', '700'))
production_chains <- as.integer(Sys.getenv('LATENT_BRMS_PRODUCTION_CHAINS', '4'))
production_iter <- as.integer(Sys.getenv('LATENT_BRMS_PRODUCTION_ITER', '2400'))
production_warmup <- as.integer(Sys.getenv('LATENT_BRMS_PRODUCTION_WARMUP', '1200'))

schemes <- c('leave_one_event_out', 'reef_blocked_2024')
requested_scheme <- trimws(Sys.getenv('LATENT_BRMS_SCHEME', ''))
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
        prepared <- prepare_latent_fold(analysis, assessment)
        cache_file <- file.path(
            model_dir,
            paste0('latent_', scheme, '_', fold, '.rds')
        )
        cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
        cache_current <- !is.null(cached) &&
            identical(cached$model_version, model_version) &&
            identical(cached$predictors, latent_mortality_predictors)

        if (cache_current) {
            fit <- cached$fit
            predictions <- cached$predictions
            diagnostics <- cached$diagnostics
            cat('Loaded latent BRMS:', scheme, fold, '\n')
        } else {
            cat('Fitting latent BRMS:', scheme, fold, '\n')
            flush.console()
            seed <- 202412L + as.integer(fold) +
                if_else(scheme == 'reef_blocked_2024', 1000L, 0L)
            fit <- fit_one(
                template, prepared$analysis, seed,
                cv_chains, cv_iter, cv_warmup
            )
            predictions <- latent_prediction_components(
                fit, prepared$assessment
            )
            diagnostics <- fit_diagnostics(fit)
            saveRDS(
                list(
                    fit = fit, predictions = predictions,
                    preprocessing = prepared$preprocessing,
                    predictors = latent_mortality_predictors,
                    diagnostics = diagnostics,
                    model_version = model_version,
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
                    log_coastal_rain30, wqc_freqcc12,
                    cyc_maxHrs4mw, DISTURBANCE_TYPE, disturbance_text,
                    observed_mortality = mortality_prop
                ) |>
                bind_cols(predictions) |>
                mutate(
                    learner = 'shared_latent_brms', scheme = scheme,
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

production_file <- file.path(model_dir, 'latent_production.rds')
production <- if (file.exists(production_file)) readRDS(production_file) else NULL
production_current <- !is.null(production) &&
    identical(production$model_version, model_version) &&
    identical(production$predictors, latent_mortality_predictors)
if (!production_current) {
    cat('Fitting latent BRMS production model\n')
    production_fit <- fit_one(
        template, prepared_all$analysis, 202412L,
        production_chains, production_iter, production_warmup
    )
    production <- list(
        fit = production_fit,
        preprocessing = prepared_all$preprocessing,
        predictors = latent_mortality_predictors,
        diagnostics = fit_diagnostics(production_fit),
        model_version = model_version,
        n_observations = nrow(data),
        n_reef_events = nrow(prepared_all$analysis_events)
    )
    saveRDS(production, production_file)
}

write_csv(all_predictions, file.path(output_dir, 'latent_predictions.csv'))
write_csv(all_diagnostics, file.path(output_dir, 'latent_diagnostics.csv'))
write_csv(
    prepared_all$analysis_events,
    file.path(output_dir, 'latent_reef_event_rows.csv')
)
write_csv(
    production$diagnostics,
    file.path(output_dir, 'latent_production_diagnostics.csv')
)

production$loo <- loo(production$fit)
production$r2 <- bayes_R2(production$fit, re_formula = NA)
saveRDS(production, production_file)

print(
    all_predictions |>
        group_by(scheme, programme_key) |>
        summarise(
            n = n(),
            rmse = sqrt(mean(residual^2)),
            mae = mean(absolute_error),
            severe_mae = mean(
                absolute_error[observed_mortality >= 0.2]
            ),
            bias = mean(residual),
            .groups = 'drop'
        )
)
print(production$diagnostics)
