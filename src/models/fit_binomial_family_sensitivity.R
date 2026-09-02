# Bayesian binomial-OLRE sensitivity model under leave-one-event-out folds.
#
# Mortality proportions are not literal colony-level binomial counts for all
# programmes. We therefore use 100 pseudo-trials only as a response-family
# sensitivity and include an observation-level random effect to absorb extra-
# binomial variation. The programme-specific beta-family model remains the
# primary bounded likelihood unless this sensitivity improves event-held-out
# prediction.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(posterior)
    library(readr)
})
source('src/lib/formal_model_helpers.R')

set.seed(202408L)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores(logical = FALSE)))

output_dir <- 'output/formal_models'
model_dir <- 'output/models/formal_binomial_sensitivity'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

chains <- as.integer(Sys.getenv('BINOMIAL_BRMS_CHAINS', '4'))
iterations <- as.integer(Sys.getenv('BINOMIAL_BRMS_ITER', '2400'))
warmup <- as.integer(Sys.getenv('BINOMIAL_BRMS_WARMUP', '1200'))
adapt_delta <- as.numeric(Sys.getenv('BINOMIAL_BRMS_ADAPT_DELTA', '0.995'))

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
        divergences = sum(
            nuts$Value[nuts$Parameter == 'divergent__'], na.rm = TRUE
        ),
        max_treedepth_hits = sum(
            nuts$Value[nuts$Parameter == 'treedepth__'] >= 14,
            na.rm = TRUE
        )
    )
}

diagnostics_pass <- function(diagnostics) {
    diagnostics$divergences == 0 &&
        diagnostics$max_treedepth_hits == 0 &&
        diagnostics$max_rhat <= 1.01 &&
        diagnostics$min_bulk_ess >= 100 &&
        diagnostics$min_tail_ess >= 100
}

all_predictions <- tibble()
all_diagnostics <- tibble()

for (programme_key in names(validation_files)) {
    data <- load_programme_rows(programme_key) |>
        mutate(
            pseudo_trials = 100L,
            mortality_success = as.integer(round(mortality_prop * 100)),
            observation_effect = factor(source_observation_id)
        )
    predictors <- programme_predictors(programme_key)
    formula <- bf(as.formula(paste(
        'mortality_success | trials(pseudo_trials) ~',
        brms_rhs(predictors), '+ (1 | observation_effect)'
    )))
    priors <- c(
        prior(normal(0, 0.5), class = 'b'),
        prior(normal(-1.5, 1), class = 'Intercept'),
        prior(exponential(1), class = 'sd')
    )
    template_fit <- NULL

    for (fold in fold_values(data, 'leave_one_event_out')) {
        assessment_index <- assessment_rows(
            data, 'leave_one_event_out', fold
        )
        analysis <- data[!assessment_index, , drop = FALSE]
        assessment <- data[assessment_index, , drop = FALSE]
        prepared <- prepare_fold_predictors(analysis, assessment, predictors)
        analysis <- prepared$analysis
        assessment <- prepared$assessment
        cache_path <- file.path(
            model_dir,
            paste0('binomial_', programme_key, '_event_', fold, '.rds')
        )
        cached <- if (file.exists(cache_path)) readRDS(cache_path) else NULL
        cache_is_current <- !is.null(cached) &&
            identical(cached$predictors, predictors) &&
            identical(cached$pseudo_trials, 100L) &&
            diagnostics_pass(cached$diagnostics)

        if (cache_is_current) {
            fit <- cached$fit
            predictions <- cached$predictions
            diagnostics <- cached$diagnostics
            cat('Loaded binomial sensitivity:', programme_key, fold, '\n')
        } else {
            cat('Fitting binomial sensitivity:', programme_key, fold, '\n')
            flush.console()
            if (is.null(template_fit)) {
                fit <- brm(
                    formula, data = analysis, family = binomial('logit'),
                    prior = priors, backend = 'rstan', chains = chains,
                    cores = min(chains, getOption('mc.cores', 1L)),
                    iter = iterations, warmup = warmup,
                    seed = 202408L + as.integer(fold),
                    control = list(
                        adapt_delta = adapt_delta, max_treedepth = 14
                    ),
                    refresh = 0, silent = 2
                )
            } else {
                fit <- update(
                    template_fit, newdata = analysis, recompile = FALSE,
                    chains = chains,
                    cores = min(chains, getOption('mc.cores', 1L)),
                    iter = iterations, warmup = warmup,
                    seed = 202408L + as.integer(fold),
                    control = list(
                        adapt_delta = adapt_delta, max_treedepth = 14
                    ),
                    refresh = 0, silent = 2
                )
            }
            probability_draws <- posterior_linpred(
                fit,
                newdata = assessment,
                transform = TRUE,
                re_formula = NA,
                allow_new_levels = TRUE
            )
            predictions <- tibble(
                predicted_mortality = colMeans(probability_draws),
                prediction_q05 = apply(
                    probability_draws, 2, quantile, 0.05
                ),
                prediction_q50 = apply(
                    probability_draws, 2, quantile, 0.50
                ),
                prediction_q95 = apply(
                    probability_draws, 2, quantile, 0.95
                )
            )
            diagnostics <- fit_diagnostics(fit)
            saveRDS(
                list(
                    fit = fit,
                    predictors = predictors,
                    pseudo_trials = 100L,
                    preprocessing = prepared$preprocessing,
                    predictions = predictions,
                    diagnostics = diagnostics
                ),
                cache_path
            )
        }
        template_fit <- fit
        all_predictions <- bind_rows(
            all_predictions,
            assessment |>
                transmute(
                    programme_key,
                    source_observation_id,
                    ReefID,
                    ReefName,
                    event_year,
                    region_block,
                    observed_mortality = mortality_prop
                ) |>
                bind_cols(predictions) |>
                mutate(
                    learner = 'brms_binomial_olre',
                    scheme = 'leave_one_event_out',
                    fold = as.character(fold)
                )
        )
        all_diagnostics <- bind_rows(
            all_diagnostics,
            diagnostics |>
                mutate(
                    programme_key,
                    scheme = 'leave_one_event_out',
                    fold = as.character(fold),
                    n_analysis = nrow(analysis),
                    n_assessment = nrow(assessment)
                )
        )
    }
}

if (any(!is.finite(all_predictions$predicted_mortality)) ||
    any(all_predictions$predicted_mortality < 0 |
        all_predictions$predicted_mortality > 1)) {
    stop('Binomial sensitivity produced invalid predictions')
}

write_csv(
    all_predictions,
    file.path(output_dir, 'binomial_blocked_predictions.csv')
)
write_csv(
    all_diagnostics,
    file.path(output_dir, 'binomial_diagnostics.csv')
)
write_csv(
    prediction_metrics(all_predictions),
    file.path(output_dir, 'binomial_metrics.csv')
)
print(prediction_metrics(all_predictions))

failed_diagnostics <- all_diagnostics |>
    filter(
        !is.finite(max_rhat) | max_rhat > 1.01 |
            !is.finite(min_bulk_ess) | min_bulk_ess < 100 |
            !is.finite(min_tail_ess) | min_tail_ess < 100 |
            divergences > 0 | max_treedepth_hits > 0
    )
if (nrow(failed_diagnostics) > 0L) {
    print(failed_diagnostics)
    stop('One or more binomial sensitivity folds failed diagnostics')
}
