# Fit the formal hierarchical Bayesian mortality learner under the fixed outer
# validation folds, then refit it on every modern event including 2024.
#
# Cross-validation fits use two chains by default to keep repeated assessment
# feasible. The settings can be tightened with the FORMAL_BRMS_CV_* environment
# variables when an individual cached fold needs a targeted diagnostic refit.
# Production fits always use four chains and longer sampling. Every fit is
# cached separately, so an interrupted run resumes without discarding work.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(posterior)
    library(readr)
})
source("scripts/formal_model_helpers.R")

set.seed(202408L)
rstan::rstan_options(auto_write = TRUE)
options(mc.cores = min(4L, parallel::detectCores(logical = FALSE)))

output_dir <- "output/formal_models"
model_dir <- "output/models/formal"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

cv_chains <- as.integer(Sys.getenv("FORMAL_BRMS_CV_CHAINS", "2"))
cv_iter <- as.integer(Sys.getenv("FORMAL_BRMS_CV_ITER", "2000"))
cv_warmup <- as.integer(Sys.getenv("FORMAL_BRMS_CV_WARMUP", "1000"))
cv_adapt_delta <- as.numeric(Sys.getenv("FORMAL_BRMS_CV_ADAPT_DELTA", "0.99"))

refit_failed_diagnostics <- identical(
    Sys.getenv('FORMAL_BRMS_REFIT_DIAGNOSTICS'), '1'
)
production_chains <- as.integer(
    Sys.getenv('FORMAL_BRMS_PRODUCTION_CHAINS', '4')
)
production_iter <- as.integer(
    Sys.getenv('FORMAL_BRMS_PRODUCTION_ITER', '2400')
)
production_warmup <- as.integer(
    Sys.getenv('FORMAL_BRMS_PRODUCTION_WARMUP', '1200')
)
production_adapt_delta <- as.numeric(
    Sys.getenv('FORMAL_BRMS_PRODUCTION_ADAPT_DELTA', '0.99')
)

diagnostics_pass <- function(diagnostics) {
    diagnostics$divergences == 0 &&
        diagnostics$max_treedepth_hits == 0 &&
        diagnostics$max_rhat <= 1.01 &&
        diagnostics$min_bulk_ess >= 100 &&
        diagnostics$min_tail_ess >= 100
}

make_model_formula <- function(programme_key, predictors) {
    mean_formula <- as.formula(
        paste("mortality_prop ~", brms_rhs(predictors))
    )
    boundary_terms <- ~ 1 + ann_maxdhw_z + cloudp_90_z + secc3m_z +
        histmDHW6_z + prop_acropora_pre_z + acropora_cover_pre_z +
        dhw10_load4_z + dhw_novelty10_z + secc3m_p10_z +
        (1 | event_effect) + (1 | region_effect)

    if (programme_key == "mmp") {
        bf(mean_formula, zoi = boundary_terms, coi = ~ 1)
    } else {
        bf(mean_formula, zi = boundary_terms)
    }
}

make_model_family <- function(programme_key) {
    if (programme_key == "mmp") {
        zero_one_inflated_beta()
    } else {
        zero_inflated_beta()
    }
}

make_model_priors <- function(programme_key) {
    priors <- c(
        prior(normal(0, 0.5), class = "b"),
        prior(normal(-1.5, 1), class = "Intercept"),
        prior(exponential(1), class = "sd"),
        prior(gamma(2, 0.1), class = "phi")
    )
    if (programme_key == "mmp") {
        c(
            priors,
            prior(normal(0, 0.7), class = "b", dpar = "zoi"),
            prior(normal(-0.5, 1), class = "Intercept", dpar = "zoi"),
            prior(normal(-2.7, 1), class = "Intercept", dpar = "coi")
        )
    } else {
        c(
            priors,
            prior(normal(0, 0.7), class = "b", dpar = "zi"),
            prior(normal(0, 1), class = "Intercept", dpar = "zi")
        )
    }
}

fit_or_update <- function(template, formula, family, priors, data, seed,
                          chains, iter, warmup, adapt_delta = 0.99) {
    arguments <- list(
        formula = formula,
        data = data,
        family = family,
        prior = priors,
        backend = "rstan",
        chains = chains,
        cores = min(chains, getOption("mc.cores", 1L)),
        iter = iter,
        warmup = warmup,
        seed = seed,
        control = list(adapt_delta = adapt_delta, max_treedepth = 12),
        refresh = 0,
        silent = 2
    )
    if (is.null(template)) {
        return(do.call(brm, arguments))
    }

    update(
        template,
        newdata = data,
        recompile = FALSE,
        chains = chains,
        cores = min(chains, getOption("mc.cores", 1L)),
        iter = iter,
        warmup = warmup,
        seed = seed,
        control = list(adapt_delta = adapt_delta, max_treedepth = 12),
        refresh = 0,
        silent = 2
    )
}

summarise_fit_diagnostics <- function(fit) {
    draw_summary <- posterior::summarise_draws(
        posterior::as_draws_array(fit),
        rhat = posterior::rhat,
        ess_bulk = posterior::ess_bulk,
        ess_tail = posterior::ess_tail
    )
    nuts <- nuts_params(fit)
    tibble(
        max_rhat = max(draw_summary$rhat, na.rm = TRUE),
        min_bulk_ess = min(draw_summary$ess_bulk, na.rm = TRUE),
        min_tail_ess = min(draw_summary$ess_tail, na.rm = TRUE),
        divergences = sum(
            nuts$Value[nuts$Parameter == "divergent__"],
            na.rm = TRUE
        ),
        max_treedepth_hits = sum(
            nuts$Value[nuts$Parameter == "treedepth__"] >= 12,
            na.rm = TRUE
        )
    )
}

random_effects_for_scheme <- function(scheme) {
    switch(
        scheme,
        leave_one_event_out = NA,
        reef_blocked_5fold = ~ (1 | event_effect) + (1 | region_effect),
        region_blocked = ~ (1 | event_effect),
        production_known_event = ~ (1 | event_effect) + (1 | region_effect),
        stop('Unknown BRMS prediction mode: ', scheme)
    )
}

predict_response <- function(fit, assessment, programme_key, scheme) {
    included_random_effects <- random_effects_for_scheme(scheme)
    draws <- posterior_epred(
        fit,
        newdata = assessment,
        re_formula = included_random_effects,
        allow_new_levels = TRUE
    )
    if (programme_key == 'mmp') {
        zoi <- posterior_linpred(
            fit, newdata = assessment, dpar = 'zoi', transform = TRUE,
            re_formula = included_random_effects, allow_new_levels = TRUE
        )
        coi <- posterior_linpred(
            fit, newdata = assessment, dpar = 'coi', transform = TRUE,
            re_formula = included_random_effects, allow_new_levels = TRUE
        )
        occurrence_draws <- 1 - zoi * (1 - coi)
    } else {
        zi <- posterior_linpred(
            fit, newdata = assessment, dpar = 'zi', transform = TRUE,
            re_formula = included_random_effects, allow_new_levels = TRUE
        )
        occurrence_draws <- 1 - zi
    }
    positive_draws <- draws / pmax(occurrence_draws, 1e-8)
    tibble(
        predicted_occurrence = colMeans(occurrence_draws),
        predicted_positive_loss = pmin(
            pmax(colMeans(positive_draws), 0), 1
        ),
        predicted_mortality = colMeans(draws),
        prediction_q05 = apply(draws, 2, quantile, 0.05),
        prediction_q50 = apply(draws, 2, quantile, 0.50),
        prediction_q95 = apply(draws, 2, quantile, 0.95)
    )
}

dry_run <- identical(Sys.getenv("FORMAL_BRMS_DRY_RUN"), "1")
if (dry_run) {
    data <- load_programme_rows("mmp")
    predictors <- programme_predictors("mmp")
    prepared <- prepare_fold_predictors(data, data, predictors)
    formula <- make_model_formula("mmp", predictors)
    family <- make_model_family("mmp")
    priors <- make_model_priors("mmp")
    get_prior(formula, prepared$analysis, family = family, prior = priors)
    make_stancode(formula, prepared$analysis, family = family, prior = priors)
    cat("BRMS dry run passed: formula, family, priors, and Stan code are valid.\n")
    quit(save = "no", status = 0)
}

all_predictions <- tibble()
all_diagnostics <- tibble()

programme_keys <- names(validation_files)
requested_programme <- trimws(Sys.getenv('FORMAL_PROGRAMME_KEY', ''))
if (nzchar(requested_programme)) {
    if (!requested_programme %in% programme_keys) {
        stop('Unknown FORMAL_PROGRAMME_KEY: ', requested_programme)
    }
    programme_keys <- requested_programme
}

for (programme_key in programme_keys) {
    data <- load_programme_rows(programme_key)
    predictors <- programme_predictors(programme_key)
    formula <- make_model_formula(programme_key, predictors)
    family <- make_model_family(programme_key)
    priors <- make_model_priors(programme_key)
    template_fit <- NULL

    for (scheme in all_validation_schemes) {
        for (fold in fold_values(data, scheme)) {
            cache_path <- file.path(
                model_dir,
                paste0(
                    "brms_", programme_key, "_", scheme, "_",
                    safe_fold_label(fold), ".rds"
                )
            )
            assessment_index <- assessment_rows(data, scheme, fold)
            analysis <- data[!assessment_index, , drop = FALSE]
            assessment <- data[assessment_index, , drop = FALSE]
            prepared <- prepare_fold_predictors(analysis, assessment, predictors)
            analysis <- prepared$analysis
            assessment <- prepared$assessment

            cached <- if (file.exists(cache_path)) readRDS(cache_path) else NULL
            cache_is_current <- !is.null(cached) &&
                identical(cached$predictors, predictors) &&
                all(c(
                    'predicted_occurrence', 'predicted_positive_loss'
                ) %in% names(cached$predictions))
            if (cache_is_current && refit_failed_diagnostics) {
                cache_is_current <- diagnostics_pass(cached$diagnostics)
            }
            if (cache_is_current) {
                fit <- cached$fit
                predictions <- cached$predictions
                diagnostics <- cached$diagnostics
                cat("Loaded cached BRMS fold:", programme_key, scheme, fold, "\n")
            } else {
                cat("Fitting BRMS fold:", programme_key, scheme, fold,
                    "with", nrow(analysis), "analysis rows\n")
                flush.console()
                fit <- fit_or_update(
                    template_fit, formula, family, priors, analysis,
                    seed = 202408L + nrow(all_diagnostics) + 1L,
                    chains = cv_chains,
                    iter = cv_iter,
                    warmup = cv_warmup,
                    adapt_delta = cv_adapt_delta
                )
                predictions <- predict_response(
                    fit, assessment, programme_key, scheme
                )
                diagnostics <- summarise_fit_diagnostics(fit)
                saveRDS(
                    list(
                        fit = fit,
                        preprocessing = prepared$preprocessing,
                        predictors = predictors,
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
                        learner = "brms",
                        scheme = scheme,
                        fold = as.character(fold)
                    )
            )
            all_diagnostics <- bind_rows(
                all_diagnostics,
                diagnostics |>
                    mutate(
                        programme_key = programme_key,
                        scheme = scheme,
                        fold = as.character(fold),
                        n_analysis = nrow(analysis),
                        n_assessment = nrow(assessment)
                    )
            )
        }
    }

    full_prepared <- prepare_fold_predictors(data, data, predictors)
    full_cache <- file.path(model_dir, paste0("brms_", programme_key, "_production.rds"))
    production_cached <- if (file.exists(full_cache)) {
        readRDS(full_cache)
    } else NULL
    production_cache_is_current <- !is.null(production_cached) &&
        identical(production_cached$predictors, predictors)
    if (production_cache_is_current && refit_failed_diagnostics) {
        production_cache_is_current <- diagnostics_pass(
            production_cached$diagnostics
        )
    }
    if (production_cache_is_current) {
        cat("Loaded cached production BRMS model:", programme_key, "\n")
    } else {
        cat("Fitting production BRMS model:", programme_key,
            "on all", nrow(data), "rows including 2024\n")
        flush.console()
        production_fit <- fit_or_update(
            template_fit, formula, family, priors, full_prepared$analysis,
            seed = 202408L,
            chains = production_chains,
            iter = production_iter,
            warmup = production_warmup,
            adapt_delta = production_adapt_delta
        )
        saveRDS(
            list(
                fit = production_fit,
                preprocessing = full_prepared$preprocessing,
                predictors = predictors,
                diagnostics = summarise_fit_diagnostics(production_fit),
                event_years = sort(unique(data$event_year)),
                n = nrow(data)
            ),
            full_cache
        )
    }
}

if (anyDuplicated(paste(
    all_predictions$programme_key,
    all_predictions$source_observation_id,
    all_predictions$scheme
))) {
    stop("A BRMS observation appears in multiple assessment folds")
}
if (any(!is.finite(all_predictions$predicted_mortality)) ||
    any(all_predictions$predicted_mortality < 0 | all_predictions$predicted_mortality > 1)) {
    stop("BRMS produced invalid mortality predictions")
}

write_csv(all_predictions, file.path(output_dir, "brms_blocked_predictions.csv"))
write_csv(all_diagnostics, file.path(output_dir, "brms_diagnostics.csv"))
write_csv(
    prediction_metrics(all_predictions),
    file.path(output_dir, "brms_metrics.csv")
)

print(prediction_metrics(all_predictions))
