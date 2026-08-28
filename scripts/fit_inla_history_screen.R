# Fast Bayesian screen of the revised thermal-history predictors.
#
# This is deliberately a screening model, not a replacement for the formal
# BRMS observation model. Mortality is decomposed into Bernoulli occurrence and
# conditional beta magnitude. The two posterior means are multiplied to obtain
# expected mortality. Existing blocked folds are reused unchanged.

suppressPackageStartupMessages({
    library(dplyr)
    library(INLA)
    library(readr)
})
source('scripts/formal_model_helpers.R')

set.seed(20260826L)
output_dir <- 'output/inla_history_screen'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

programmes <- c('manta', 'ltmp', 'mmp')
candidates <- c('legacy_history', 'paper_history')
schemes <- c('reef_blocked_5fold', 'leave_one_event_out')

# Shared ecological and observation-context terms. Keeping this set compact
# makes the history comparison identifiable in programmes with fewer than 100
# observations.
shared_raw <- c(
    'ann_maxdhw', 'prop_acropora_pre', 'observed_pre_cover',
    'dhw10_load4', 'dhw_novelty10', 'secc3m_p10', 'cloudp_90'
)
legacy_raw <- c('histmDHW6', 'yrsince6')
paper_raw <- c(
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8'
)

make_features <- function(analysis, assessment, candidate) {
    raw <- c(
        shared_raw,
        if (candidate == 'legacy_history') legacy_raw else paper_raw
    )
    rules <- tibble(predictor = raw, median = NA_real_, mean = NA_real_, sd = NA_real_)

    for (i in seq_along(raw)) {
        name <- raw[[i]]
        observed <- analysis[[name]][is.finite(analysis[[name]])]
        if (length(observed) == 0L) stop('No analysis values for ', name)
        rules$median[[i]] <- median(observed)
        analysis[[name]][!is.finite(analysis[[name]])] <- rules$median[[i]]
        assessment[[name]][!is.finite(assessment[[name]])] <- rules$median[[i]]
        rules$mean[[i]] <- mean(analysis[[name]])
        rules$sd[[i]] <- sd(analysis[[name]])
        if (!is.finite(rules$sd[[i]]) || rules$sd[[i]] == 0) rules$sd[[i]] <- 1
        z_name <- paste0(name, '_z')
        analysis[[z_name]] <- (analysis[[name]] - rules$mean[[i]]) / rules$sd[[i]]
        assessment[[z_name]] <- (assessment[[name]] - rules$mean[[i]]) / rules$sd[[i]]
    }

    add_derived <- function(rows) {
        rows |>
            mutate(
                dhw_excess4_z = pmax(ann_maxdhw - 4, 0) /
                    max(sd(pmax(analysis$ann_maxdhw - 4, 0)), 1e-6),
                dhw_excess8_z = pmax(ann_maxdhw - 8, 0) /
                    max(sd(pmax(analysis$ann_maxdhw - 8, 0)), 1e-6),
                dhw4_x_acropora = dhw_excess4_z * prop_acropora_pre_z,
                dhw8_x_acropora = dhw_excess8_z * prop_acropora_pre_z,
                dhw_no_prior_n6 = as.numeric(dhw_no_prior_n6)
            )
    }
    analysis <- add_derived(analysis)
    assessment <- add_derived(assessment)

    if (candidate == 'legacy_history') {
        analysis <- analysis |>
            mutate(
                dhw_x_history_count = ann_maxdhw_z * histmDHW6_z,
                dhw_x_recovery = ann_maxdhw_z * yrsince6_z
            )
        assessment <- assessment |>
            mutate(
                dhw_x_history_count = ann_maxdhw_z * histmDHW6_z,
                dhw_x_recovery = ann_maxdhw_z * yrsince6_z
            )
        history_terms <- c('histmDHW6_z', 'yrsince6_z')
    } else {
        analysis <- analysis |>
            mutate(
                dhw_x_history_count = ann_maxdhw_z * dhw_events_since2016_n6_z,
                dhw_x_recovery = ann_maxdhw_z * dhw_years_since_last_n6_capped8_z
            )
        assessment <- assessment |>
            mutate(
                dhw_x_history_count = ann_maxdhw_z * dhw_events_since2016_n6_z,
                dhw_x_recovery = ann_maxdhw_z * dhw_years_since_last_n6_capped8_z
            )
        history_terms <- c(
            'dhw_events_since2016_n6_z',
            'dhw_years_since_last_n6_capped8_z', 'dhw_no_prior_n6'
        )
    }

    fixed_terms <- c(
        'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z',
        'prop_acropora_pre_z', 'observed_pre_cover_z',
        'dhw10_load4_z', 'dhw_novelty10_z', 'secc3m_p10_z',
        'cloudp_90_z', history_terms,
        'dhw4_x_acropora', 'dhw8_x_acropora',
        'dhw_x_history_count', 'dhw_x_recovery'
    )
    list(analysis = analysis, assessment = assessment, terms = fixed_terms, rules = rules)
}

random_terms <- paste(
    "f(reef_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.7,0.05))))",
    "f(event_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.7,0.05))))",
    "f(region_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.5,0.05))))",
    sep = ' + '
)

fit_components <- function(analysis, assessment, terms, compute_criteria = FALSE) {
    all_levels <- bind_rows(analysis, assessment) |>
        transmute(ReefID, event_year, region_block)
    reef_levels <- unique(all_levels$ReefID)
    event_levels <- sort(unique(all_levels$event_year))
    region_levels <- unique(all_levels$region_block)

    index_rows <- function(rows) {
        rows |>
            mutate(
                reef_index = match(ReefID, reef_levels),
                event_index = match(event_year, event_levels),
                region_index = match(region_block, region_levels)
            )
    }
    analysis <- index_rows(analysis)
    assessment <- index_rows(assessment)

    occurrence_data <- bind_rows(
        analysis |> mutate(response = as.numeric(mortality_prop > 0)),
        assessment |> mutate(response = NA_real_)
    )
    fixed_rhs <- paste(terms, collapse = ' + ')
    occurrence_formula <- as.formula(paste(
        'response ~', fixed_rhs, '+', random_terms
    ))
    occurrence_fit <- inla(
        occurrence_formula, family = 'binomial', Ntrials = 1,
        data = occurrence_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
        ),
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(
            waic = compute_criteria, dic = compute_criteria,
            mlik = compute_criteria
        ),
        verbose = FALSE
    )

    positive <- analysis |> filter(mortality_prop > 0)
    positive_n <- nrow(positive)
    # The beta likelihood is open on (0,1). This adjustment mainly affects the
    # two exact-one MMP observations; zeros are handled by occurrence.
    positive$response <- (
        positive$mortality_prop * (positive_n - 1) + 0.5
    ) / positive_n
    magnitude_data <- bind_rows(
        positive,
        assessment |> mutate(response = NA_real_)
    )
    magnitude_formula <- as.formula(paste(
        'response ~', fixed_rhs, '+', random_terms
    ))
    magnitude_fit <- inla(
        magnitude_formula, family = 'beta', data = magnitude_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
        ),
        control.family = list(
            hyper = list(theta = list(prior = 'loggamma', param = c(2, 0.1)))
        ),
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(
            waic = compute_criteria, dic = compute_criteria,
            mlik = compute_criteria
        ),
        verbose = FALSE
    )

    occurrence_index <- nrow(analysis) + seq_len(nrow(assessment))
    magnitude_index <- positive_n + seq_len(nrow(assessment))
    occurrence <- occurrence_fit$summary.fitted.values[occurrence_index, ]
    magnitude <- magnitude_fit$summary.fitted.values[magnitude_index, ]
    predictions <- assessment |>
        transmute(
            source_observation_id, ReefID, ReefName, event_year,
            observed_mortality = mortality_prop,
            observed_occurrence = as.numeric(mortality_prop > 0),
            predicted_occurrence = occurrence$mean,
            predicted_positive_mortality = magnitude$mean,
            predicted_mortality = occurrence$mean * magnitude$mean
        )
    list(
        predictions = predictions,
        occurrence = occurrence_fit,
        magnitude = magnitude_fit
    )
}

all_predictions <- tibble()
fit_runtime <- tibble()

for (programme in programmes) {
    rows <- load_programme_rows(programme)
    for (candidate in candidates) {
        for (scheme in schemes) {
            folds <- fold_values(rows, scheme)
            for (fold in folds) {
                assessment_index <- assessment_rows(rows, scheme, fold)
                prepared <- make_features(
                    rows[!assessment_index, , drop = FALSE],
                    rows[assessment_index, , drop = FALSE],
                    candidate
                )
                started <- proc.time()[['elapsed']]
                fitted <- fit_components(
                    prepared$analysis, prepared$assessment, prepared$terms
                )
                elapsed <- proc.time()[['elapsed']] - started
                all_predictions <- bind_rows(
                    all_predictions,
                    fitted$predictions |>
                        mutate(
                            programme_key = programme,
                            candidate = candidate,
                            scheme = scheme,
                            fold = as.character(fold)
                        )
                )
                fit_runtime <- bind_rows(
                    fit_runtime,
                    tibble(
                        programme_key = programme, candidate = candidate,
                        scheme = scheme, fold = as.character(fold),
                        analysis_rows = nrow(prepared$analysis),
                        assessment_rows = nrow(prepared$assessment),
                        elapsed_seconds = elapsed
                    )
                )
                message(
                    programme, ' / ', candidate, ' / ', scheme,
                    ' / ', fold, ': ', round(elapsed, 1), ' s'
                )
            }
        }
    }
}

metrics <- all_predictions |>
    group_by(programme_key, candidate, scheme) |>
    summarise(
        rows = n(),
        rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
        mae = mean(abs(observed_mortality - predicted_mortality)),
        predictive_r2 = 1 -
            sum((observed_mortality - predicted_mortality)^2) /
            sum((observed_mortality - mean(observed_mortality))^2),
        bias = mean(predicted_mortality - observed_mortality),
        occurrence_brier = mean(
            (observed_occurrence - predicted_occurrence)^2
        ),
        severe_rows = sum(observed_mortality >= 0.5),
        severe_bias = if_else(
            severe_rows > 0,
            mean(predicted_mortality[observed_mortality >= 0.5] -
                     observed_mortality[observed_mortality >= 0.5]),
            NA_real_
        ),
        .groups = 'drop'
    )

metrics_2024 <- all_predictions |>
    filter(scheme == 'reef_blocked_5fold', event_year == 2024) |>
    group_by(programme_key, candidate) |>
    summarise(
        rows = n(),
        rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
        mae = mean(abs(observed_mortality - predicted_mortality)),
        predictive_r2 = 1 -
            sum((observed_mortality - predicted_mortality)^2) /
            sum((observed_mortality - mean(observed_mortality))^2),
        bias = mean(predicted_mortality - observed_mortality),
        severe_rows = sum(observed_mortality >= 0.5),
        severe_bias = if_else(
            severe_rows > 0,
            mean(predicted_mortality[observed_mortality >= 0.5] -
                     observed_mortality[observed_mortality >= 0.5]),
            NA_real_
        ),
        .groups = 'drop'
    )

# Production fits provide approximate Bayesian information criteria and
# interpretable posterior summaries. Predictive selection remains based on the
# blocked results above.
production_criteria <- tibble()
fixed_effects <- tibble()
for (programme in programmes) {
    rows <- load_programme_rows(programme)
    for (candidate in candidates) {
        prepared <- make_features(rows, rows[0, , drop = FALSE], candidate)
        fitted <- fit_components(
            prepared$analysis, prepared$assessment, prepared$terms,
            compute_criteria = TRUE
        )
        production_criteria <- bind_rows(
            production_criteria,
            tibble(
                programme_key = programme, candidate = candidate,
                occurrence_waic = fitted$occurrence$waic$waic,
                magnitude_waic = fitted$magnitude$waic$waic,
                total_waic = occurrence_waic + magnitude_waic,
                occurrence_dic = fitted$occurrence$dic$dic,
                magnitude_dic = fitted$magnitude$dic$dic,
                total_dic = occurrence_dic + magnitude_dic
            )
        )
        fixed_effects <- bind_rows(
            fixed_effects,
            as.data.frame(fitted$occurrence$summary.fixed) |>
                tibble::rownames_to_column('term') |>
                as_tibble() |>
                transmute(
                    programme_key = programme, candidate = candidate,
                    component = 'occurrence', term,
                    mean, sd, q025 = `0.025quant`, q50 = `0.5quant`,
                    q975 = `0.975quant`
                ),
            as.data.frame(fitted$magnitude$summary.fixed) |>
                tibble::rownames_to_column('term') |>
                as_tibble() |>
                transmute(
                    programme_key = programme, candidate = candidate,
                    component = 'positive magnitude', term,
                    mean, sd, q025 = `0.025quant`, q50 = `0.5quant`,
                    q975 = `0.975quant`
                )
        )
    }
}

write_csv(all_predictions, file.path(output_dir, 'predictions.csv'), na = '')
write_csv(metrics, file.path(output_dir, 'metrics.csv'), na = '')
write_csv(metrics_2024, file.path(output_dir, 'metrics_2024.csv'), na = '')
write_csv(production_criteria, file.path(output_dir, 'production_criteria.csv'), na = '')
write_csv(fixed_effects, file.path(output_dir, 'fixed_effects.csv'), na = '')
write_csv(fit_runtime, file.path(output_dir, 'runtime.csv'), na = '')

print(metrics)
print(metrics_2024)
print(production_criteria)
