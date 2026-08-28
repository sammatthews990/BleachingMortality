# Test how zero relative-mortality observations affect the two-part model.
#
# Zeros are downsampled from training folds only. All zeros remain in held-out
# assessment rows. Raw downsampled predictions show leverage; prevalence-
# corrected predictions remove the artificial class-balance shift. The
# conditional beta prediction is also retained to show the hypothetical result
# if mortality occurrence were treated as certain.

suppressPackageStartupMessages({
    library(dplyr)
    library(INLA)
    library(readr)
})
source('scripts/formal_model_helpers.R')

set.seed(20260826L)
output_dir <- 'output/zero_mortality_sensitivity'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

programmes <- c('manta', 'ltmp', 'mmp')
zero_retention <- c(1.00, 0.75, 0.50, 0.25)
raw_predictors <- c(
    'ann_maxdhw', 'prop_acropora_pre', 'observed_pre_cover',
    'dhw10_load4', 'dhw_novelty10', 'secc3m_p10', 'cloudp_90',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8'
)

prepare_features <- function(analysis, assessment) {
    for (name in raw_predictors) {
        observed <- analysis[[name]][is.finite(analysis[[name]])]
        fill <- median(observed)
        analysis[[name]][!is.finite(analysis[[name]])] <- fill
        assessment[[name]][!is.finite(assessment[[name]])] <- fill
        centre <- mean(analysis[[name]])
        scale <- sd(analysis[[name]])
        if (!is.finite(scale) || scale == 0) scale <- 1
        zname <- paste0(name, '_z')
        analysis[[zname]] <- (analysis[[name]] - centre) / scale
        assessment[[zname]] <- (assessment[[name]] - centre) / scale
    }
    excess4_scale <- max(sd(pmax(analysis$ann_maxdhw - 4, 0)), 1e-6)
    excess8_scale <- max(sd(pmax(analysis$ann_maxdhw - 8, 0)), 1e-6)
    add_derived <- function(rows) {
        rows |>
            mutate(
                dhw_excess4_z = pmax(ann_maxdhw - 4, 0) / excess4_scale,
                dhw_excess8_z = pmax(ann_maxdhw - 8, 0) / excess8_scale,
                dhw_no_prior_n6 = as.numeric(dhw_no_prior_n6),
                dhw4_x_acropora = dhw_excess4_z * prop_acropora_pre_z,
                dhw8_x_acropora = dhw_excess8_z * prop_acropora_pre_z,
                dhw_x_history_count = ann_maxdhw_z *
                    dhw_events_since2016_n6_z,
                dhw_x_recovery = ann_maxdhw_z *
                    dhw_years_since_last_n6_capped8_z
            )
    }
    list(
        analysis = add_derived(analysis),
        assessment = add_derived(assessment)
    )
}

fixed_terms <- c(
    'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z',
    'prop_acropora_pre_z', 'observed_pre_cover_z',
    'dhw10_load4_z', 'dhw_novelty10_z', 'secc3m_p10_z', 'cloudp_90_z',
    'dhw_events_since2016_n6_z',
    'dhw_years_since_last_n6_capped8_z', 'dhw_no_prior_n6',
    'dhw4_x_acropora', 'dhw8_x_acropora',
    'dhw_x_history_count', 'dhw_x_recovery'
)
random_rhs <- paste(
    "f(reef_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.7,0.05))))",
    "f(event_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.7,0.05))))",
    "f(region_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.5,0.05))))",
    sep = ' + '
)
model_formula <- as.formula(paste(
    'response ~', paste(c(fixed_terms, random_rhs), collapse = ' + ')
))

add_indices <- function(analysis, assessment) {
    all_rows <- bind_rows(analysis, assessment)
    reef_levels <- unique(all_rows$ReefID)
    event_levels <- sort(unique(all_rows$event_year))
    region_levels <- unique(all_rows$region_block)
    index <- function(rows) {
        rows |>
            mutate(
                reef_index = match(ReefID, reef_levels),
                event_index = match(event_year, event_levels),
                region_index = match(region_block, region_levels)
            )
    }
    list(analysis = index(analysis), assessment = index(assessment))
}

fit_positive_beta <- function(analysis, assessment) {
    positive <- analysis |> filter(mortality_prop > 0)
    positive_n <- nrow(positive)
    positive$response <- (
        positive$mortality_prop * (positive_n - 1) + 0.5
    ) / positive_n
    model_data <- bind_rows(
        positive,
        assessment |> mutate(response = NA_real_)
    )
    fit <- inla(
        model_formula, family = 'beta', data = model_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
        ),
        control.family = list(
            hyper = list(theta = list(prior = 'loggamma', param = c(2, 0.1)))
        ),
        control.predictor = list(compute = TRUE, link = 1),
        verbose = FALSE
    )
    fit$summary.fitted.values[
        positive_n + seq_len(nrow(assessment)), 'mean'
    ]
}

retain_training_zeros <- function(rows, retention, seed) {
    positives <- rows |> filter(mortality_prop > 0)
    zeros <- rows |> filter(mortality_prop == 0)
    if (retention >= 1) return(bind_rows(positives, zeros))
    set.seed(seed)
    kept_zeros <- zeros |>
        group_by(event_year) |>
        group_modify(~ slice_sample(
            .x, n = max(1L, round(nrow(.x) * retention))
        )) |>
        ungroup()
    bind_rows(positives, kept_zeros) |>
        arrange(event_year, ReefID, source_observation_id)
}

fit_occurrence <- function(analysis, assessment) {
    model_data <- bind_rows(
        analysis |> mutate(response = as.numeric(mortality_prop > 0)),
        assessment |> mutate(response = NA_real_)
    )
    fit <- inla(
        model_formula, family = 'binomial', Ntrials = 1, data = model_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
        ),
        control.predictor = list(compute = TRUE, link = 1),
        verbose = FALSE
    )
    fit$summary.fitted.values[
        nrow(analysis) + seq_len(nrow(assessment)), 'mean'
    ]
}

predictions <- tibble()
runtime <- tibble()
for (programme in programmes) {
    rows <- load_programme_rows(programme)
    for (fold in fold_values(rows, 'reef_blocked_5fold')) {
        held_out <- assessment_rows(rows, 'reef_blocked_5fold', fold)
        raw_analysis <- rows[!held_out, , drop = FALSE]
        raw_assessment <- rows[held_out, , drop = FALSE]
        prepared <- prepare_features(raw_analysis, raw_assessment)
        indexed <- add_indices(prepared$analysis, prepared$assessment)

        started <- proc.time()[['elapsed']]
        positive_prediction <- fit_positive_beta(
            indexed$analysis, indexed$assessment
        )
        beta_elapsed <- proc.time()[['elapsed']] - started

        for (retention in zero_retention) {
            occurrence_analysis <- retain_training_zeros(
                indexed$analysis, retention,
                seed = 20260826L + fold * 100L +
                    match(programme, programmes) * 1000L +
                    round(retention * 100)
            )
            started <- proc.time()[['elapsed']]
            occurrence_raw <- fit_occurrence(
                occurrence_analysis, indexed$assessment
            )
            occurrence_corrected <- if (retention >= 1) {
                occurrence_raw
            } else {
                # Case-control correction: sampling zeros at rate r inflates
                # the fitted odds of mortality occurrence by approximately 1/r.
                plogis(qlogis(occurrence_raw) + log(retention))
            }
            elapsed <- proc.time()[['elapsed']] - started

            base <- indexed$assessment |>
                transmute(
                    programme_key = .env$programme,
                    fold = as.character(.env$fold),
                    source_observation_id, ReefID, ReefName, event_year,
                    ann_maxdhw, observed_mortality = mortality_prop,
                    zero_retention = retention,
                    training_rows = nrow(occurrence_analysis),
                    training_zeros = sum(occurrence_analysis$mortality_prop == 0),
                    predicted_positive_mortality = positive_prediction
                )
            predictions <- bind_rows(
                predictions,
                base |>
                    mutate(
                        prediction_type = 'raw_downsampled',
                        predicted_occurrence = occurrence_raw,
                        predicted_mortality = occurrence_raw *
                            predicted_positive_mortality
                    ),
                base |>
                    mutate(
                        prediction_type = 'prevalence_corrected',
                        predicted_occurrence = occurrence_corrected,
                        predicted_mortality = occurrence_corrected *
                            predicted_positive_mortality
                    ),
                base |>
                    mutate(
                        prediction_type = 'conditional_positive',
                        predicted_occurrence = 1,
                        predicted_mortality = predicted_positive_mortality
                    )
            )
            runtime <- bind_rows(
                runtime,
                tibble(
                    programme_key = programme, fold = as.character(fold),
                    zero_retention = retention,
                    beta_seconds = beta_elapsed,
                    occurrence_seconds = elapsed
                )
            )
            message(
                programme, ' / fold ', fold, ' / retain ',
                retention * 100, '% zeros: ', round(elapsed, 1), ' s'
            )
        }
    }
}

predictions <- predictions |>
    mutate(
        dhw_band = cut(
            ann_maxdhw, c(-Inf, 4, 8, 12, Inf), right = FALSE,
            labels = c('<4', '4-<8', '8-<12', '>=12')
        ),
        severe_observed = observed_mortality >= 0.50
    )

metric_summary <- function(rows) {
    rows |>
        summarise(
            n = n(),
            observed_mean = mean(observed_mortality),
            predicted_mean = mean(predicted_mortality),
            occurrence_mean = mean(predicted_occurrence),
            positive_mean = mean(predicted_positive_mortality),
            bias = mean(predicted_mortality - observed_mortality),
            mae = mean(abs(predicted_mortality - observed_mortality)),
            rmse = sqrt(mean((predicted_mortality - observed_mortality)^2)),
            predictive_r2 = 1 -
                sum((observed_mortality - predicted_mortality)^2) /
                sum((observed_mortality - mean(observed_mortality))^2),
            .groups = 'drop'
        )
}

overall_metrics <- predictions |>
    group_by(programme_key, zero_retention, prediction_type) |>
    metric_summary()
band_metrics <- predictions |>
    group_by(programme_key, zero_retention, prediction_type, dhw_band) |>
    metric_summary()
severe_metrics <- predictions |>
    filter(severe_observed) |>
    group_by(programme_key, zero_retention, prediction_type) |>
    metric_summary()

zero_audit <- bind_rows(lapply(programmes, function(programme) {
    load_programme_rows(programme) |>
        mutate(programme_key = .env$programme)
})) |>
    mutate(dhw_band = cut(
        ann_maxdhw, c(-Inf, 4, 8, 12, Inf), right = FALSE,
        labels = c('<4', '4-<8', '8-<12', '>=12')
    )) |>
    group_by(programme_key, event_year, dhw_band) |>
    summarise(
        n = n(), zeros = sum(mortality_prop == 0),
        zero_fraction = mean(mortality_prop == 0),
        observed_mean = mean(mortality_prop),
        observed_positive_mean = if_else(
            any(mortality_prop > 0),
            mean(mortality_prop[mortality_prop > 0]), NA_real_
        ),
        .groups = 'drop'
    )

write_csv(predictions, file.path(output_dir, 'cross_validated_predictions.csv'), na = '')
write_csv(overall_metrics, file.path(output_dir, 'overall_metrics.csv'), na = '')
write_csv(band_metrics, file.path(output_dir, 'dhw_band_metrics.csv'), na = '')
write_csv(severe_metrics, file.path(output_dir, 'severe_metrics.csv'), na = '')
write_csv(zero_audit, file.path(output_dir, 'zero_audit.csv'), na = '')
write_csv(runtime, file.path(output_dir, 'runtime.csv'), na = '')

print(overall_metrics)
print(severe_metrics)
