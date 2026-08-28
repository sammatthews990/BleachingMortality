# Sensitivity test for the hypothesis that high observed mortality at very low
# DHW flattens the beta-family mortality response.
#
# Important: suspect rows are removed from each analysis (training) fold only.
# Every held-out observation remains in the assessment fold, including the rows
# labelled as suspect. This prevents apparent gains caused by deleting difficult
# observations from the validation target.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(INLA)
    library(readr)
    library(tidyr)
})
source('scripts/formal_model_helpers.R')

set.seed(20260826L)
output_dir <- 'output/low_dhw_sensitivity'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

programmes <- c('manta', 'ltmp')
variants <- tibble(
    variant = c('all_rows', 'exclude_low_dhw_ge05',
                'exclude_low_dhw_ge10', 'exclude_low_dhw_ge20'),
    loss_threshold = c(Inf, 0.05, 0.10, 0.20),
    label = c('All training rows', 'Exclude DHW < 4 and loss >= 5%',
              'Exclude DHW < 4 and loss >= 10%',
              'Exclude DHW < 4 and loss >= 20%')
)

# Keep the model deliberately aligned with the fast INLA thermal-history screen:
# Bernoulli mortality occurrence multiplied by conditional beta magnitude.
raw_predictors <- c(
    'ann_maxdhw', 'prop_acropora_pre', 'observed_pre_cover',
    'dhw10_load4', 'dhw_novelty10', 'secc3m_p10', 'cloudp_90',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8'
)

prepare_features <- function(analysis, assessment) {
    for (name in raw_predictors) {
        observed <- analysis[[name]][is.finite(analysis[[name]])]
        if (length(observed) == 0L) stop('No analysis values for ', name)
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

    excess4_scale <- sd(pmax(analysis$ann_maxdhw - 4, 0))
    excess8_scale <- sd(pmax(analysis$ann_maxdhw - 8, 0))
    if (!is.finite(excess4_scale) || excess4_scale == 0) excess4_scale <- 1
    if (!is.finite(excess8_scale) || excess8_scale == 0) excess8_scale <- 1

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

fit_two_part <- function(analysis, assessment) {
    all_levels <- bind_rows(analysis, assessment)
    reef_levels <- unique(all_levels$ReefID)
    event_levels <- sort(unique(all_levels$event_year))
    region_levels <- unique(all_levels$region_block)
    add_indices <- function(rows) {
        rows |>
            mutate(
                reef_index = match(ReefID, reef_levels),
                event_index = match(event_year, event_levels),
                region_index = match(region_block, region_levels)
            )
    }
    analysis <- add_indices(analysis)
    assessment <- add_indices(assessment)
    rhs <- paste(c(fixed_terms, random_rhs), collapse = ' + ')

    occurrence_data <- bind_rows(
        analysis |> mutate(response = as.numeric(mortality_prop > 0)),
        assessment |> mutate(response = NA_real_)
    )
    occurrence_fit <- inla(
        as.formula(paste('response ~', rhs)), family = 'binomial',
        Ntrials = 1, data = occurrence_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
        ),
        control.predictor = list(compute = TRUE, link = 1),
        verbose = FALSE
    )

    positive <- analysis |> filter(mortality_prop > 0)
    positive_n <- nrow(positive)
    positive$response <- (
        positive$mortality_prop * (positive_n - 1) + 0.5
    ) / positive_n
    magnitude_data <- bind_rows(
        positive,
        assessment |> mutate(response = NA_real_)
    )
    magnitude_fit <- inla(
        as.formula(paste('response ~', rhs)), family = 'beta',
        data = magnitude_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
        ),
        control.family = list(
            hyper = list(theta = list(prior = 'loggamma', param = c(2, 0.1)))
        ),
        control.predictor = list(compute = TRUE, link = 1),
        verbose = FALSE
    )

    occurrence_index <- nrow(analysis) + seq_len(nrow(assessment))
    magnitude_index <- positive_n + seq_len(nrow(assessment))
    occurrence <- occurrence_fit$summary.fitted.values[occurrence_index, 'mean']
    magnitude <- magnitude_fit$summary.fitted.values[magnitude_index, 'mean']

    assessment |>
        transmute(
            source_observation_id, ReefID, ReefName, event_year,
            ann_maxdhw, mortality_prop,
            suspected_low_dhw_loss = ann_maxdhw < 4 & mortality_prop >= 0.10,
            predicted_occurrence = occurrence,
            predicted_positive_mortality = magnitude,
            predicted_mortality = occurrence * magnitude
        )
}

predictions <- tibble()
runtime <- tibble()
for (programme in programmes) {
    rows <- load_programme_rows(programme)
    for (fold in fold_values(rows, 'reef_blocked_5fold')) {
        held_out <- assessment_rows(rows, 'reef_blocked_5fold', fold)
        assessment <- rows[held_out, , drop = FALSE]
        unfiltered_analysis <- rows[!held_out, , drop = FALSE]

        for (i in seq_len(nrow(variants))) {
            threshold <- variants$loss_threshold[[i]]
            exclude <- if (is.infinite(threshold)) {
                rep(FALSE, nrow(unfiltered_analysis))
            } else {
                unfiltered_analysis$ann_maxdhw < 4 &
                    unfiltered_analysis$mortality_prop >= threshold
            }
            analysis <- unfiltered_analysis[!exclude, , drop = FALSE]
            prepared <- prepare_features(analysis, assessment)
            started <- proc.time()[['elapsed']]
            fitted <- fit_two_part(prepared$analysis, prepared$assessment)
            elapsed <- proc.time()[['elapsed']] - started

            predictions <- bind_rows(
                predictions,
                fitted |>
                    mutate(
                        programme_key = programme,
                        fold = as.character(fold),
                        variant = variants$variant[[i]],
                        variant_label = variants$label[[i]],
                        excluded_training_rows = sum(exclude)
                    )
            )
            runtime <- bind_rows(
                runtime,
                tibble(
                    programme_key = programme, fold = as.character(fold),
                    variant = variants$variant[[i]],
                    analysis_rows = nrow(analysis),
                    excluded_training_rows = sum(exclude),
                    assessment_rows = nrow(assessment),
                    elapsed_seconds = elapsed
                )
            )
            message(
                programme, ' / fold ', fold, ' / ', variants$variant[[i]],
                ': excluded ', sum(exclude), ' rows; ', round(elapsed, 1), ' s'
            )
        }
    }
}

predictions <- predictions |>
    mutate(
        dhw_band = cut(
            ann_maxdhw, breaks = c(-Inf, 4, 8, 12, Inf), right = FALSE,
            labels = c('<4', '4-<8', '8-<12', '>=12')
        ),
        residual = mortality_prop - predicted_mortality
    )

metric_summary <- function(rows) {
    rows |>
        summarise(
            n = n(),
            observed_mean = mean(mortality_prop),
            predicted_mean = mean(predicted_mortality),
            bias_predicted_minus_observed = mean(
                predicted_mortality - mortality_prop
            ),
            mae = mean(abs(predicted_mortality - mortality_prop)),
            rmse = sqrt(mean((predicted_mortality - mortality_prop)^2)),
            predictive_r2 = 1 -
                sum((mortality_prop - predicted_mortality)^2) /
                sum((mortality_prop - mean(mortality_prop))^2),
            .groups = 'drop'
        )
}

overall_metrics <- predictions |>
    group_by(programme_key, variant, variant_label) |>
    metric_summary()
band_metrics <- predictions |>
    group_by(programme_key, variant, variant_label, dhw_band) |>
    metric_summary()
severe_metrics <- predictions |>
    filter(mortality_prop >= 0.50) |>
    group_by(programme_key, variant, variant_label) |>
    metric_summary()
suspect_metrics <- predictions |>
    filter(suspected_low_dhw_loss) |>
    group_by(programme_key, variant, variant_label) |>
    metric_summary()

# Audit the exact rows that motivate the sensitivity analysis. Metadata are
# descriptive only; no disturbance label is used to choose exclusions.
observation_audit <- bind_rows(lapply(c('manta', 'ltmp', 'mmp'), function(p) {
    load_programme_rows(p) |> mutate(programme_key = p)
})) |>
    filter(ann_maxdhw < 4, mortality_prop >= 0.05) |>
    select(
        programme_key, source_observation_id, ReefID, ReefName, event_year,
        ann_maxdhw, mortality_prop, observed_pre_cover, prop_acropora_pre,
        DISTURBANCE_TYPE, disturbance_text, disturbance_has_bleaching,
        disturbance_has_cyclone, disturbance_has_flood, disturbance_has_cots,
        lon, lat
    ) |>
    arrange(desc(mortality_prop))

write_csv(predictions, file.path(output_dir, 'cross_validated_predictions.csv'), na = '')
write_csv(overall_metrics, file.path(output_dir, 'overall_metrics.csv'), na = '')
write_csv(band_metrics, file.path(output_dir, 'dhw_band_metrics.csv'), na = '')
write_csv(severe_metrics, file.path(output_dir, 'severe_mortality_metrics.csv'), na = '')
write_csv(suspect_metrics, file.path(output_dir, 'suspect_row_metrics.csv'), na = '')
write_csv(observation_audit, file.path(output_dir, 'low_dhw_observation_audit.csv'), na = '')
write_csv(runtime, file.path(output_dir, 'runtime.csv'), na = '')

print(overall_metrics)
print(band_metrics)
print(severe_metrics)
