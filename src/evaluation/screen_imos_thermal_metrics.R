# Screen IMOS night-SST thermal metrics in cause-restricted two-part BRTs.
#
# The source-control comparison replaces only current NOAA/local-first DHW with
# an NOAA-equivalent reconstruction from IMOS nighttime SST and SSTAARS MMM.
# Mechanism candidates add continuity, acute peak, MHW, percentile-gated heat
# dose, cooling and multiscale ID/IDF fields. Full-data importance is
# descriptive; leave-event-out and
# reef-blocked performance decides whether formal assessment is merited.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
    library(tidyr)
})
source('src/lib/formal_model_helpers.R')

set.seed(20260914L)
output_dir <- 'output/imos_thermal_screen'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

feature_file <- 'data/processed/imos_thermal_metrics_validation.csv'
noaa_idf_file <- 'data/processed/noaa_idf_metrics_validation.csv'
features <- read_csv(feature_file, show_col_types = FALSE) |>
    left_join(
        read_csv(noaa_idf_file, show_col_types = FALSE) |>
            select(
                ReefID, year, noaa_idf_max_return_period_years,
                noaa_idf_critical_duration_days,
                noaa_idf_mean_log_return_period,
                noaa_idf_persistence_logrp_difference
            ),
        by = c('ReefID', 'year'), relationship = 'one-to-one'
    )
required_features <- c(
    'ReefID', 'year', 'sstaars_mmm_c', 'imos_dhw1_84_max',
    'imos_dhw1_84_max_observed', 'imos_dhw1_84_peak_coverage',
    'imos_event_coverage', 'imos_hotspell_mmm1_max_days',
    'imos_hotspot_3d_max', 'imos_extreme_dhd_mmm2',
    'imos_hot_nights_mmm1_count',
    'imos_mhw_max_duration_days', 'imos_mhw_max_intensity_c',
    'imos_mhw_cumulative_intensity_c_days', 'imos_mhw_category2_days',
    'imos_pbd10_c_weeks', 'imos_pbd11_c_weeks', 'imos_pbd12_c_weeks',
    'imos_pbd12_c_weeks_observed', 'imos_pbd12_coverage',
    'imos_day_to_following_night_drop_c', 'imos_night_relief_fraction',
    'imos_hot_day_night_pair_count', 'imos_hot_day_night_pair_coverage',
    'imos_id_hotspot_level_c',
    'imos_id_hotspot_persistence_slope_c_per_log_day',
    'noaa_idf_max_return_period_years',
    'noaa_idf_critical_duration_days',
    'noaa_idf_mean_log_return_period',
    'noaa_idf_persistence_logrp_difference'
)
missing_features <- setdiff(required_features, names(features))
if (length(missing_features) > 0) {
    stop('IMOS thermal table lacks: ', paste(missing_features, collapse = ', '))
}
if (anyDuplicated(features[c('ReefID', 'year')])) {
    stop('IMOS thermal table has duplicate ReefID-year keys')
}
features <- features |>
    select(all_of(required_features))

is_explicit_nonthermal <- function(rows) {
    coalesce(rows$disturbance_has_cots, FALSE) |
        coalesce(rows$disturbance_has_cyclone, FALSE) |
        coalesce(rows$disturbance_has_flood, FALSE)
}

load_screen_rows <- function(programme) {
    rows <- load_programme_rows(programme) |>
        left_join(
            features,
            by = c('ReefID', 'event_year' = 'year'),
            relationship = 'many-to-one'
        ) |>
        filter(!is_explicit_nonthermal(pick(everything())))
    if (nrow(rows) == 0) stop('No cause-restricted rows for ', programme)
    rows
}

history_for_programme <- c(manta = 'paper', ltmp = 'legacy', mmp = 'legacy')
history_variables <- list(
    legacy = c('histmDHW6', 'yrsince6'),
    paper = c(
        'dhw_events_since2016_n6',
        'dhw_years_since_last_n6_capped8', 'dhw_no_prior_n6'
    )
)
ecological_anchors <- c('prop_acropora_pre', 'observed_pre_cover')
existing_modifiers <- c(
    'dhw10_load4', 'dhw_novelty10', 'secc3m_p10', 'cloudp_90'
)
structure_metrics <- c(
    'imos_hotspell_mmm1_max_days',
    'imos_hotspot_3d_max',
    'imos_hot_nights_mmm1_count'
)
mechanism_pool <- c(structure_metrics, 'imos_extreme_dhd_mmm2')
mhw_metrics <- c(
    'imos_mhw_max_duration_days', 'imos_mhw_max_intensity_c'
)
pbd_metric <- 'imos_pbd12_c_weeks'
cooling_metrics <- c(
    'imos_day_to_following_night_drop_c', 'imos_night_relief_fraction'
)
id_profile_metrics <- c(
    'imos_id_hotspot_level_c',
    'imos_id_hotspot_persistence_slope_c_per_log_day'
)
idf_frequency_metrics <- c(
    'noaa_idf_mean_log_return_period',
    'noaa_idf_persistence_logrp_difference'
)
mechanism_pool <- c(
    mechanism_pool, mhw_metrics, pbd_metric,
    'imos_mhw_cumulative_intensity_c_days', 'imos_mhw_category2_days',
    cooling_metrics, id_profile_metrics, idf_frequency_metrics
)
validation_candidates <- c(
    'noaa_dhw', 'imos_dhw', 'imos_dhw_observed', 'imos_pbd12',
    'noaa_plus_pbd12', 'noaa_plus_structure', 'imos_plus_structure',
    'noaa_plus_mhw', 'imos_plus_mhw', 'noaa_plus_cooling',
    'noaa_plus_compact', 'noaa_plus_id_profile',
    'noaa_plus_idf_frequency'
)
full_data_candidates <- c(validation_candidates, 'both_dhw_mechanism_pool')
schemes <- c('leave_one_event_out', 'reef_blocked_5fold')

candidate_predictors <- function(programme, candidate) {
    history <- history_variables[[history_for_programme[[programme]]]]
    context <- c(ecological_anchors, existing_modifiers, history)
    switch(
        candidate,
        noaa_dhw = c('ann_maxdhw', context),
        imos_dhw = c('imos_dhw1_84_max', context),
        imos_dhw_observed = c('imos_dhw1_84_max_observed', context),
        imos_pbd12 = c(pbd_metric, context),
        noaa_plus_pbd12 = c('ann_maxdhw', context, pbd_metric),
        noaa_plus_structure = c('ann_maxdhw', context, structure_metrics),
        imos_plus_structure = c(
            'imos_dhw1_84_max', context, structure_metrics
        ),
        noaa_plus_mhw = c('ann_maxdhw', context, mhw_metrics),
        imos_plus_mhw = c('imos_dhw1_84_max', context, mhw_metrics),
        noaa_plus_cooling = c('ann_maxdhw', context, cooling_metrics),
        noaa_plus_compact = c(
            'ann_maxdhw', context, mhw_metrics, cooling_metrics
        ),
        noaa_plus_id_profile = c(
            'ann_maxdhw', context, id_profile_metrics
        ),
        noaa_plus_idf_frequency = c(
            'ann_maxdhw', context, idf_frequency_metrics
        ),
        both_dhw_mechanism_pool = c(
            'ann_maxdhw', 'imos_dhw1_84_max', context, mechanism_pool
        ),
        stop('Unknown candidate: ', candidate)
    )
}

impute_from_analysis <- function(analysis, assessment, predictors) {
    rules <- tibble(
        predictor = predictors, median = NA_real_,
        analysis_missing = NA_integer_, assessment_missing = NA_integer_
    )
    for (i in seq_along(predictors)) {
        predictor <- predictors[[i]]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        if (length(observed) == 0) {
            stop('No finite analysis values for ', predictor)
        }
        value <- median(observed)
        analysis_missing <- !is.finite(analysis[[predictor]])
        assessment_missing <- !is.finite(assessment[[predictor]])
        analysis[[predictor]][analysis_missing] <- value
        assessment[[predictor]][assessment_missing] <- value
        rules$median[[i]] <- value
        rules$analysis_missing[[i]] <- sum(analysis_missing)
        rules$assessment_missing[[i]] <- sum(assessment_missing)
    }
    list(analysis = analysis, assessment = assessment, rules = rules)
}

brt_settings <- function(rows) {
    list(
        n.trees = 1200L,
        interaction.depth = 2L,
        shrinkage = 0.02,
        n.minobsinnode = min(5L, max(2L, floor(rows / 20L))),
        bag.fraction = 0.7,
        train.fraction = 1,
        keep.data = FALSE,
        verbose = FALSE
    )
}

fit_brt_component <- function(data, response, predictors, distribution) {
    do.call(gbm, c(
        list(
            formula = reformulate(predictors, response = response),
            data = data,
            distribution = distribution
        ),
        brt_settings(nrow(data))
    ))
}

fit_brt <- function(analysis, predictors, include_direct = FALSE) {
    analysis$has_loss <- as.numeric(analysis$mortality_prop > 0)
    positive <- analysis |> filter(has_loss == 1)
    if (nrow(positive) < 8) stop('Too few positive rows for magnitude BRT')
    epsilon <- min(0.01, 0.5 / nrow(positive))
    positive$positive_logit <- qlogis(pmin(
        pmax(positive$mortality_prop, epsilon), 1 - epsilon
    ))
    fitted <- list(
        occurrence = fit_brt_component(
            analysis, 'has_loss', predictors, 'bernoulli'
        ),
        magnitude = fit_brt_component(
            positive, 'positive_logit', predictors, 'gaussian'
        )
    )
    if (include_direct) {
        fitted$direct <- fit_brt_component(
            analysis, 'mortality_prop', predictors, 'gaussian'
        )
    }
    fitted
}

predict_brt <- function(model, assessment) {
    occurrence <- predict(
        model$occurrence, assessment, n.trees = 1200, type = 'response'
    )
    magnitude <- plogis(predict(
        model$magnitude, assessment, n.trees = 1200, type = 'response'
    ))
    tibble(
        predicted_occurrence = pmin(pmax(occurrence, 0), 1),
        predicted_positive_mortality = pmin(pmax(magnitude, 0), 1),
        predicted_mortality = predicted_occurrence *
            predicted_positive_mortality
    )
}

importance_table <- function(model, programme, candidate) {
    bind_rows(lapply(names(model), function(component) {
        summary(model[[component]], plotit = FALSE) |>
            as_tibble() |>
            transmute(
                programme_key = programme,
                candidate,
                component,
                predictor = var,
                relative_influence = rel.inf,
                rank = rank(-rel.inf, ties.method = 'min')
            )
    }))
}

safe_rmse <- function(observed, predicted, selected) {
    if (!any(selected)) return(NA_real_)
    sqrt(mean((observed[selected] - predicted[selected])^2))
}

metric_summary <- function(rows, groups) {
    rows |>
        group_by(across(all_of(groups))) |>
        summarise(
            rows = n(),
            rmse = sqrt(mean(
                (observed_mortality - predicted_mortality)^2
            )),
            mae = mean(abs(observed_mortality - predicted_mortality)),
            predictive_r2 = 1 -
                sum((observed_mortality - predicted_mortality)^2) /
                sum((observed_mortality - mean(observed_mortality))^2),
            bias = mean(predicted_mortality - observed_mortality),
            occurrence_brier = mean(
                (observed_occurrence - predicted_occurrence)^2
            ),
            severe_rows = sum(observed_mortality >= 0.2),
            severe_rmse = safe_rmse(
                observed_mortality, predicted_mortality,
                observed_mortality >= 0.2
            ),
            false_extreme_rate = mean(
                predicted_mortality >= 0.5 & observed_mortality < 0.2
            ),
            .groups = 'drop'
        )
}

# Descriptive full-data fits: direct and two-part component importance.
full_importance <- tibble()
apparent_predictions <- tibble()
for (programme in names(validation_files)) {
    rows <- load_screen_rows(programme)
    for (candidate in full_data_candidates) {
        set.seed(
            20260914L + 100L * match(programme, names(validation_files)) +
                match(candidate, full_data_candidates)
        )
        predictors <- candidate_predictors(programme, candidate)
        prepared <- impute_from_analysis(
            rows, rows[0, , drop = FALSE], predictors
        )
        model <- fit_brt(
            prepared$analysis, predictors, include_direct = TRUE
        )
        full_importance <- bind_rows(
            full_importance,
            importance_table(model, programme, candidate)
        )
        predicted <- predict_brt(model, prepared$analysis)
        direct <- predict(
            model$direct, prepared$analysis,
            n.trees = 1200, type = 'response'
        )
        apparent_predictions <- bind_rows(
            apparent_predictions,
            bind_cols(
                prepared$analysis |>
                    transmute(
                        programme_key = .env$programme,
                        source_observation_id, ReefID, event_year,
                        observed_mortality = mortality_prop,
                        observed_occurrence = as.numeric(mortality_prop > 0)
                    ),
                predicted
            ) |>
                mutate(
                    predicted_direct = pmin(pmax(direct, 0), 1),
                    candidate
                )
        )
        message('Full data: ', programme, ' / ', candidate)
    }
}

# Fixed-setting blocked ablations. Preprocessing is learned inside each fold.
predictions <- tibble()
preprocessing <- tibble()
for (programme in names(validation_files)) {
    rows <- load_screen_rows(programme)
    for (candidate in validation_candidates) {
        predictors <- candidate_predictors(programme, candidate)
        for (scheme in schemes) {
            for (fold in fold_values(rows, scheme)) {
                held_out <- assessment_rows(rows, scheme, fold)
                prepared <- impute_from_analysis(
                    rows[!held_out, , drop = FALSE],
                    rows[held_out, , drop = FALSE],
                    predictors
                )
                set.seed(
                    20260914L +
                        1000L * match(programme, names(validation_files)) +
                        100L * match(scheme, schemes) +
                        match(
                            as.character(fold),
                            as.character(fold_values(rows, scheme))
                        )
                )
                predicted <- predict_brt(
                    fit_brt(prepared$analysis, predictors),
                    prepared$assessment
                )
                predictions <- bind_rows(
                    predictions,
                    bind_cols(
                        prepared$assessment |>
                            transmute(
                                programme_key = .env$programme,
                                source_observation_id, ReefID, event_year,
                                observed_mortality = mortality_prop,
                                observed_occurrence =
                                    as.numeric(mortality_prop > 0),
                                noaa_dhw = ann_maxdhw,
                                imos_dhw = imos_dhw1_84_max
                            ),
                        predicted
                    ) |>
                        mutate(
                            candidate,
                            scheme,
                            fold = as.character(fold)
                        )
                )
                preprocessing <- bind_rows(
                    preprocessing,
                    prepared$rules |>
                        mutate(
                            programme_key = programme,
                            candidate, scheme,
                            fold = as.character(fold)
                        )
                )
                message(
                    programme, ' / ', candidate, ' / ', scheme, ' / ', fold
                )
            }
        }
    }
}

metrics <- metric_summary(
    predictions, c('programme_key', 'candidate', 'scheme')
)
metrics_by_event <- metric_summary(
    predictions,
    c('programme_key', 'candidate', 'scheme', 'event_year')
)

# Mechanism sensitivities on matched rows with observed candidate blocks.
# These condition on satellite availability and are not operational results,
# but separate signal from fold-wise median-imputation dilution.
complete_case_predictions <- tibble()
complete_case_skipped <- tibble()
complete_case_sets <- list(
    structure = list(
        metrics = structure_metrics,
        candidates = c(
            'noaa_dhw', 'imos_dhw',
            'noaa_plus_structure', 'imos_plus_structure'
        )
    ),
    mhw = list(
        metrics = mhw_metrics,
        candidates = c(
            'noaa_dhw', 'imos_dhw', 'noaa_plus_mhw', 'imos_plus_mhw'
        )
    ),
    pbd = list(
        metrics = pbd_metric,
        candidates = c(
            'noaa_dhw', 'imos_dhw', 'imos_pbd12', 'noaa_plus_pbd12'
        )
    ),
    cooling = list(
        metrics = cooling_metrics,
        candidates = c('noaa_dhw', 'noaa_plus_cooling')
    ),
    compact = list(
        metrics = c(mhw_metrics, cooling_metrics),
        candidates = c('noaa_dhw', 'noaa_plus_compact')
    ),
    id_profile = list(
        metrics = id_profile_metrics,
        candidates = c('noaa_dhw', 'noaa_plus_id_profile')
    ),
    idf_frequency = list(
        metrics = idf_frequency_metrics,
        candidates = c('noaa_dhw', 'noaa_plus_idf_frequency')
    )
)
for (comparison_set in names(complete_case_sets)) {
  specification <- complete_case_sets[[comparison_set]]
  for (programme in names(validation_files)) {
    rows <- load_screen_rows(programme) |>
        filter(if_all(all_of(specification$metrics), is.finite))
    for (candidate in specification$candidates) {
        predictors <- candidate_predictors(programme, candidate)
        for (scheme in schemes) {
            supported <- all(vapply(
                fold_values(rows, scheme),
                function(test_fold) {
                    training <- rows[
                        !assessment_rows(rows, scheme, test_fold),
                        , drop = FALSE
                    ]
                    sum(training$mortality_prop > 0) >= 8 &&
                        n_distinct(training$mortality_prop > 0) == 2
                },
                logical(1)
            ))
            if (!supported) {
                complete_case_skipped <- bind_rows(
                    complete_case_skipped,
                    tibble(
                        comparison_set, programme_key = programme,
                        candidate, scheme,
                        reason = paste(
                            'at least one fold has fewer than eight positive',
                            'rows or one occurrence class'
                        )
                    )
                )
                next
            }
            for (fold in fold_values(rows, scheme)) {
                held_out <- assessment_rows(rows, scheme, fold)
                prepared <- impute_from_analysis(
                    rows[!held_out, , drop = FALSE],
                    rows[held_out, , drop = FALSE],
                    predictors
                )
                set.seed(
                    20261914L +
                        10000L * match(
                            comparison_set, names(complete_case_sets)
                        ) +
                        1000L * match(programme, names(validation_files)) +
                        100L * match(scheme, schemes) +
                        match(
                            as.character(fold),
                            as.character(fold_values(rows, scheme))
                        )
                )
                predicted <- predict_brt(
                    fit_brt(prepared$analysis, predictors),
                    prepared$assessment
                )
                complete_case_predictions <- bind_rows(
                    complete_case_predictions,
                    bind_cols(
                        prepared$assessment |>
                            transmute(
                                programme_key = .env$programme,
                                source_observation_id, ReefID, event_year,
                                observed_mortality = mortality_prop,
                                observed_occurrence =
                                    as.numeric(mortality_prop > 0)
                            ),
                        predicted
                    ) |>
                        mutate(
                            comparison_set, candidate, scheme,
                            fold = as.character(fold)
                        )
                )
            }
        }
    }
  }
}
complete_case_metrics <- metric_summary(
    complete_case_predictions,
    c('comparison_set', 'programme_key', 'candidate', 'scheme')
)
apparent_metrics <- apparent_predictions |>
    group_by(programme_key, candidate) |>
    summarise(
        rows = n(),
        two_part_rmse = sqrt(mean(
            (observed_mortality - predicted_mortality)^2
        )),
        direct_rmse = sqrt(mean(
            (observed_mortality - predicted_direct)^2
        )),
        .groups = 'drop'
    )

source_comparison <- metrics |>
    filter(candidate %in% c(
        'noaa_dhw', 'imos_dhw', 'imos_dhw_observed', 'imos_pbd12'
    )) |>
    select(
        programme_key, scheme, candidate,
        rmse, predictive_r2, occurrence_brier, severe_rmse,
        false_extreme_rate
    ) |>
    pivot_wider(
        names_from = candidate,
        values_from = c(
            rmse, predictive_r2, occurrence_brier, severe_rmse,
            false_extreme_rate
        )
    ) |>
    mutate(
        rmse_change_imos_minus_noaa = rmse_imos_dhw - rmse_noaa_dhw,
        r2_change_imos_minus_noaa =
            predictive_r2_imos_dhw - predictive_r2_noaa_dhw,
        brier_change_imos_minus_noaa =
            occurrence_brier_imos_dhw - occurrence_brier_noaa_dhw,
        severe_rmse_change_imos_minus_noaa =
            severe_rmse_imos_dhw - severe_rmse_noaa_dhw,
        rmse_change_observed_imos_minus_noaa =
            rmse_imos_dhw_observed - rmse_noaa_dhw,
        r2_change_observed_imos_minus_noaa =
            predictive_r2_imos_dhw_observed - predictive_r2_noaa_dhw,
        brier_change_observed_imos_minus_noaa =
            occurrence_brier_imos_dhw_observed - occurrence_brier_noaa_dhw,
        severe_rmse_change_observed_imos_minus_noaa =
            severe_rmse_imos_dhw_observed - severe_rmse_noaa_dhw,
        rmse_change_pbd_minus_noaa = rmse_imos_pbd12 - rmse_noaa_dhw,
        r2_change_pbd_minus_noaa =
            predictive_r2_imos_pbd12 - predictive_r2_noaa_dhw,
        brier_change_pbd_minus_noaa =
            occurrence_brier_imos_pbd12 - occurrence_brier_noaa_dhw,
        severe_rmse_change_pbd_minus_noaa =
            severe_rmse_imos_pbd12 - severe_rmse_noaa_dhw
    )

coverage_audit <- features |>
    group_by(year) |>
    summarise(
        reef_events = n(),
        finite_imos_dhw = sum(is.finite(imos_dhw1_84_max)),
        finite_structure = sum(
            is.finite(imos_hotspell_mmm1_max_days) &
                is.finite(imos_hotspot_3d_max) &
                is.finite(imos_hot_nights_mmm1_count)
        ),
        finite_mhw = sum(
            is.finite(imos_mhw_max_duration_days) &
                is.finite(imos_mhw_max_intensity_c)
        ),
        finite_pbd = sum(is.finite(imos_pbd12_c_weeks)),
        finite_cooling = sum(
            is.finite(imos_day_to_following_night_drop_c) &
                is.finite(imos_night_relief_fraction)
        ),
        finite_id_profile = sum(
            is.finite(imos_id_hotspot_level_c) &
                is.finite(imos_id_hotspot_persistence_slope_c_per_log_day)
        ),
        finite_idf_frequency = sum(
            is.finite(noaa_idf_mean_log_return_period) &
                is.finite(noaa_idf_persistence_logrp_difference)
        ),
        median_hot_pair_coverage = median(
            imos_hot_day_night_pair_coverage, na.rm = TRUE
        ),
        median_peak_coverage = median(
            imos_dhw1_84_peak_coverage, na.rm = TRUE
        ),
        minimum_peak_coverage = min(
            imos_dhw1_84_peak_coverage, na.rm = TRUE
        ),
        median_event_coverage = median(imos_event_coverage, na.rm = TRUE),
        minimum_event_coverage = min(imos_event_coverage, na.rm = TRUE),
        .groups = 'drop'
    )

feature_rows <- bind_rows(lapply(names(validation_files), function(programme) {
    load_screen_rows(programme) |>
        mutate(programme_key = programme)
}))
correlation_variables <- c(
    'ann_maxdhw', 'imos_dhw1_84_max', mechanism_pool
)
correlations <- feature_rows |>
    select(all_of(correlation_variables)) |>
    cor(use = 'pairwise.complete.obs', method = 'spearman') |>
    as.data.frame() |>
    tibble::rownames_to_column('feature_1') |>
    pivot_longer(
        -feature_1, names_to = 'feature_2', values_to = 'spearman_rho'
    )

write_csv(full_importance, file.path(output_dir, 'full_data_importance.csv'))
write_csv(apparent_metrics, file.path(output_dir, 'apparent_metrics.csv'))
write_csv(predictions, file.path(output_dir, 'blocked_predictions.csv'))
write_csv(metrics, file.path(output_dir, 'blocked_metrics.csv'))
write_csv(metrics_by_event, file.path(output_dir, 'blocked_metrics_by_event.csv'))
write_csv(
    complete_case_predictions,
    file.path(output_dir, 'complete_case_predictions.csv')
)
write_csv(
    complete_case_metrics,
    file.path(output_dir, 'complete_case_metrics.csv')
)
write_csv(
    complete_case_skipped,
    file.path(output_dir, 'complete_case_skipped.csv')
)
write_csv(source_comparison, file.path(output_dir, 'source_comparison.csv'))
write_csv(preprocessing, file.path(output_dir, 'fold_preprocessing.csv'))
write_csv(coverage_audit, file.path(output_dir, 'coverage_audit.csv'))
write_csv(correlations, file.path(output_dir, 'correlations.csv'))
write_csv(
    tibble(
        seed = 20260914L,
        trees = 1200L,
        interaction_depth = 2L,
        shrinkage = 0.02,
        bag_fraction = 0.7,
        cause_restriction =
            'exclude explicit COTS, cyclone and flood observations',
        selection_rule =
            'leave-one-event-out primary; reef-blocked supporting',
        pbd_specification =
            'SSTAARS p90-gated full anomaly; 12-week primary window'
    ),
    file.path(output_dir, 'screen_manifest.csv')
)

print(coverage_audit)
print(source_comparison)
