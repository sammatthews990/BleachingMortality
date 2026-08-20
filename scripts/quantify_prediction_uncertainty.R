# Calibrate predictive intervals from genuinely held-out residuals and screen
# the 2024 GBR environmental grid against each programme's training domain.
# This script does not make a GBR mortality map: required pre-event community
# covariates are not yet available for every registry feature.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})
source('scripts/formal_model_helpers.R')

formal_dir <- 'output/formal_models'
output_dir <- 'output/prediction_uncertainty'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

comparison_file <- file.path(formal_dir, 'model_comparison_predictions.csv')
brms_file <- file.path(formal_dir, 'brms_blocked_predictions.csv')
environment_file <- 'data/processed/cheung_recreated_gbr_full.csv'
required_files <- c(comparison_file, brms_file, environment_file)
if (any(!file.exists(required_files))) {
    stop('Formal predictions or the full environmental table are missing.')
}

comparison <- read_csv(comparison_file, show_col_types = FALSE)
brms_predictions <- read_csv(brms_file, show_col_types = FALSE)

# These routes are the prespecified choices documented in the formal report.
production_routes <- tribble(
    ~programme_key, ~learner,
    'ltmp', 'stacked_ensemble',
    'manta', 'formal_ensemble',
    'mmp', 'formal_ensemble'
)
calibration_schemes <- c('leave_one_event_out', 'reef_blocked_5fold')

selected <- comparison |>
    inner_join(production_routes, by = c('programme_key', 'learner')) |>
    filter(scheme %in% calibration_schemes)

prediction_keys <- c(
    'programme_key', 'source_observation_id', 'scheme', 'fold'
)
base_predictions <- comparison |>
    filter(
        learner %in% c('brms', 'brt'),
        scheme %in% calibration_schemes
    ) |>
    select(all_of(prediction_keys), learner, predicted_mortality) |>
    pivot_wider(
        names_from = learner,
        values_from = predicted_mortality,
        names_prefix = 'prediction_'
    )

selected <- selected |>
    left_join(base_predictions, by = prediction_keys) |>
    mutate(
        model_disagreement = abs(prediction_brms - prediction_brt),
        residual = observed_mortality - predicted_mortality
    )

# For a 90% split-conformal interval, use the finite-sample corrected empirical
# quantile of absolute residuals. Calibration rows always come from other outer
# folds, so an assessment outcome never calibrates its own interval.
conformal_quantile <- function(absolute_residuals, coverage = 0.90) {
    n <- length(absolute_residuals)
    probability <- min(1, ceiling((n + 1) * coverage) / n)
    as.numeric(quantile(
        absolute_residuals, probs = probability, type = 1, names = FALSE
    ))
}

interval_rows <- list()
row_index <- 1L
for (current_programme in unique(selected$programme_key)) {
    for (current_scheme in calibration_schemes) {
        programme_scheme <- selected |>
            filter(
                .data$programme_key == current_programme,
                .data$scheme == current_scheme
            )
        for (held_out_fold in unique(programme_scheme$fold)) {
            assessment <- filter(programme_scheme, fold == held_out_fold)
            calibration <- filter(programme_scheme, fold != held_out_fold)
            radius <- conformal_quantile(abs(calibration$residual), 0.90)
            severe_calibration <- calibration |>
                filter(observed_mortality >= 0.2)
            severe_radius <- conformal_quantile(
                abs(severe_calibration$residual), 0.90
            )
            interval_rows[[row_index]] <- assessment |>
                mutate(
                    interval_lower_90 = pmax(predicted_mortality - radius, 0),
                    interval_upper_90 = pmin(predicted_mortality + radius, 1),
                    interval_radius_90 = radius,
                    calibration_n = nrow(calibration),
                    severe_upper_guard = pmin(
                        predicted_mortality + severe_radius, 1
                    ),
                    severe_guard_radius = severe_radius,
                    severe_calibration_n = nrow(severe_calibration),
                    interval_covers = observed_mortality >= interval_lower_90 &
                        observed_mortality <= interval_upper_90,
                    severe_upper_guard_covers =
                        observed_mortality <= severe_upper_guard
                )
            row_index <- row_index + 1L
        }
    }
}
interval_predictions <- bind_rows(interval_rows)

coverage_rows <- list()
row_index <- 1L
subgroups <- list(
    all = rep(TRUE, nrow(interval_predictions)),
    zero = interval_predictions$observed_mortality == 0,
    positive = interval_predictions$observed_mortality > 0,
    severe = interval_predictions$observed_mortality >= 0.2
)
for (subgroup_name in names(subgroups)) {
    subgroup <- interval_predictions[subgroups[[subgroup_name]], , drop = FALSE]
    coverage_rows[[row_index]] <- subgroup |>
        group_by(programme_key, scheme, learner) |>
        summarise(
            subgroup = subgroup_name,
            n = n(),
            empirical_coverage = mean(interval_covers),
            mean_interval_width = mean(interval_upper_90 - interval_lower_90),
            mean_interval_radius = mean(interval_radius_90),
            severe_upper_guard_coverage = mean(severe_upper_guard_covers),
            mean_severe_upper_guard = mean(severe_upper_guard),
            mean_model_disagreement = mean(model_disagreement),
            .groups = 'drop'
        )
    row_index <- row_index + 1L
}
coverage_summary <- bind_rows(coverage_rows)

# The saved BRMS limits are posterior intervals for expected mortality because
# they were generated with posterior_epred. Their empirical coverage is useful,
# but they are not observation-level predictive intervals.
brms_expected_interval_coverage <- brms_predictions |>
    filter(scheme %in% calibration_schemes) |>
    group_by(programme_key, scheme) |>
    summarise(
        n = n(),
        expected_interval_coverage = mean(
            observed_mortality >= prediction_q05 &
                observed_mortality <= prediction_q95
        ),
        severe_n = sum(observed_mortality >= 0.2),
        severe_expected_interval_coverage = mean(
            observed_mortality[observed_mortality >= 0.2] >=
                prediction_q05[observed_mortality >= 0.2] &
            observed_mortality[observed_mortality >= 0.2] <=
                prediction_q95[observed_mortality >= 0.2]
        ),
        mean_expected_interval_width = mean(prediction_q95 - prediction_q05),
        .groups = 'drop'
    )

write_csv(
    interval_predictions,
    file.path(output_dir, 'selected_route_cross_calibrated_intervals.csv')
)
write_csv(
    coverage_summary,
    file.path(output_dir, 'selected_route_interval_coverage.csv')
)
write_csv(
    brms_expected_interval_coverage,
    file.path(output_dir, 'brms_expected_interval_coverage.csv')
)

# Environmental support for every 2024 GBR registry feature. Community
# covariates are intentionally excluded here and reported as a separate gate.
registry_2024 <- read_csv(environment_file, show_col_types = FALSE) |>
    filter(year == 2024)
if (nrow(registry_2024) != 7063L) {
    stop('Expected 7,063 registry features for 2024, found ', nrow(registry_2024))
}

environment_only <- setdiff(
    environmental_predictors,
    c('prop_acropora_pre', 'acropora_cover_pre')
)
support_rows <- list()
threshold_rows <- list()
predictor_rows <- list()

for (programme_key in names(validation_files)) {
    training <- load_programme_rows(programme_key)
    prepared <- prepare_fold_predictors(
        training, registry_2024, environment_only
    )
    scaled_names <- paste0(environment_only, '_z')
    training_scaled <- as.matrix(prepared$analysis[scaled_names])
    registry_scaled <- as.matrix(prepared$assessment[scaled_names])

    training_distances <- as.matrix(dist(training_scaled))
    diag(training_distances) <- Inf
    training_nearest <- apply(training_distances, 1, min)
    threshold_95 <- as.numeric(quantile(training_nearest, 0.95, type = 8))
    threshold_99 <- as.numeric(quantile(training_nearest, 0.99, type = 8))

    nearest_distance <- numeric(nrow(registry_scaled))
    block_starts <- seq(1L, nrow(registry_scaled), by = 250L)
    for (block_start in block_starts) {
        block_end <- min(block_start + 249L, nrow(registry_scaled))
        block <- registry_scaled[block_start:block_end, , drop = FALSE]
        squared_distance <- matrix(
            0, nrow = nrow(block), ncol = nrow(training_scaled)
        )
        for (predictor_index in seq_along(environment_only)) {
            squared_distance <- squared_distance + outer(
                block[, predictor_index],
                training_scaled[, predictor_index],
                '-'
            )^2
        }
        nearest_distance[block_start:block_end] <- sqrt(
            apply(squared_distance, 1, min)
        )
    }

    outside_matrix <- matrix(
        FALSE,
        nrow = nrow(registry_2024),
        ncol = length(environment_only),
        dimnames = list(NULL, environment_only)
    )
    missing_matrix <- outside_matrix
    for (predictor in environment_only) {
        training_min <- min(training[[predictor]], na.rm = TRUE)
        training_max <- max(training[[predictor]], na.rm = TRUE)
        missing_matrix[, predictor] <- !is.finite(registry_2024[[predictor]])
        outside_matrix[, predictor] <-
            registry_2024[[predictor]] < training_min |
            registry_2024[[predictor]] > training_max
        outside_matrix[missing_matrix[, predictor], predictor] <- FALSE

        predictor_rows[[length(predictor_rows) + 1L]] <- tibble(
            programme_key = programme_key,
            predictor = predictor,
            training_min = training_min,
            training_max = training_max,
            registry_min = min(registry_2024[[predictor]], na.rm = TRUE),
            registry_max = max(registry_2024[[predictor]], na.rm = TRUE),
            missing_features = sum(missing_matrix[, predictor]),
            below_training_range = sum(
                registry_2024[[predictor]] < training_min, na.rm = TRUE
            ),
            above_training_range = sum(
                registry_2024[[predictor]] > training_max, na.rm = TRUE
            )
        )
    }

    outside_count <- rowSums(outside_matrix)
    missing_count <- rowSums(missing_matrix)
    outside_predictors <- apply(outside_matrix, 1, function(flag) {
        paste(environment_only[flag], collapse = ';')
    })
    missing_predictors <- apply(missing_matrix, 1, function(flag) {
        paste(environment_only[flag], collapse = ';')
    })

    support_class <- case_when(
        missing_count > 0 ~ 'missing environmental input',
        outside_count >= 2 | nearest_distance > threshold_99 ~
            'extreme extrapolation',
        outside_count == 1 | nearest_distance > threshold_95 ~
            'caution',
        TRUE ~ 'within environmental support'
    )
    missing_required <- if (programme_key == 'mmp') {
        'prop_acropora_pre;acropora_cover_pre;depth'
    } else {
        'prop_acropora_pre;acropora_cover_pre'
    }

    support_rows[[programme_key]] <- registry_2024 |>
        transmute(
            programme_key = programme_key,
            LABEL_ID, LOC_NAME_S, lon, lat, year,
            environmental_nearest_distance = nearest_distance,
            environmental_distance_q95 = threshold_95,
            environmental_distance_q99 = threshold_99,
            environmental_outside_range_count = outside_count,
            environmental_outside_range_predictors = outside_predictors,
            environmental_missing_count = missing_count,
            environmental_missing_predictors = missing_predictors,
            environmental_support_class = support_class,
            missing_required_model_predictors = missing_required,
            full_model_prediction_ready = FALSE
        )
    threshold_rows[[programme_key]] <- tibble(
        programme_key = programme_key,
        training_rows = nrow(training),
        environmental_predictors = length(environment_only),
        nearest_distance_q95 = threshold_95,
        nearest_distance_q99 = threshold_99
    )
}

environmental_support <- bind_rows(support_rows)
environmental_thresholds <- bind_rows(threshold_rows)
environmental_predictor_summary <- bind_rows(predictor_rows)

write_csv(
    environmental_support,
    file.path(output_dir, 'environmental_support_2024.csv')
)
write_csv(
    environmental_thresholds,
    file.path(output_dir, 'environmental_support_thresholds.csv')
)
write_csv(
    environmental_predictor_summary,
    file.path(output_dir, 'environmental_support_predictor_summary.csv')
)

cat('\nCross-calibrated interval coverage:\n')
print(coverage_summary)
cat('\n2024 environmental support:\n')
print(environmental_support |>
    count(programme_key, environmental_support_class))
cat('\nGBR mortality map gate: required community covariates are not complete.\n')
