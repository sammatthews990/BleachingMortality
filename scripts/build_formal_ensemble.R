# Combine matched BRMS, BRT, and benchmark outer-fold predictions without using
# an assessment row's outcome to estimate its ensemble weight. We retain a
# BRMS--BRT-only ensemble for framework comparison and a three-member stack for
# prediction. Production weights use all event- and reef-blocked predictions;
# region-blocked results remain extrapolation stress tests.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})
source("scripts/formal_model_helpers.R")

output_dir <- "output/formal_models"
required <- file.path(
    output_dir,
    c("brms_blocked_predictions.csv", "brt_blocked_predictions.csv")
)
selected_file <- "output/model_validation/selected_models.csv"
benchmark_file <- "output/model_validation/blocked_predictions.csv"
if (any(!file.exists(c(required, selected_file, benchmark_file)))) {
    stop(
        "Base-learner predictions are incomplete. Run the benchmark, BRMS, ",
        "and BRT scripts before building the ensemble."
    )
}

brms_predictions <- read_csv(required[[1]], show_col_types = FALSE)
brt_predictions <- read_csv(required[[2]], show_col_types = FALSE)
selected <- read_csv(selected_file, show_col_types = FALSE) |>
    select(programme_key, model)
benchmark_predictions <- read_csv(benchmark_file, show_col_types = FALSE) |>
    inner_join(selected, by = c("programme_key", "model"))

keys <- c(
    "programme_key", "source_observation_id", "ReefID", "ReefName",
    "event_year", "region_block", "observed_mortality", "scheme", "fold"
)
matched <- inner_join(
    brms_predictions |>
        select(
            all_of(keys), brms_prediction = predicted_mortality,
            brms_occurrence = predicted_occurrence,
            brms_positive = predicted_positive_loss
        ),
    brt_predictions |>
        select(
            all_of(keys), brt_prediction = predicted_mortality,
            brt_occurrence = predicted_occurrence,
            brt_positive = predicted_positive_loss
        ),
    by = keys
) |>
    inner_join(
        benchmark_predictions |>
            select(
                all_of(keys), benchmark_prediction = predicted_mortality,
                benchmark_occurrence = predicted_occurrence,
                benchmark_positive = predicted_positive_loss
            ),
        by = keys
)
if (nrow(matched) != nrow(brms_predictions) ||
    nrow(matched) != nrow(brt_predictions) ||
    nrow(matched) != nrow(benchmark_predictions)) {
    stop("Benchmark, BRMS, and BRT predictions do not match one-to-one")
}

convex_weight <- function(observed, brms_prediction, brt_prediction) {
    difference <- brms_prediction - brt_prediction
    denominator <- sum(difference^2)
    if (!is.finite(denominator) || denominator < 1e-12) return(0.5)
    weight <- sum(difference * (observed - brt_prediction)) / denominator
    min(max(weight, 0), 1)
}

simplex_weights <- function(observed, brms_prediction, brt_prediction,
                            benchmark_prediction) {
    objective <- function(weight) {
        weight_benchmark <- 1 - sum(weight)
        prediction <- weight[[1]] * brms_prediction +
            weight[[2]] * brt_prediction +
            weight_benchmark * benchmark_prediction
        mean((observed - prediction)^2)
    }
    fit <- constrOptim(
        theta = c(1 / 3, 1 / 3),
        f = objective,
        grad = NULL,
        ui = rbind(c(1, 0), c(0, 1), c(-1, -1)),
        ci = c(0, 0, -1)
    )
    weights <- c(
        weight_brms = fit$par[[1]],
        weight_brt = fit$par[[2]],
        weight_benchmark = 1 - sum(fit$par)
    )
    weights[weights < 1e-7] <- 0
    weights / sum(weights)
}

# Cross-fit the meta-learner: each fold's weight is learned exclusively from
# the other folds under the same programme and validation scheme.
ensemble_predictions <- tibble()
fold_weights <- tibble()
groups <- split(matched, interaction(
    matched$programme_key, matched$scheme, drop = TRUE
))
for (group in groups) {
    for (held_out_fold in unique(group$fold)) {
        meta_analysis <- group[group$fold != held_out_fold, , drop = FALSE]
        assessment <- group[group$fold == held_out_fold, , drop = FALSE]
        weight_brms <- convex_weight(
            meta_analysis$observed_mortality,
            meta_analysis$brms_prediction,
            meta_analysis$brt_prediction
        )
        formal_predictions <- weight_brms * assessment$brms_prediction +
            (1 - weight_brms) * assessment$brt_prediction
        weights <- simplex_weights(
            meta_analysis$observed_mortality,
            meta_analysis$brms_prediction,
            meta_analysis$brt_prediction,
            meta_analysis$benchmark_prediction
        )
        stacked_predictions <- weights[["weight_brms"]] *
            assessment$brms_prediction +
            weights[["weight_brt"]] * assessment$brt_prediction +
            weights[["weight_benchmark"]] * assessment$benchmark_prediction

        formal_occurrence <- weight_brms * assessment$brms_occurrence +
            (1 - weight_brms) * assessment$brt_occurrence
        formal_positive <- weight_brms * assessment$brms_positive +
            (1 - weight_brms) * assessment$brt_positive
        stacked_occurrence <- weights[['weight_brms']] *
            assessment$brms_occurrence +
            weights[['weight_brt']] * assessment$brt_occurrence +
            weights[['weight_benchmark']] * assessment$benchmark_occurrence
        stacked_positive <- weights[['weight_brms']] *
            assessment$brms_positive +
            weights[['weight_brt']] * assessment$brt_positive +
            weights[['weight_benchmark']] * assessment$benchmark_positive

        ensemble_predictions <- bind_rows(
            ensemble_predictions,
            assessment |>
                transmute(
                    across(all_of(keys)),
                    learner = "formal_ensemble",
                    predicted_occurrence = formal_occurrence,
                    predicted_positive_loss = formal_positive,
                    predicted_mortality = formal_predictions,
                    weight_brms = weight_brms,
                    weight_brt = 1 - weight_brms,
                    weight_benchmark = 0
                ),
            assessment |>
                transmute(
                    across(all_of(keys)),
                    learner = "stacked_ensemble",
                    predicted_occurrence = stacked_occurrence,
                    predicted_positive_loss = stacked_positive,
                    predicted_mortality = stacked_predictions,
                    weight_brms = weights[["weight_brms"]],
                    weight_brt = weights[["weight_brt"]],
                    weight_benchmark = weights[["weight_benchmark"]]
                )
        )
        fold_weights <- bind_rows(
            fold_weights,
            tibble(
                programme_key = first(group$programme_key),
                scheme = first(group$scheme),
                fold = as.character(held_out_fold),
                learner = "formal_ensemble",
                meta_analysis_rows = nrow(meta_analysis),
                weight_brms = weight_brms,
                weight_brt = 1 - weight_brms,
                weight_benchmark = 0
            ),
            tibble(
                programme_key = first(group$programme_key),
                scheme = first(group$scheme),
                fold = as.character(held_out_fold),
                learner = "stacked_ensemble",
                meta_analysis_rows = nrow(meta_analysis),
                weight_brms = weights[["weight_brms"]],
                weight_brt = weights[["weight_brt"]],
                weight_benchmark = weights[["weight_benchmark"]]
            )
        )
    }
}

# Production weights use only the two prespecified selection schemes. Each row
# contributes once per scheme; region-blocked stress tests do not tune them.
formal_production_weights <- matched |>
    filter(scheme %in% selection_schemes) |>
    group_by(programme_key) |>
    group_modify(~ {
        weight_brms <- convex_weight(
            .x$observed_mortality,
            .x$brms_prediction,
            .x$brt_prediction
        )
        tibble(
            n_cross_fitted_predictions = nrow(.x),
            learner = "formal_ensemble",
            weight_brms = weight_brms,
            weight_brt = 1 - weight_brms,
            weight_benchmark = 0
        )
    }) |>
    ungroup()

stacked_production_weights <- matched |>
    filter(scheme %in% selection_schemes) |>
    group_by(programme_key) |>
    group_modify(~ {
        weights <- simplex_weights(
            .x$observed_mortality,
            .x$brms_prediction,
            .x$brt_prediction,
            .x$benchmark_prediction
        )
        tibble(
            n_cross_fitted_predictions = nrow(.x),
            learner = "stacked_ensemble",
            weight_brms = weights[["weight_brms"]],
            weight_brt = weights[["weight_brt"]],
            weight_benchmark = weights[["weight_benchmark"]]
        )
    }) |>
    ungroup()

production_weights <- bind_rows(
    formal_production_weights,
    stacked_production_weights
)

base_predictions <- bind_rows(
    brms_predictions |>
        select(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, region_block, observed_mortality, scheme, fold,
            learner, predicted_occurrence, predicted_positive_loss,
            predicted_mortality
        ),
    brt_predictions |>
        select(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, region_block, observed_mortality, scheme, fold,
            learner, predicted_occurrence, predicted_positive_loss,
            predicted_mortality
        ),
    benchmark_predictions |>
        transmute(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, region_block, observed_mortality, scheme, fold,
            learner = "benchmark",
            predicted_occurrence, predicted_positive_loss,
            predicted_mortality
        )
)

comparison_predictions <- bind_rows(
    base_predictions,
    ensemble_predictions |>
        select(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, region_block, observed_mortality, scheme, fold,
            learner, predicted_occurrence, predicted_positive_loss,
            predicted_mortality
        )
)
comparison_metrics <- prediction_metrics(comparison_predictions)

write_csv(ensemble_predictions, file.path(output_dir, "ensemble_blocked_predictions.csv"))
write_csv(fold_weights, file.path(output_dir, "ensemble_fold_weights.csv"))
write_csv(production_weights, file.path(output_dir, "production_ensemble_weights.csv"))
write_csv(comparison_predictions, file.path(output_dir, "model_comparison_predictions.csv"))
write_csv(comparison_metrics, file.path(output_dir, "model_comparison_metrics.csv"))

event_component_calibration <- comparison_predictions |>
    filter(
        learner %in% c('brms', 'brt', 'formal_ensemble'),
        scheme %in% c('leave_one_event_out', 'reef_blocked_5fold'),
        event_year %in% c(2022L, 2024L)
    ) |>
    group_by(programme_key, learner, scheme, event_year) |>
    summarise(
        n = n(),
        obs = mean(observed_mortality),
        pred = mean(predicted_mortality),
        obs_occ = mean(observed_mortality > 0),
        pred_occ = mean(predicted_occurrence),
        obs_pos = if_else(
            any(observed_mortality > 0),
            mean(observed_mortality[observed_mortality > 0]),
            NA_real_
        ),
        pred_pos = if_else(
            any(observed_mortality > 0),
            mean(predicted_positive_loss[observed_mortality > 0]),
            NA_real_
        ),
        .groups = 'drop'
    )
write_csv(
    event_component_calibration,
    file.path(output_dir, 'event_component_calibration.csv')
)

event_context <- bind_rows(lapply(names(validation_files), function(key) {
    readRDS(validation_files[[key]])
})) |>
    group_by(programme_key, event_year, region_block) |>
    summarise(
        n = n(),
        mortality = mean(mortality_prop),
        occurrence = mean(mortality_prop > 0),
        dhw = mean(ann_maxdhw, na.rm = TRUE),
        load10 = mean(dhw10_load4, na.rm = TRUE),
        novelty10 = mean(dhw_novelty10, na.rm = TRUE),
        secchi_mean = mean(secc3m, na.rm = TRUE),
        secchi_p10 = mean(secc3m_p10, na.rm = TRUE),
        .groups = 'drop'
    )
write_csv(
    event_context,
    file.path(output_dir, 'event_context_predictor_audit.csv')
)

print(production_weights)
print(comparison_metrics)
