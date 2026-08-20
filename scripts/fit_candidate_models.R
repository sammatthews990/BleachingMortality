# Fit and compare a small, prespecified set of coral-mortality models.
#
# Exact zero and one outcomes are retained. Each candidate is a two-part model:
# a binomial model for whether cover loss occurs and a fractional-logit model
# for loss magnitude conditional on positive loss. All preprocessing is learned
# from the analysis portion of each blocked fold.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(splines)
})

set.seed(202408L)

validation_files <- c(
    manta = "data/processed/validation_rows_manta.rds",
    ltmp = "data/processed/validation_rows_ltmp.rds",
    mmp = "data/processed/validation_rows_mmp.rds"
)
missing_inputs <- validation_files[!file.exists(validation_files)]
if (length(missing_inputs) > 0L) {
    stop("Missing validation tables: ", paste(missing_inputs, collapse = ", "),
         ". Run scripts/build_validation_splits.R first.")
}

output_dir <- "output/model_validation"
model_dir <- "output/models"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

# These are nested scientific hypotheses, not a broad automated model search.
# Complexity is used by the parsimony rule after blocked validation.
candidate_specs <- tibble(
    model = c(
        'thermal_composition',
        "thermal_linear",
        "thermal_nonlinear",
        "thermal_optical",
        "thermal_context"
    ),
    complexity = c(5L, 1:4),
    hypothesis = c(
        paste(
            'All attainable thermal, optical, physical and pre-event',
            'Acropora context improves prediction'
        ),
        "Mortality changes linearly with event DHW",
        "The DHW response is smooth and nonlinear",
        "Cloud and water clarity modify the nonlinear DHW response",
        "Current exposure and heat-stress history add transferable information"
    ),
    terms = c(
        paste(
            'splines::ns(ann_maxdhw_z, df = 3)',
            '+ cloudp_90_z + secc3m_z + mcur_90_z + histmDHW6_z',
            '+ prop_acropora_pre_z + acropora_cover_pre_z',
            '+ ann_maxdhw_z:prop_acropora_pre_z'
        ),
        "ann_maxdhw_z",
        "splines::ns(ann_maxdhw_z, df = 3)",
        paste(
            "splines::ns(ann_maxdhw_z, df = 3)",
            "+ cloudp_90_z + secc3m_z"
        ),
        paste(
            "splines::ns(ann_maxdhw_z, df = 3)",
            "+ cloudp_90_z + secc3m_z + mcur_90_z + histmDHW6_z"
        )
    )
)

predictors <- c("ann_maxdhw", "cloudp_90", "secc3m", "mcur_90", "histmDHW6")
schemes <- c("leave_one_event_out", "reef_blocked_5fold", "region_blocked")

predictors <- c('prop_acropora_pre', 'acropora_cover_pre', predictors)

prepare_predictors <- function(analysis, assessment, predictor_names) {
    preprocessing <- tibble(
        predictor = predictor_names,
        median = NA_real_,
        mean = NA_real_,
        sd = NA_real_
    )

    for (i in seq_along(predictor_names)) {
        predictor <- predictor_names[i]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        if (length(observed) == 0L) {
            stop("No finite training values for ", predictor)
        }
        impute_value <- median(observed)
        analysis[[predictor]][!is.finite(analysis[[predictor]])] <- impute_value
        assessment[[predictor]][!is.finite(assessment[[predictor]])] <- impute_value

        centre <- mean(analysis[[predictor]])
        scale <- sd(analysis[[predictor]])
        if (!is.finite(scale) || scale == 0) scale <- 1

        analysis[[paste0(predictor, "_z")]] <-
            (analysis[[predictor]] - centre) / scale
        assessment[[paste0(predictor, "_z")]] <-
            (assessment[[predictor]] - centre) / scale

        preprocessing$median[i] <- impute_value
        preprocessing$mean[i] <- centre
        preprocessing$sd[i] <- scale
    }

    list(analysis = analysis, assessment = assessment, preprocessing = preprocessing)
}

fit_two_part_model <- function(data, terms) {
    data$has_loss <- data$mortality_prop > 0
    occurrence_formula <- as.formula(paste("has_loss ~", terms))
    magnitude_formula <- as.formula(paste("mortality_prop ~", terms))

    occurrence <- glm(occurrence_formula, data = data, family = binomial())
    positive_rows <- data[data$has_loss, , drop = FALSE]
    magnitude <- glm(
        magnitude_formula,
        data = positive_rows,
        family = quasibinomial(link = "logit")
    )

    list(occurrence = occurrence, magnitude = magnitude)
}

predict_two_part_model <- function(model, new_data) {
    occurrence <- predict(model$occurrence, newdata = new_data, type = "response")
    magnitude <- predict(model$magnitude, newdata = new_data, type = "response")
    expected <- pmin(pmax(occurrence * magnitude, 0), 1)
    tibble(
        predicted_occurrence = occurrence,
        predicted_positive_loss = magnitude,
        predicted_mortality = expected
    )
}

all_predictions <- tibble()
warning_log <- tibble()

for (programme_name in names(validation_files)) {
    data <- readRDS(validation_files[[programme_name]]) |>
        arrange(event_year, ReefID, source_observation_id)

    for (scheme in schemes) {
        fold_values <- switch(
            scheme,
            leave_one_event_out = sort(unique(data$event_year)),
            reef_blocked_5fold = sort(unique(data$reef_fold)),
            region_blocked = sort(unique(data$region_block))
        )

        for (fold in fold_values) {
            assessment_index <- switch(
                scheme,
                leave_one_event_out = data$event_year == as.integer(fold),
                reef_blocked_5fold = data$reef_fold == as.integer(fold),
                region_blocked = data$region_block == as.character(fold)
            )
            analysis <- data[!assessment_index, , drop = FALSE]
            assessment <- data[assessment_index, , drop = FALSE]

            prepared <- prepare_predictors(analysis, assessment, predictors)
            analysis <- prepared$analysis
            assessment <- prepared$assessment

            for (candidate in seq_len(nrow(candidate_specs))) {
                specification <- candidate_specs[candidate, ]
                fit_warnings <- character()
                result <- withCallingHandlers(
                    {
                        fitted <- fit_two_part_model(analysis, specification$terms)
                        predictions <- predict_two_part_model(fitted, assessment)
                        list(fitted = fitted, predictions = predictions)
                    },
                    warning = function(warning) {
                        fit_warnings <<- c(fit_warnings, conditionMessage(warning))
                        invokeRestart("muffleWarning")
                    }
                )
                fitted <- result$fitted
                predictions <- result$predictions

                if (length(fit_warnings) > 0L) {
                    warning_log <- bind_rows(warning_log, tibble(
                        programme_key = programme_name,
                        scheme = scheme,
                        fold = as.character(fold),
                        model = specification$model,
                        warning = unique(fit_warnings)
                    ))
                }

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
                            ann_maxdhw,
                            observed_mortality = mortality_prop,
                            observed_loss = mortality_prop > 0
                        ) |>
                        bind_cols(predictions) |>
                        mutate(
                            scheme = scheme,
                            fold = as.character(fold),
                            model = specification$model,
                            complexity = specification$complexity
                        )
                )
            }
        }
    }
}

if (any(!is.finite(all_predictions$predicted_mortality)) ||
    any(all_predictions$predicted_mortality < 0 | all_predictions$predicted_mortality > 1)) {
    stop("Blocked validation produced invalid mortality predictions")
}

# Every observation must appear exactly once per programme/model/scheme.
prediction_keys <- paste(
    all_predictions$programme_key,
    all_predictions$source_observation_id,
    all_predictions$model,
    all_predictions$scheme
)
if (anyDuplicated(prediction_keys)) {
    stop("An observation appears in multiple assessment folds for the same scheme")
}

calibration <- function(observed, predicted) {
    clipped <- pmin(pmax(predicted, 1e-6), 1 - 1e-6)
    fit <- tryCatch(
        suppressWarnings(glm(
            observed ~ qlogis(clipped),
            family = quasibinomial(link = "logit")
        )),
        error = function(e) NULL
    )
    if (is.null(fit) || length(coef(fit)) < 2L) {
        return(c(intercept = NA_real_, slope = NA_real_))
    }
    c(intercept = unname(coef(fit)[1]), slope = unname(coef(fit)[2]))
}

metric_groups <- split(
    all_predictions,
    interaction(
        all_predictions$programme_key,
        all_predictions$model,
        all_predictions$scheme,
        drop = TRUE
    )
)

metrics <- bind_rows(lapply(metric_groups, function(data) {
    cal <- calibration(data$observed_mortality, data$predicted_mortality)
    tibble(
        programme_key = first(data$programme_key),
        model = first(data$model),
        complexity = first(data$complexity),
        scheme = first(data$scheme),
        n = nrow(data),
        mae = mean(abs(data$observed_mortality - data$predicted_mortality)),
        rmse = sqrt(mean((data$observed_mortality - data$predicted_mortality)^2)),
        predictive_r2 = 1 -
            sum((data$observed_mortality - data$predicted_mortality)^2) /
            sum((data$observed_mortality - mean(data$observed_mortality))^2),
        bias = mean(data$predicted_mortality - data$observed_mortality),
        observed_mean = mean(data$observed_mortality),
        predicted_mean = mean(data$predicted_mortality),
        occurrence_brier = mean(
            (data$observed_loss - data$predicted_occurrence)^2
        ),
        occurrence_log_loss = -mean(
            data$observed_loss * log(pmax(data$predicted_occurrence, 1e-8)) +
                (1 - data$observed_loss) *
                log(pmax(1 - data$predicted_occurrence, 1e-8))
        ),
        zero_mae = mean(abs(
            data$observed_mortality[data$observed_mortality == 0] -
                data$predicted_mortality[data$observed_mortality == 0]
        )),
        positive_mae = mean(abs(
            data$observed_mortality[data$observed_mortality > 0] -
                data$predicted_mortality[data$observed_mortality > 0]
        )),
        balanced_mae = mean(c(zero_mae, positive_mae)),
        severe_mae = mean(abs(
            data$observed_mortality[data$observed_mortality >= 0.2] -
                data$predicted_mortality[data$observed_mortality >= 0.2]
        )),
        calibration_intercept = cal["intercept"],
        calibration_slope = cal["slope"]
    )
}))

fold_metric_groups <- split(
    all_predictions,
    interaction(
        all_predictions$programme_key,
        all_predictions$model,
        all_predictions$scheme,
        all_predictions$fold,
        drop = TRUE
    )
)
fold_metrics <- bind_rows(lapply(fold_metric_groups, function(data) {
    tibble(
        programme_key = first(data$programme_key),
        model = first(data$model),
        complexity = first(data$complexity),
        scheme = first(data$scheme),
        fold = first(data$fold),
        n = nrow(data),
        mae = mean(abs(data$observed_mortality - data$predicted_mortality)),
        rmse = sqrt(mean((data$observed_mortality - data$predicted_mortality)^2)),
        predictive_r2 = 1 -
            sum((data$observed_mortality - data$predicted_mortality)^2) /
            sum((data$observed_mortality - mean(data$observed_mortality))^2),
        bias = mean(data$predicted_mortality - data$observed_mortality),
        observed_mean = mean(data$observed_mortality),
        predicted_mean = mean(data$predicted_mortality),
        occurrence_brier = mean(
            (data$observed_loss - data$predicted_occurrence)^2
        ),
        positive_mae = mean(abs(
            data$observed_mortality[data$observed_mortality > 0] -
                data$predicted_mortality[data$observed_mortality > 0]
        )),
        severe_mae = mean(abs(
            data$observed_mortality[data$observed_mortality >= 0.2] -
                data$predicted_mortality[data$observed_mortality >= 0.2]
        ))
    )
}))

# Use event and reef transfer for model selection. Region blocking remains a
# deliberately difficult stress test. Select the simplest candidate whose mean
# RMSE is within 2% of the best candidate for that programme.
selection_scores <- metrics |>
    filter(scheme %in% c("leave_one_event_out", "reef_blocked_5fold")) |>
    group_by(programme_key, model, complexity) |>
    summarise(
        selection_rmse = mean(rmse),
        selection_mae = mean(mae),
        selection_predictive_r2 = mean(predictive_r2),
        selection_balanced_mae = mean(balanced_mae),
        selection_severe_mae = mean(severe_mae),
        .groups = "drop"
    ) |>
    group_by(programme_key) |>
    mutate(
        best_rmse = min(selection_rmse),
        within_two_percent = selection_rmse <= best_rmse * 1.02
    ) |>
    arrange(programme_key, desc(within_two_percent), complexity, selection_rmse) |>
    group_by(programme_key) |>
    mutate(selected_for_parsimony = row_number() == 1L) |>
    ungroup() |>
    group_by(programme_key) |>
    arrange(selection_rmse, selection_mae, complexity, .by_group = TRUE) |>
    mutate(selected_for_prediction = row_number() == 1L) |>
    ungroup()

selected_models <- selection_scores |>
    filter(selected_for_prediction) |>
    select(
        programme_key, model, complexity, selection_rmse, selection_mae,
        selection_predictive_r2, selection_balanced_mae, selection_severe_mae
    )

# Refit the selected benchmark on every modern event, including 2024. These
# remain comparison baselines; the BRMS--BRT ensemble is the production target.
benchmark_models <- list()
coefficient_summary <- tibble()
for (programme_name in names(validation_files)) {
    data <- readRDS(validation_files[[programme_name]]) |>
        arrange(event_year, ReefID, source_observation_id)
    selected_name <- selected_models$model[
        selected_models$programme_key == programme_name
    ]
    specification <- candidate_specs[candidate_specs$model == selected_name, ]

    prepared <- prepare_predictors(data, data, predictors)
    fitted <- fit_two_part_model(prepared$analysis, specification$terms)
    benchmark_models[[programme_name]] <- list(
        programme_key = programme_name,
        model = selected_name,
        hypothesis = specification$hypothesis,
        terms = specification$terms,
        predictors = predictors,
        preprocessing = prepared$preprocessing,
        occurrence = fitted$occurrence,
        magnitude = fitted$magnitude,
        event_years = sort(unique(data$event_year)),
        n = nrow(data)
    )
    saveRDS(
        benchmark_models[[programme_name]],
        file.path(model_dir, paste0("mortality_benchmark_two_part_", programme_name, ".rds"))
    )

    for (component in c("occurrence", "magnitude")) {
        coefficient_matrix <- coef(summary(fitted[[component]]))
        coefficient_summary <- bind_rows(coefficient_summary, tibble(
            programme_key = programme_name,
            model = selected_name,
            component = component,
            term = rownames(coefficient_matrix),
            estimate = coefficient_matrix[, 1],
            std_error = coefficient_matrix[, 2],
            statistic = coefficient_matrix[, 3],
            p_value = coefficient_matrix[, 4]
        ))
    }
}

write_csv(all_predictions, file.path(output_dir, "blocked_predictions.csv"))
write_csv(metrics, file.path(output_dir, "candidate_metrics.csv"))
write_csv(fold_metrics, file.path(output_dir, "candidate_fold_metrics.csv"))
write_csv(selection_scores, file.path(output_dir, "candidate_selection.csv"))
write_csv(selected_models, file.path(output_dir, "selected_models.csv"))
write_csv(warning_log, file.path(output_dir, "fit_warnings.csv"))
write_csv(coefficient_summary, file.path(output_dir, "production_coefficients.csv"))
saveRDS(benchmark_models, file.path(model_dir, "mortality_two_part_benchmarks.rds"))

# Apparent fit statistics support ecological interpretation. AIC/AICc are
# reported for the Bernoulli occurrence component only; the positive-magnitude
# model is quasibinomial and therefore has no likelihood-based AIC.
explanatory_metrics <- tibble()
for (programme_name in names(validation_files)) {
    data <- readRDS(validation_files[[programme_name]]) |>
        arrange(event_year, ReefID, source_observation_id)
    prepared <- prepare_predictors(data, data, predictors)
    analysis <- prepared$analysis
    analysis$has_loss <- analysis$mortality_prop > 0

    for (candidate in seq_len(nrow(candidate_specs))) {
        specification <- candidate_specs[candidate, ]
        fitted <- suppressWarnings(
            fit_two_part_model(analysis, specification$terms)
        )
        predictions <- predict_two_part_model(fitted, analysis)
        occurrence_null <- glm(
            has_loss ~ 1, data = analysis, family = binomial()
        )
        parameter_count <- attr(logLik(fitted$occurrence), 'df')
        observation_count <- nrow(analysis)
        occurrence_aic <- AIC(fitted$occurrence)

        explanatory_metrics <- bind_rows(
            explanatory_metrics,
            tibble(
                programme_key = programme_name,
                model = specification$model,
                complexity = specification$complexity,
                n = observation_count,
                occurrence_aic = occurrence_aic,
                occurrence_aicc = if (
                    observation_count > parameter_count + 1
                ) occurrence_aic +
                    2 * parameter_count * (parameter_count + 1) /
                    (observation_count - parameter_count - 1) else NA_real_,
                occurrence_mcfadden_r2 = 1 -
                    as.numeric(logLik(fitted$occurrence)) /
                    as.numeric(logLik(occurrence_null)),
                positive_deviance_explained = 1 -
                    fitted$magnitude$deviance /
                    fitted$magnitude$null.deviance,
                apparent_response_r2 = 1 -
                    sum((analysis$mortality_prop -
                         predictions$predicted_mortality)^2) /
                    sum((analysis$mortality_prop -
                         mean(analysis$mortality_prop))^2),
                magnitude_aic = NA_real_,
                magnitude_aic_note =
                    'Undefined for the quasibinomial magnitude component'
            )
        )
    }
}
write_csv(
    explanatory_metrics,
    file.path(output_dir, 'candidate_explanatory_metrics.csv')
)

print(
    left_join(selected_models, candidate_specs, by = c("model", "complexity")) |>
        select(programme_key, model, selection_rmse, selection_mae, hypothesis)
)
