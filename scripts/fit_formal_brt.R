# Fit the formal two-part boosted regression tree learner under the fixed outer
# validation folds, then refit it on all modern events including 2024.
#
# Occurrence is Bernoulli. Positive mortality magnitude compares the original
# Gaussian BRT with a logit-Gaussian BRT that is intrinsically bounded after
# back-transformation. Response scale and tree settings are selected inside
# reef-blocked inner folds contained wholly in each outer analysis set.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
})
source("scripts/formal_model_helpers.R")

set.seed(202408L)

output_dir <- "output/formal_models"
model_dir <- "output/models/formal"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

tuning_grid <- expand.grid(
    interaction_depth = c(1L, 2L, 3L),
    shrinkage = c(0.01, 0.03),
    minobs = c(5L, 10L),
    magnitude_mode = c('gaussian', 'logit_gaussian'),
    stringsAsFactors = FALSE
)
tree_grid <- c(100L, 250L, 500L, 1000L, 1500L)
maximum_trees <- max(tree_grid)

fit_component <- function(formula, data, response, distribution, settings,
                          n_trees, seed) {
    outcome <- data[[response]]
    if (length(unique(outcome)) < 2L || sd(outcome) == 0) {
        return(list(type = "constant", value = mean(outcome)))
    }

    effective_minobs <- min(
        settings$minobs,
        max(1L, floor(nrow(data) / 5L))
    )
    set.seed(seed)
    model <- gbm(
        formula,
        data = data,
        distribution = distribution,
        n.trees = n_trees,
        interaction.depth = settings$interaction_depth,
        shrinkage = settings$shrinkage,
        n.minobsinnode = effective_minobs,
        bag.fraction = 0.70,
        train.fraction = 1,
        keep.data = FALSE,
        verbose = FALSE
    )
    list(type = "gbm", model = model)
}

predict_component <- function(component, new_data, n_trees, type) {
    if (component$type == "constant") {
        return(rep(component$value, nrow(new_data)))
    }
    predict(
        component$model,
        newdata = new_data,
        n.trees = n_trees,
        type = type
    )
}

fit_two_part_brt <- function(data, predictors, settings, n_trees, seed) {
    data$has_loss <- as.numeric(data$mortality_prop > 0)
    occurrence_formula <- as.formula(
        paste("has_loss ~", paste(predictors, collapse = " + "))
    )
    magnitude_formula <- as.formula(
        paste("mortality_prop ~", paste(predictors, collapse = " + "))
    )
    positive <- data[data$has_loss == 1, , drop = FALSE]
    magnitude_mode <- as.character(settings$magnitude_mode)
    if (magnitude_mode == 'logit_gaussian') {
        epsilon <- min(0.01, 0.5 / nrow(positive))
        positive$mortality_prop <- qlogis(
            pmin(pmax(positive$mortality_prop, epsilon), 1 - epsilon)
        )
    }

    model <- list(
        occurrence = fit_component(
            occurrence_formula, data, "has_loss", "bernoulli",
            settings, n_trees, seed
        ),
        magnitude = fit_component(
            magnitude_formula, positive, "mortality_prop", "gaussian",
            settings, n_trees, seed + 1L
        )
    )
    model$magnitude_mode <- magnitude_mode
    model
}

predict_two_part_brt <- function(model, new_data, n_trees) {
    occurrence <- predict_component(
        model$occurrence, new_data, n_trees, "response"
    )
    magnitude <- predict_component(
        model$magnitude, new_data, n_trees, "response"
    )
    if (model$magnitude_mode == 'logit_gaussian') {
        magnitude <- plogis(magnitude)
    }
    occurrence <- pmin(pmax(occurrence, 0), 1)
    magnitude <- pmin(pmax(magnitude, 0), 1)
    tibble(
        predicted_occurrence = occurrence,
        predicted_positive_loss = magnitude,
        predicted_mortality = occurrence * magnitude
    )
}

tune_brt <- function(data, predictors, seed) {
    inner_folds <- make_balanced_inner_folds(data, folds = 4L)
    tuning_results <- tibble()

    for (setting_index in seq_len(nrow(tuning_grid))) {
        settings <- tuning_grid[setting_index, ]
        for (inner_fold in sort(unique(inner_folds))) {
            inner_assessment <- inner_folds == inner_fold
            inner_analysis_data <- data[!inner_assessment, , drop = FALSE]
            inner_assessment_data <- data[inner_assessment, , drop = FALSE]
            prepared <- prepare_fold_predictors(
                inner_analysis_data,
                inner_assessment_data,
                predictors
            )
            inner_model <- fit_two_part_brt(
                prepared$analysis,
                predictors,
                settings,
                maximum_trees,
                seed + setting_index * 100L + inner_fold
            )

            for (n_trees in tree_grid) {
                predictions <- predict_two_part_brt(
                    inner_model, prepared$assessment, n_trees
                )
                tuning_results <- bind_rows(
                    tuning_results,
                    tibble(
                        interaction_depth = settings$interaction_depth,
                        shrinkage = settings$shrinkage,
                        minobs = settings$minobs,
                        magnitude_mode = settings$magnitude_mode,
                        n_trees = n_trees,
                        inner_fold = inner_fold,
                        n = nrow(prepared$assessment),
                        squared_error = sum(
                            (prepared$assessment$mortality_prop -
                             predictions$predicted_mortality)^2
                        ),
                        absolute_error = sum(abs(
                            prepared$assessment$mortality_prop -
                            predictions$predicted_mortality
                        ))
                    )
                )
            }
        }
    }

    scores <- tuning_results |>
        group_by(
            interaction_depth, shrinkage, minobs, magnitude_mode, n_trees
        ) |>
        summarise(
            n = sum(n),
            rmse = sqrt(sum(squared_error) / n),
            mae = sum(absolute_error) / n,
            .groups = "drop"
        ) |>
        arrange(
            rmse, mae, interaction_depth, n_trees,
            desc(shrinkage), minobs, magnitude_mode
        )

    list(
        best = scores[1, , drop = FALSE],
        scores = scores,
        fold_results = tuning_results
    )
}

settings_from_best <- function(best) {
    list(
        interaction_depth = as.integer(best$interaction_depth[[1]]),
        shrinkage = best$shrinkage[[1]],
        minobs = as.integer(best$minobs[[1]]),
        magnitude_mode = as.character(best$magnitude_mode[[1]])
    )
}

all_predictions <- tibble()
all_tuning <- tibble()

for (programme_key in names(validation_files)) {
    data <- load_programme_rows(programme_key)
    predictors <- programme_predictors(programme_key)

    for (scheme in all_validation_schemes) {
        for (fold in fold_values(data, scheme)) {
            cache_path <- file.path(
                model_dir,
                paste0(
                    "brt_", programme_key, "_", scheme, "_",
                    safe_fold_label(fold), ".rds"
                )
            )
            assessment_index <- assessment_rows(data, scheme, fold)
            analysis <- data[!assessment_index, , drop = FALSE]
            assessment <- data[assessment_index, , drop = FALSE]
            prepared <- prepare_fold_predictors(analysis, assessment, predictors)

            cached <- if (file.exists(cache_path)) readRDS(cache_path) else NULL
            cache_is_current <- !is.null(cached) &&
                identical(cached$predictors, predictors)
            if (cache_is_current) {
                predictions <- cached$predictions
                tuning <- cached$tuning
                cat("Loaded cached BRT fold:", programme_key, scheme, fold, "\n")
            } else {
                cat("Tuning BRT fold:", programme_key, scheme, fold,
                    "with", nrow(analysis), "analysis rows\n")
                flush.console()
                tuning <- tune_brt(
                    analysis,
                    predictors,
                    seed = 202408L + nrow(all_tuning) + 1L
                )
                best <- tuning$best
                model <- fit_two_part_brt(
                    prepared$analysis,
                    predictors,
                    settings_from_best(best),
                    as.integer(best$n_trees[[1]]),
                    seed = 202408L
                )
                predictions <- predict_two_part_brt(
                    model,
                    prepared$assessment,
                    as.integer(best$n_trees[[1]])
                )
                saveRDS(
                    list(
                        model = model,
                        preprocessing = prepared$preprocessing,
                        predictors = predictors,
                        tuning = tuning,
                        predictions = predictions
                    ),
                    cache_path
                )
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
                        observed_mortality = mortality_prop
                    ) |>
                    bind_cols(predictions) |>
                    mutate(
                        learner = "brt",
                        scheme = scheme,
                        fold = as.character(fold)
                    )
            )
            all_tuning <- bind_rows(
                all_tuning,
                tuning$best |>
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
    production_cache <- file.path(
        model_dir, paste0("brt_", programme_key, "_production.rds")
    )
    production_cached <- if (file.exists(production_cache)) {
        readRDS(production_cache)
    } else NULL
    production_cache_is_current <- !is.null(production_cached) &&
        identical(production_cached$predictors, predictors)
    if (production_cache_is_current) {
        cat("Loaded cached production BRT model:", programme_key, "\n")
    } else {
        cat("Tuning production BRT model:", programme_key,
            "on all", nrow(data), "rows including 2024\n")
        flush.console()
        tuning <- tune_brt(data, predictors, seed = 202408L)
        best <- tuning$best
        production_model <- fit_two_part_brt(
            full_prepared$analysis,
            predictors,
            settings_from_best(best),
            as.integer(best$n_trees[[1]]),
            seed = 202408L
        )
        saveRDS(
            list(
                model = production_model,
                preprocessing = full_prepared$preprocessing,
                predictors = predictors,
                tuning = tuning,
                event_years = sort(unique(data$event_year)),
                n = nrow(data)
            ),
            production_cache
        )
    }
}

if (anyDuplicated(paste(
    all_predictions$programme_key,
    all_predictions$source_observation_id,
    all_predictions$scheme
))) {
    stop("A BRT observation appears in multiple assessment folds")
}
if (any(!is.finite(all_predictions$predicted_mortality)) ||
    any(all_predictions$predicted_mortality < 0 | all_predictions$predicted_mortality > 1)) {
    stop("BRT produced invalid mortality predictions")
}

write_csv(all_predictions, file.path(output_dir, "brt_blocked_predictions.csv"))
write_csv(all_tuning, file.path(output_dir, "brt_tuning.csv"))
write_csv(
    prediction_metrics(all_predictions),
    file.path(output_dir, "brt_metrics.csv")
)

print(prediction_metrics(all_predictions))
