# Fit a monotone, composition-aware BRT candidate under leave-one-event-out
# validation. Tuning balances overall, event-balanced, and high-DHW errors.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
})
source('src/lib/shape_aware_model_helpers.R')

set.seed(202408L)
output_dir <- 'output/shape_aware_models'
model_dir <- 'output/models/shape_aware'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

tuning_grid <- expand.grid(
    interaction_depth = c(2L, 3L),
    shrinkage = c(0.01, 0.03),
    minobs = c(5L, 10L),
    magnitude_mode = c('gaussian', 'logit_gaussian'),
    stringsAsFactors = FALSE
)
tree_grid <- c(100L, 250L, 500L, 1000L)
maximum_trees <- max(tree_grid)

shape_brt_predictors <- function(programme_key) {
    c(
        'ann_maxdhw',
        setdiff(
            shape_predictors(programme_key),
            c('dhw_excess4', 'dhw_excess8')
        )
    )
}

reef_year_weights <- function(data) {
    counts <- data |>
        count(event_year, ReefID, name = 'reef_year_rows')
    weights <- data |>
        select(event_year, ReefID) |>
        left_join(counts, by = c('event_year', 'ReefID')) |>
        transmute(weight = 1 / reef_year_rows) |>
        pull(weight)
    weights / mean(weights)
}

fit_component <- function(formula, data, response, distribution, settings,
                          n_trees, weights, monotone, seed) {
    outcome <- data[[response]]
    if (length(unique(outcome)) < 2L || sd(outcome) == 0) {
        return(list(type = 'constant', value = weighted.mean(outcome, weights)))
    }
    effective_minobs <- min(
        settings$minobs,
        max(1L, floor(nrow(data) / 5L))
    )
    set.seed(seed)
    model <- do.call(
        gbm,
        list(
            formula = formula, data = data, distribution = distribution,
            weights = as.numeric(weights), var.monotone = monotone,
            n.trees = n_trees,
            interaction.depth = settings$interaction_depth,
            shrinkage = settings$shrinkage,
            n.minobsinnode = effective_minobs,
            bag.fraction = 0.70, train.fraction = 1,
            keep.data = FALSE, verbose = FALSE
        )
    )
    list(type = 'gbm', model = model)
}

fit_shape_brt <- function(data, predictors, settings, n_trees, seed) {
    data$has_loss <- as.numeric(data$mortality_prop > 0)
    weights <- reef_year_weights(data)
    monotone <- ifelse(
        predictors %in% c('ann_maxdhw', 'prop_acropora_pre'), 1L, 0L
    )
    occurrence_formula <- as.formula(
        paste('has_loss ~', paste(predictors, collapse = ' + '))
    )
    magnitude_formula <- as.formula(
        paste('mortality_prop ~', paste(predictors, collapse = ' + '))
    )
    positive <- data$has_loss == 1
    positive_data <- data[positive, , drop = FALSE]
    positive_weights <- weights[positive]
    magnitude_mode <- as.character(settings$magnitude_mode)
    if (magnitude_mode == 'logit_gaussian') {
        epsilon <- min(0.01, 0.5 / nrow(positive_data))
        positive_data$mortality_prop <- qlogis(pmin(pmax(
            positive_data$mortality_prop, epsilon
        ), 1 - epsilon))
    }

    list(
        occurrence = fit_component(
            occurrence_formula, data, 'has_loss', 'bernoulli', settings,
            n_trees, weights, monotone, seed
        ),
        magnitude = fit_component(
            magnitude_formula, positive_data, 'mortality_prop', 'gaussian',
            settings, n_trees, positive_weights, monotone, seed + 1L
        ),
        magnitude_mode = magnitude_mode
    )
}

predict_component <- function(component, new_data, n_trees, type) {
    if (component$type == 'constant') {
        return(rep(component$value, nrow(new_data)))
    }
    predict(
        component$model, newdata = new_data,
        n.trees = n_trees, type = type
    )
}

predict_shape_brt <- function(model, new_data, n_trees) {
    occurrence <- predict_component(
        model$occurrence, new_data, n_trees, 'response'
    )
    magnitude <- predict_component(
        model$magnitude, new_data, n_trees, 'response'
    )
    if (model$magnitude_mode == 'logit_gaussian') magnitude <- plogis(magnitude)
    occurrence <- pmin(pmax(occurrence, 0), 1)
    magnitude <- pmin(pmax(magnitude, 0), 1)
    tibble(
        predicted_occurrence = occurrence,
        predicted_positive_loss = magnitude,
        predicted_mortality = occurrence * magnitude
    )
}

tune_shape_brt <- function(data, predictors, seed) {
    inner_folds <- make_balanced_inner_folds(data, folds = 4L)
    results <- tibble()

    for (setting_index in seq_len(nrow(tuning_grid))) {
        settings <- tuning_grid[setting_index, ]
        for (inner_fold in sort(unique(inner_folds))) {
            assessment_index <- inner_folds == inner_fold
            analysis <- data[!assessment_index, , drop = FALSE]
            assessment <- data[assessment_index, , drop = FALSE]
            prepared <- prepare_fold_predictors(
                analysis, assessment, predictors
            )
            model <- fit_shape_brt(
                prepared$analysis, predictors, settings,
                maximum_trees, seed + 100L * setting_index + inner_fold
            )

            for (n_trees in tree_grid) {
                prediction <- predict_shape_brt(
                    model, prepared$assessment, n_trees
                )$predicted_mortality
                squared_error <- (
                    prepared$assessment$mortality_prop - prediction
                )^2
                event_sizes <- table(prepared$assessment$event_year)
                event_weights <- 1 / as.numeric(event_sizes[
                    as.character(prepared$assessment$event_year)
                ])
                high_dhw <- prepared$assessment$ann_maxdhw >= 10
                squared_error_sum <- sum(squared_error)
                event_weighted_error_sum <- sum(
                    squared_error * event_weights
                )
                high_dhw_error_sum <- sum(squared_error[high_dhw])
                results <- bind_rows(
                    results,
                    tibble(
                        interaction_depth = settings$interaction_depth,
                        shrinkage = settings$shrinkage,
                        minobs = settings$minobs,
                        magnitude_mode = settings$magnitude_mode,
                        n_trees, inner_fold,
                        n = length(squared_error),
                        squared_error = squared_error_sum,
                        event_weighted_error = event_weighted_error_sum,
                        event_groups = length(event_sizes),
                        high_dhw_error = high_dhw_error_sum,
                        high_dhw_n = sum(high_dhw)
                    )
                )
            }
        }
    }

    scores <- results |>
        group_by(
            interaction_depth, shrinkage, minobs,
            magnitude_mode, n_trees
        ) |>
        summarise(
            n = sum(n),
            rmse = sqrt(sum(squared_error) / n),
            event_rmse = sqrt(sum(event_weighted_error) / sum(event_groups)),
            high_dhw_rmse = if_else(
                sum(high_dhw_n) > 0,
                sqrt(sum(high_dhw_error) / sum(high_dhw_n)),
                NA_real_
            ),
            selection_score = rowMeans(
                cbind(rmse, event_rmse, high_dhw_rmse), na.rm = TRUE
            ),
            .groups = 'drop'
        ) |>
        arrange(
            selection_score, rmse, interaction_depth,
            n_trees, desc(shrinkage), minobs
        )
    list(best = scores[1, , drop = FALSE], scores = scores, folds = results)
}

settings_from_best <- function(best) {
    list(
        interaction_depth = as.integer(best$interaction_depth[[1]]),
        shrinkage = best$shrinkage[[1]],
        minobs = as.integer(best$minobs[[1]]),
        magnitude_mode = as.character(best$magnitude_mode[[1]])
    )
}

if (identical(Sys.getenv('SHAPE_BRT_DRY_RUN'), '1')) {
    data <- load_shape_rows('ltmp')
    predictors <- shape_brt_predictors('ltmp')
    prepared <- prepare_fold_predictors(data, data, predictors)
    settings <- list(
        interaction_depth = 2L, shrinkage = 0.03,
        minobs = 5L, magnitude_mode = 'logit_gaussian'
    )
    model <- fit_shape_brt(
        prepared$analysis, predictors, settings, 25L, 202408L
    )
    stopifnot(all(is.finite(
        predict_shape_brt(model, prepared$assessment, 25L)$predicted_mortality
    )))
    cat('Shape-aware BRT dry run passed.\n')
    quit(save = 'no', status = 0)
}

programme_keys <- names(validation_files)
requested_programme <- trimws(Sys.getenv('SHAPE_PROGRAMME_KEY', ''))
if (nzchar(requested_programme)) programme_keys <- requested_programme
all_predictions <- tibble()
all_tuning <- tibble()

for (programme_key in programme_keys) {
    data <- load_shape_rows(programme_key)
    predictors <- shape_brt_predictors(programme_key)

    for (fold in sort(unique(data$event_year))) {
        assessment_index <- data$event_year == fold
        analysis <- data[!assessment_index, , drop = FALSE]
        assessment <- data[assessment_index, , drop = FALSE]
        prepared <- prepare_fold_predictors(analysis, assessment, predictors)
        cache_file <- file.path(
            model_dir,
            paste0('shape_brt_', programme_key, '_event_', fold, '.rds')
        )
        cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
        cache_current <- !is.null(cached) &&
            identical(cached$predictors, predictors) &&
            identical(cached$model_version, 'monotone_weighted_v2')

        if (cache_current) {
            model <- cached$model
            tuning <- cached$tuning
            predictions <- cached$predictions
            cat('Loaded shape-aware BRT:', programme_key, fold, '\n')
        } else {
            cat('Tuning shape-aware BRT:', programme_key, fold, '\n')
            tuning <- tune_shape_brt(
                prepared$analysis, predictors, 202408L + fold
            )
            best <- tuning$best
            model <- fit_shape_brt(
                prepared$analysis, predictors, settings_from_best(best),
                as.integer(best$n_trees[[1]]), 202408L + fold
            )
            predictions <- predict_shape_brt(
                model, prepared$assessment, as.integer(best$n_trees[[1]])
            )
            saveRDS(
                list(
                    model = model, tuning = tuning,
                    preprocessing = prepared$preprocessing,
                    predictors = predictors, predictions = predictions,
                    model_version = 'monotone_weighted_v2'
                ),
                cache_file
            )
        }

        all_predictions <- bind_rows(
            all_predictions,
            assessment |>
                transmute(
                    programme_key, source_observation_id, ReefID, ReefName,
                    event_year, region_block,
                    ann_maxdhw, prop_acropora_pre, observed_pre_cover,
                    observed_mortality = mortality_prop
                ) |>
                bind_cols(predictions) |>
                mutate(
                    learner = 'shape_brt',
                    scheme = 'leave_one_event_out', fold = as.character(fold)
                )
        )
        all_tuning <- bind_rows(
            all_tuning,
            tuning$best |>
                mutate(
                    programme_key, fold = as.character(fold),
                    n_analysis = nrow(analysis),
                    n_assessment = nrow(assessment)
                )
        )
    }

    prepared <- prepare_fold_predictors(data, data, predictors)
    cache_file <- file.path(
        model_dir, paste0('shape_brt_', programme_key, '_production.rds')
    )
    cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
    cache_current <- !is.null(cached) &&
        identical(cached$predictors, predictors) &&
        identical(cached$model_version, 'monotone_weighted_v2')
    if (!cache_current) {
        cat('Tuning production shape-aware BRT:', programme_key, '\n')
        tuning <- tune_shape_brt(data, predictors, 202408L)
        best <- tuning$best
        model <- fit_shape_brt(
            prepared$analysis, predictors, settings_from_best(best),
            as.integer(best$n_trees[[1]]), 202408L
        )
        saveRDS(
            list(
                model = model, tuning = tuning,
                preprocessing = prepared$preprocessing,
                predictors = predictors,
                event_years = sort(unique(data$event_year)), n = nrow(data),
                model_version = 'monotone_weighted_v2'
            ),
            cache_file
        )
    }
}

if (any(!is.finite(all_predictions$predicted_mortality))) {
    stop('Shape-aware BRT produced invalid predictions.')
}
write_csv(
    all_predictions,
    file.path(output_dir, 'shape_brt_event_predictions.csv')
)
write_csv(
    all_tuning,
    file.path(output_dir, 'shape_brt_tuning.csv')
)
write_csv(
    prediction_metrics(all_predictions),
    file.path(output_dir, 'shape_brt_metrics.csv')
)
print(prediction_metrics(all_predictions))
