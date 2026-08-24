# Joint multi-program BRT companion using the same response decomposition,
# compound-weather variables, outer folds, and reef-event weighting as BRMS.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
})
source('scripts/joint_compound_model_helpers.R')

set.seed(202410L)
output_dir <- 'output/joint_compound_models'
model_dir <- 'output/models/joint_compound'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

model_variant <- trimws(Sys.getenv('JOINT_BRT_VARIANT', 'compound'))
allowed_variants <- c(
    'compound', 'core', 'rrn_wq', 'rrn_cyclone', 'rrn', 'compound_rrn',
    'rrn_relative', 'rrn_full',
    'compound_rrn_relative', 'compound_rrn_full'
)
if (!model_variant %in% allowed_variants) {
    stop('JOINT_BRT_VARIANT must be one of: ', paste(allowed_variants, collapse = ', '))
}
selected_predictors <- switch(
    model_variant,
    compound = joint_compound_predictors,
    core = joint_core_predictors,
    rrn_wq = joint_rrn_wq_predictors,
    rrn_cyclone = joint_rrn_cyclone_predictors,
    rrn = joint_rrn_predictors,
    compound_rrn = joint_compound_rrn_predictors,
    rrn_relative = joint_rrn_relative_predictors,
    rrn_full = joint_rrn_full_predictors,
    compound_rrn_relative = joint_compound_rrn_relative_predictors,
    compound_rrn_full = joint_compound_rrn_full_predictors
)
joint_brt_predictors <- c(
    'programme_factor',
    setdiff(selected_predictors, c('dhw_excess4', 'dhw_excess8'))
)

tuning_grid <- expand.grid(
    interaction_depth = c(2L, 3L),
    shrinkage = c(0.01, 0.03),
    minobs = c(5L, 10L),
    magnitude_mode = c('gaussian', 'logit_gaussian'),
    stringsAsFactors = FALSE
)
tree_grid <- c(100L, 250L, 500L, 1000L)
maximum_trees <- max(tree_grid)
model_version <- paste0('joint_', model_variant, '_brt_v3')
file_tag <- switch(
    model_variant,
    compound = 'joint_brt',
    core = 'joint_brt_core',
    paste0('joint_brt_', model_variant)
)

reef_event_weights <- function(data) {
    sizes <- data |>
        count(reef_event_effect, name = 'rows')
    weights <- data |>
        select(reef_event_effect) |>
        left_join(sizes, by = 'reef_event_effect') |>
        transmute(weight = 1 / rows) |>
        pull(weight)
    weights / mean(weights)
}

fit_component <- function(formula, data, response, distribution, settings,
                          n_trees, weights, monotone, seed) {
    outcome <- data[[response]]
    if (length(unique(outcome)) < 2L || sd(outcome) == 0) {
        return(list(type = 'constant', value = weighted.mean(outcome, weights)))
    }
    effective_minobs <- min(settings$minobs, max(1L, floor(nrow(data) / 5L)))
    set.seed(seed)
    model <- do.call(gbm, list(
        formula = formula, data = data, distribution = distribution,
        weights = as.numeric(weights), var.monotone = monotone,
        n.trees = n_trees, interaction.depth = settings$interaction_depth,
        shrinkage = settings$shrinkage,
        n.minobsinnode = effective_minobs,
        bag.fraction = 0.70, train.fraction = 1,
        keep.data = FALSE, verbose = FALSE
    ))
    list(type = 'gbm', model = model)
}

fit_joint_brt <- function(data, settings, n_trees, seed) {
    data$has_loss <- as.numeric(data$mortality_prop > 0)
    weights <- reef_event_weights(data)
    monotone <- ifelse(
        joint_brt_predictors %in% c('ann_maxdhw', 'prop_acropora_pre'),
        1L, 0L
    )
    formula_occurrence <- as.formula(paste(
        'has_loss ~', paste(joint_brt_predictors, collapse = ' + ')
    ))
    formula_magnitude <- as.formula(paste(
        'mortality_prop ~', paste(joint_brt_predictors, collapse = ' + ')
    ))
    positive <- data$has_loss == 1
    magnitude_data <- data[positive, , drop = FALSE]
    magnitude_weights <- weights[positive]
    magnitude_mode <- as.character(settings$magnitude_mode)
    if (magnitude_mode == 'logit_gaussian') {
        epsilon <- min(0.01, 0.5 / nrow(magnitude_data))
        magnitude_data$mortality_prop <- qlogis(pmin(pmax(
            magnitude_data$mortality_prop, epsilon
        ), 1 - epsilon))
    }
    list(
        occurrence = fit_component(
            formula_occurrence, data, 'has_loss', 'bernoulli', settings,
            n_trees, weights, monotone, seed
        ),
        magnitude = fit_component(
            formula_magnitude, magnitude_data, 'mortality_prop', 'gaussian',
            settings, n_trees, magnitude_weights, monotone, seed + 1L
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

predict_joint_brt <- function(model, new_data, n_trees) {
    occurrence <- predict_component(
        model$occurrence, new_data, n_trees, 'response'
    )
    positive <- predict_component(
        model$magnitude, new_data, n_trees, 'response'
    )
    if (model$magnitude_mode == 'logit_gaussian') positive <- plogis(positive)
    occurrence <- pmin(pmax(occurrence, 0), 1)
    positive <- pmin(pmax(positive, 0), 1)
    tibble(
        predicted_occurrence = occurrence,
        predicted_positive_loss = positive,
        predicted_mortality = occurrence * positive
    )
}

tune_joint_brt <- function(data, seed) {
    inner_folds <- assign_balanced_group_folds(
        data, 'reef_event_effect', folds = 4L
    )
    results <- tibble()
    for (setting_index in seq_len(nrow(tuning_grid))) {
        settings <- tuning_grid[setting_index, ]
        for (inner_fold in sort(unique(inner_folds))) {
            assessment_index <- inner_folds == inner_fold
            analysis <- data[!assessment_index, , drop = FALSE]
            assessment <- data[assessment_index, , drop = FALSE]
            prepared <- prepare_joint_predictors(analysis, assessment)
            model <- fit_joint_brt(
                prepared$analysis, settings, maximum_trees,
                seed + 100L * setting_index + inner_fold
            )
            for (n_trees in tree_grid) {
                prediction <- predict_joint_brt(
                    model, prepared$assessment, n_trees
                )$predicted_mortality
                squared_error <- (prepared$assessment$mortality_prop - prediction)^2
                severe <- prepared$assessment$mortality_prop >= 0.2
                programme_sizes <- table(prepared$assessment$programme_key)
                programme_weights <- 1 / as.numeric(programme_sizes[
                    prepared$assessment$programme_key
                ])
                event_sizes <- table(prepared$assessment$event_year)
                event_weights <- 1 / as.numeric(event_sizes[
                    as.character(prepared$assessment$event_year)
                ])
                results <- bind_rows(results, tibble(
                    interaction_depth = settings$interaction_depth,
                    shrinkage = settings$shrinkage,
                    minobs = settings$minobs,
                    magnitude_mode = settings$magnitude_mode,
                    n_trees, inner_fold, n = length(squared_error),
                    squared_error = sum(squared_error),
                    programme_error = sum(squared_error * programme_weights),
                    programme_groups = length(programme_sizes),
                    event_error = sum(squared_error * event_weights),
                    event_groups = length(event_sizes),
                    severe_error = sum(squared_error[severe]),
                    severe_n = sum(severe)
                ))
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
            programme_rmse = sqrt(
                sum(programme_error) / sum(programme_groups)
            ),
            event_rmse = sqrt(sum(event_error) / sum(event_groups)),
            severe_rmse = sqrt(sum(severe_error) / sum(severe_n)),
            selection_score = mean(c(
                rmse, programme_rmse, event_rmse, severe_rmse
            )),
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

data <- load_joint_compound_rows()
if (identical(Sys.getenv('JOINT_BRT_DRY_RUN'), '1')) {
    prepared <- prepare_joint_predictors(data, data)
    settings <- list(
        interaction_depth = 2L, shrinkage = 0.03,
        minobs = 5L, magnitude_mode = 'logit_gaussian'
    )
    model <- fit_joint_brt(prepared$analysis, settings, 25L, 202410L)
    prediction <- predict_joint_brt(model, prepared$assessment, 25L)
    stopifnot(all(is.finite(prediction$predicted_mortality)))
    cat('Joint ', model_variant, ' BRT dry run passed.\n', sep = '')
    quit(save = 'no', status = 0)
}

schemes <- c('leave_one_event_out', 'reef_blocked_2024')
requested_scheme <- trimws(Sys.getenv('JOINT_BRT_SCHEME', ''))
if (nzchar(requested_scheme)) schemes <- requested_scheme
all_predictions <- tibble()
all_tuning <- tibble()

for (scheme in schemes) {
    folds <- if (scheme == 'leave_one_event_out') {
        sort(unique(data$event_year))
    } else {
        sort(unique(data$joint_reef_fold[data$event_year == 2024L]))
    }
    for (fold in folds) {
        assessment_index <- joint_assessment_index(data, scheme, fold)
        analysis_index <- joint_analysis_index(data, scheme, fold)
        analysis <- data[analysis_index, , drop = FALSE]
        assessment <- data[assessment_index, , drop = FALSE]
        prepared <- prepare_joint_predictors(analysis, assessment)
        cache_file <- file.path(
            model_dir, paste0(file_tag, '_', scheme, '_', fold, '.rds')
        )
        cached <- if (file.exists(cache_file)) readRDS(cache_file) else NULL
        cache_current <- !is.null(cached) &&
            identical(cached$model_version, model_version) &&
            identical(cached$predictors, joint_brt_predictors)
        if (cache_current) {
            prediction <- cached$predictions
            tuning <- cached$tuning
            cat('Loaded joint BRT:', scheme, fold, '\n')
        } else {
            cat('Tuning joint BRT:', scheme, fold, '\n')
            tuning <- tune_joint_brt(
                prepared$analysis, 202410L + as.integer(fold)
            )
            settings <- settings_from_best(tuning$best)
            model <- fit_joint_brt(
                prepared$analysis, settings,
                as.integer(tuning$best$n_trees[[1]]),
                203410L + as.integer(fold)
            )
            prediction <- predict_joint_brt(
                model, prepared$assessment,
                as.integer(tuning$best$n_trees[[1]])
            )
            cached <- list(
                model = model, predictions = prediction,
                preprocessing = prepared$preprocessing,
                predictors = joint_brt_predictors, tuning = tuning,
                model_version = model_version, scheme = scheme, fold = fold
            )
            saveRDS(cached, cache_file)
        }
        all_predictions <- bind_rows(
            all_predictions,
            assessment |>
                transmute(
                    programme_key, source_observation_id, ReefID, ReefName,
                    event_year, region_block, joint_reef_fold,
                    ann_maxdhw, prop_acropora_pre, observed_pre_cover,
                    log_coastal_rain30, era5_wind_mean,
                    era5_wind_calm_fraction, era5_coastal_distance_km,
                    wqc_freqcc12, wqc_prior10_percentile,
                    wqc_prior10_delta, cyc_maxHrs4mw,
                    log1p_cyc_maxHrs4mw, cot_meanpertow,
                    cot_idwmeanpertow, log1p_cot_idwmeanpertow,
                    DISTURBANCE_TYPE, storm_name, description,
                    tooltip, disturbance_text,
                    observed_mortality = mortality_prop
                ) |>
                bind_cols(prediction) |>
                mutate(
                    learner = paste0('joint_', model_variant, '_brt'),
                    scheme = scheme,
                    fold = as.character(fold),
                    residual = observed_mortality - predicted_mortality,
                    absolute_error = abs(residual)
                )
        )
        all_tuning <- bind_rows(
            all_tuning,
            tuning$best |> mutate(scheme = scheme, fold = as.character(fold))
        )
    }
}

prepared <- prepare_joint_predictors(data, data)
production_file <- file.path(model_dir, paste0(file_tag, '_production.rds'))
production <- if (file.exists(production_file)) readRDS(production_file) else NULL
production_current <- !is.null(production) &&
    identical(production$model_version, model_version) &&
    identical(production$predictors, joint_brt_predictors)
if (!production_current) {
    cat('Tuning joint BRT production model\n')
    tuning <- tune_joint_brt(prepared$analysis, 202410L)
    settings <- settings_from_best(tuning$best)
    model <- fit_joint_brt(
        prepared$analysis, settings,
        as.integer(tuning$best$n_trees[[1]]), 203410L
    )
    production <- list(
        model = model, preprocessing = prepared$preprocessing,
        predictors = joint_brt_predictors, tuning = tuning,
        model_version = model_version, n = nrow(data)
    )
    saveRDS(production, production_file)
}

prediction_name <- paste0(file_tag, '_predictions.csv')
tuning_name <- paste0(file_tag, '_tuning.csv')
write_csv(all_predictions, file.path(output_dir, prediction_name))
write_csv(all_tuning, file.path(output_dir, tuning_name))
print(all_predictions |>
    group_by(scheme, programme_key) |>
    summarise(
        n = n(), rmse = sqrt(mean(residual^2)),
        mae = mean(absolute_error),
        severe_mae = mean(absolute_error[observed_mortality >= 0.2]),
        bias = mean(residual), .groups = 'drop'
    ))
