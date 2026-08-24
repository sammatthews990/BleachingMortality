# Initial operational benchmarks for annual coral-cover dynamics. These models
# retain gains and losses and predict both next cover and percentage-point
# change. They are intentionally compact; the purpose is to establish whether
# a dynamic target improves prediction before adding a full Bayesian state
# model.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
    library(splines)
})

set.seed(202413L)
output_dir <- 'output/annual_coral_change'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- readRDS('data/processed/annual_coral_transitions.rds') |>
    mutate(
        programme_factor = factor(
            programme_key, levels = c('ltmp', 'manta', 'mmp')
        ),
        region_factor = factor(
            region_block,
            levels = c('Northern GBR', 'Central GBR', 'Southern GBR')
        ),
        previous_change_available = as.numeric(previous_change_available),
        wqc_available = as.numeric(wqc_available),
        cots_available = as.numeric(cots_available)
    )

numeric_predictors <- c(
    'pre_cover', 'interval_years', 'prior_change_pp',
    'prop_acropora_pre', 'depth', 'sst_maxdhw',
    'wqc_freqcc12', 'wqc_prior10_percentile',
    'log1p_cyc_maxHrs4mw', 'log1p_cot_idwmeanpertow',
    'previous_change_available', 'wqc_available', 'cots_available'
)
model_predictors <- c(
    'programme_factor', 'region_factor', numeric_predictors
)

assign_reef_folds <- function(rows, folds = 5L) {
    sizes <- rows |>
        count(ReefID, name = 'observations') |>
        arrange(desc(observations), ReefID)
    load <- rep(0L, folds)
    sizes$fold <- NA_integer_
    for (i in seq_len(nrow(sizes))) {
        selected <- which.min(load)
        sizes$fold[[i]] <- selected
        load[[selected]] <- load[[selected]] + sizes$observations[[i]]
    }
    sizes$fold[match(rows$ReefID, sizes$ReefID)]
}
data$reef_fold <- assign_reef_folds(data)

prepare_fold <- function(analysis, assessment) {
    rules <- tibble(predictor = numeric_predictors, median = NA_real_)
    for (i in seq_along(numeric_predictors)) {
        predictor <- numeric_predictors[[i]]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        if (length(observed) == 0L) stop('No training values for ', predictor)
        value <- median(observed)
        analysis[[predictor]][!is.finite(analysis[[predictor]])] <- value
        assessment[[predictor]][!is.finite(assessment[[predictor]])] <- value
        rules$median[[i]] <- value
    }
    analysis$programme_factor <- factor(
        analysis$programme_key, levels = c('ltmp', 'manta', 'mmp')
    )
    assessment$programme_factor <- factor(
        assessment$programme_key, levels = c('ltmp', 'manta', 'mmp')
    )
    analysis$region_factor <- factor(
        analysis$region_block,
        levels = c('Northern GBR', 'Central GBR', 'Southern GBR')
    )
    assessment$region_factor <- factor(
        assessment$region_block,
        levels = c('Northern GBR', 'Central GBR', 'Southern GBR')
    )
    list(analysis = analysis, assessment = assessment, rules = rules)
}

balanced_weights <- function(rows) {
    reef_n <- table(rows$ReefID)
    programme_n <- table(rows$programme_key)
    weight <- 1 / as.numeric(reef_n[rows$ReefID]) *
        1 / as.numeric(programme_n[rows$programme_key])
    weight / mean(weight)
}

fit_brt <- function(rows, response, seed) {
    set.seed(seed)
    gbm(
        as.formula(paste(response, '~', paste(model_predictors, collapse = '+'))),
        data = rows, distribution = 'gaussian',
        weights = balanced_weights(rows),
        n.trees = 1200L, interaction.depth = 3L,
        shrinkage = 0.02, n.minobsinnode = 10L,
        bag.fraction = 0.7, train.fraction = 1,
        keep.data = FALSE, verbose = FALSE
    )
}

schemes <- c('recent_year_holdout', 'reef_blocked_5fold', 'forward_2024_2025')
all_predictions <- tibble()
preprocessing <- tibble()

for (scheme in schemes) {
    folds <- switch(
        scheme,
        recent_year_holdout = 2016:2025,
        reef_blocked_5fold = 1:5,
        forward_2024_2025 = '2024_2025'
    )
    for (fold in folds) {
        assessment_index <- switch(
            scheme,
            recent_year_holdout = data$event_year == as.integer(fold),
            reef_blocked_5fold = data$reef_fold == as.integer(fold),
            forward_2024_2025 = data$event_year >= 2024L
        )
        analysis_index <- switch(
            scheme,
            recent_year_holdout = !assessment_index,
            reef_blocked_5fold = !assessment_index,
            forward_2024_2025 = data$event_year <= 2023L
        )
        analysis <- data[analysis_index, , drop = FALSE]
        assessment <- data[assessment_index, , drop = FALSE]
        prepared <- prepare_fold(analysis, assessment)
        seed <- 202413L + match(scheme, schemes) * 100L +
            match(fold, folds)

        post_model <- fit_brt(prepared$analysis, 'post_cover', seed)
        change_model <- fit_brt(
            prepared$analysis, 'cover_change_pp', seed + 1L
        )
        growth_model <- lm(
            cover_change_pp ~ programme_factor + region_factor +
                splines::ns(pre_cover, df = 3) + interval_years +
                prior_change_pp + prop_acropora_pre + sst_maxdhw +
                wqc_freqcc12 + log1p_cyc_maxHrs4mw +
                log1p_cot_idwmeanpertow,
            data = prepared$analysis,
            weights = balanced_weights(prepared$analysis)
        )

        prediction_post_brt <- pmin(pmax(predict(
            post_model, prepared$assessment,
            n.trees = 1200L, type = 'response'
        ), 0), 1)
        prediction_change_brt <- predict(
            change_model, prepared$assessment,
            n.trees = 1200L, type = 'response'
        )
        prediction_post_from_change <- pmin(pmax(
            prepared$assessment$pre_cover + prediction_change_brt / 100,
            0
        ), 1)
        prediction_growth <- predict(growth_model, prepared$assessment)
        prediction_post_growth <- pmin(pmax(
            prepared$assessment$pre_cover + prediction_growth / 100,
            0
        ), 1)

        fold_predictions <- bind_rows(
            tibble(
                learner = 'Persistence',
                predicted_post_cover = prepared$assessment$pre_cover
            ),
            tibble(
                learner = 'Growth/environment LM',
                predicted_post_cover = prediction_post_growth
            ),
            tibble(
                learner = 'Post-cover BRT',
                predicted_post_cover = prediction_post_brt
            ),
            tibble(
                learner = 'Change BRT',
                predicted_post_cover = prediction_post_from_change
            )
        ) |>
            group_by(learner) |>
            mutate(row_index = row_number()) |>
            ungroup() |>
            left_join(
                prepared$assessment |>
                    mutate(row_index = row_number()) |>
                    select(
                        row_index, transition_id, programme_key,
                        ReefID, ReefName, event_year, survey_date,
                        pre_cover, post_cover, cover_change_pp,
                        sst_maxdhw, prop_acropora_pre,
                        wqc_freqcc12, cyc_maxHrs4mw,
                        cot_idwmeanpertow, DISTURBANCE_TYPE,
                        disturbance_text
                    ),
                by = 'row_index', relationship = 'many-to-one'
            ) |>
            mutate(
                predicted_change_pp = 100 * (
                    predicted_post_cover - pre_cover
                ),
                post_cover_error = post_cover - predicted_post_cover,
                change_error_pp = cover_change_pp - predicted_change_pp,
                scheme = scheme, fold = as.character(fold)
            ) |>
            select(-row_index)

        all_predictions <- bind_rows(all_predictions, fold_predictions)
        preprocessing <- bind_rows(
            preprocessing,
            prepared$rules |>
                mutate(scheme = scheme, fold = as.character(fold))
        )
        cat('Completed annual transition fold:', scheme, fold, '\n')
    }
}

metric_summary <- function(rows, by_programme = FALSE) {
    groups <- c('learner', 'scheme')
    if (by_programme) groups <- c(groups, 'programme_key')
    rows |>
        group_by(across(all_of(groups))) |>
        summarise(
            transitions = n(),
            post_cover_rmse_pp = 100 * sqrt(mean(post_cover_error^2)),
            post_cover_mae_pp = 100 * mean(abs(post_cover_error)),
            change_rmse_pp = sqrt(mean(change_error_pp^2)),
            change_mae_pp = mean(abs(change_error_pp)),
            change_bias_pp = mean(change_error_pp),
            change_predictive_r2 = 1 -
                sum(change_error_pp^2) /
                sum((cover_change_pp - mean(cover_change_pp))^2),
            direction_accuracy = mean(
                sign(predicted_change_pp) == sign(cover_change_pp)
            ),
            severe_loss_mae_pp = mean(
                abs(change_error_pp)[cover_change_pp <= -10]
            ),
            .groups = 'drop'
        )
}

metrics <- bind_rows(
    metric_summary(all_predictions) |> mutate(programme_key = 'all'),
    metric_summary(all_predictions, by_programme = TRUE)
)
write_csv(
    all_predictions,
    file.path(output_dir, 'annual_change_predictions.csv')
)
write_csv(metrics, file.path(output_dir, 'annual_change_metrics.csv'))
write_csv(
    preprocessing,
    file.path(output_dir, 'annual_change_preprocessing.csv')
)

# All-data models are retained for interpretation and later operational use.
# Their apparent fit is not reported as validation performance.
prepared_all <- prepare_fold(data, data)
production_post_brt <- fit_brt(
    prepared_all$analysis, 'post_cover', 202413L
)
production_change_brt <- fit_brt(
    prepared_all$analysis, 'cover_change_pp', 202414L
)
production_growth_lm <- lm(
    cover_change_pp ~ programme_factor + region_factor +
        splines::ns(pre_cover, df = 3) + interval_years +
        prior_change_pp + prop_acropora_pre + sst_maxdhw +
        wqc_freqcc12 + log1p_cyc_maxHrs4mw +
        log1p_cot_idwmeanpertow,
    data = prepared_all$analysis,
    weights = balanced_weights(prepared_all$analysis)
)
saveRDS(
    list(
        post_cover_brt = production_post_brt,
        change_brt = production_change_brt,
        growth_lm = production_growth_lm,
        preprocessing = prepared_all$rules,
        predictors = model_predictors,
        n_trees = 1200L,
        model_version = 'annual_change_benchmarks_v1'
    ),
    file.path(output_dir, 'annual_change_production_models.rds')
)

importance <- bind_rows(
    summary(production_post_brt, plotit = FALSE) |>
        transmute(learner = 'Post-cover BRT', variable = var, influence = rel.inf),
    summary(production_change_brt, plotit = FALSE) |>
        transmute(learner = 'Change BRT', variable = var, influence = rel.inf)
)
write_csv(importance, file.path(output_dir, 'annual_change_importance.csv'))

outliers <- all_predictions |>
    filter(
        scheme == 'forward_2024_2025',
        learner %in% c('Growth/environment LM', 'Change BRT')
    ) |>
    mutate(absolute_change_error_pp = abs(change_error_pp)) |>
    group_by(learner) |>
    arrange(desc(absolute_change_error_pp), .by_group = TRUE) |>
    slice_head(n = 30L) |>
    ungroup()
write_csv(outliers, file.path(output_dir, 'annual_change_outliers.csv'))

cover_bins <- data |>
    mutate(
        pre_cover_bin = cut(
            pre_cover,
            breaks = c(-Inf, 0.05, 0.10, 0.20, 0.40, Inf),
            labels = c('<5%', '5-10%', '10-20%', '20-40%', '>40%')
        ),
        relative_change = cover_change_pp / pmax(100 * pre_cover, 1e-8)
    ) |>
    group_by(pre_cover_bin) |>
    summarise(
        transitions = n(),
        median_absolute_change_pp = median(abs(cover_change_pp)),
        relative_change_q05 = quantile(relative_change, 0.05),
        relative_change_q50 = quantile(relative_change, 0.50),
        relative_change_q95 = quantile(relative_change, 0.95),
        .groups = 'drop'
    )
write_csv(cover_bins, file.path(output_dir, 'change_by_initial_cover.csv'))

print(filter(metrics, programme_key == 'all'))
