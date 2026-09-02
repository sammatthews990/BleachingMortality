# Cross-validated severe-loss expert for annual coral-cover transitions.
# The expert is deliberately separate from the general decomposed model. A
# Bernoulli gate estimates the probability of an annual loss >=10 percentage
# points, while a magnitude learner estimates change within that severe tail.
# Soft mixtures test whether the specialist improves extremes without causing
# unacceptable false alarms in ordinary years.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
})
source('src/lib/joint_compound_model_helpers.R')

set.seed(202416L)
output_dir <- 'output/extreme_loss_ensemble'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
n_trees <- 800L

data <- readRDS('data/processed/annual_coral_transitions.rds') |>
    mutate(
        severe_loss = as.numeric(cover_change_pp <= -10),
        tail_change_pp = cover_change_pp,
        dhw_excess4 = pmax(sst_maxdhw - 4, 0),
        dhw_excess8 = pmax(sst_maxdhw - 8, 0),
        acropora_cover_pre = pre_cover * prop_acropora_pre,
        previous_change_available = as.numeric(previous_change_available),
        wqc_available = as.numeric(wqc_available),
        cots_available = as.numeric(cots_available),
        mortality_environment_available = as.numeric(
            mortality_environment_available
        ),
        weather_available = as.numeric(weather_available),
        model_source = 'annual transition'
    )

mortality_auxiliary <- load_joint_compound_rows() |>
    mutate(auxiliary_change_pp = -100 * observed_pre_cover * mortality_prop) |>
    transmute(
        programme_key, ReefID, ReefName, event_year, region_block,
        pre_cover = observed_pre_cover,
        interval_years = interval_days / 365.25,
        prior_change_pp = NA_real_, prop_acropora_pre,
        acropora_cover_pre = observed_pre_cover * prop_acropora_pre,
        depth, sst_maxdhw = ann_maxdhw, dhw_excess4, dhw_excess8,
        histmDHW6, yrsince6, histmDHW4 = NA_real_, yrsince4 = NA_real_,
        dhw10_load4, dhw_novelty10, ann_maxsst = NA_real_,
        winyear_mean = NA_real_, winyear_sd = NA_real_,
        mcur_90 = NA_real_, dist_to_er_km = NA_real_,
        secc3m = NA_real_, secc3m_p10, cloudp_90,
        wqc_freqcc12, wqc_prior10_percentile,
        log1p_cyc_maxHrs4mw, log1p_cot_idwmeanpertow,
        log_coastal_rain30, era5_wind_mean,
        era5_wind_calm_fraction, era5_coastal_distance_km,
        previous_change_available = 0,
        wqc_available = as.numeric(is.finite(wqc_freqcc12)),
        cots_available = as.numeric(is.finite(log1p_cot_idwmeanpertow)),
        mortality_environment_available = 1,
        weather_available = 1,
        severe_loss = as.numeric(auxiliary_change_pp <= -10),
        tail_change_pp = auxiliary_change_pp,
        model_source = 'restricted bleaching mortality'
    )

numeric_predictors <- c(
    'pre_cover', 'interval_years', 'prior_change_pp',
    'prop_acropora_pre', 'acropora_cover_pre', 'depth',
    'sst_maxdhw', 'dhw_excess4', 'dhw_excess8',
    'histmDHW6', 'yrsince6', 'histmDHW4', 'yrsince4',
    'dhw10_load4', 'dhw_novelty10', 'ann_maxsst',
    'winyear_mean', 'winyear_sd', 'mcur_90', 'dist_to_er_km',
    'secc3m', 'secc3m_p10', 'cloudp_90',
    'wqc_freqcc12', 'wqc_prior10_percentile',
    'log1p_cyc_maxHrs4mw', 'log1p_cot_idwmeanpertow',
    'log_coastal_rain30', 'era5_wind_mean',
    'era5_wind_calm_fraction', 'era5_coastal_distance_km',
    'previous_change_available', 'wqc_available', 'cots_available',
    'mortality_environment_available', 'weather_available'
)

interaction_predictors <- c(
    'dhw4_x_acropora', 'dhw8_x_acropora',
    'dhw4_x_wqc', 'dhw4_x_relative_wqc',
    'dhw4_x_rain', 'dhw4_x_calm',
    'dhw4_x_secchi', 'dhw4_x_cloud'
)

model_predictors <- c(
    'programme_factor', 'region_factor', numeric_predictors,
    interaction_predictors
)

add_interactions <- function(rows) {
    rows |>
        mutate(
            dhw4_x_acropora = dhw_excess4 * prop_acropora_pre,
            dhw8_x_acropora = dhw_excess8 * prop_acropora_pre,
            dhw4_x_wqc = dhw_excess4 * wqc_freqcc12,
            dhw4_x_relative_wqc =
                dhw_excess4 * wqc_prior10_percentile,
            dhw4_x_rain = dhw_excess4 * log_coastal_rain30,
            dhw4_x_calm = dhw_excess4 * era5_wind_calm_fraction,
            dhw4_x_secchi = dhw_excess4 * secc3m_p10,
            dhw4_x_cloud = dhw_excess4 * cloudp_90,
            programme_factor = factor(
                programme_key, levels = c('ltmp', 'manta', 'mmp')
            ),
            region_factor = factor(
                region_block,
                levels = c('Northern GBR', 'Central GBR', 'Southern GBR')
            )
        )
}

assign_reef_folds <- function(rows, folds = 5L) {
    sizes <- rows |>
        count(ReefID, name = 'observations') |>
        arrange(desc(observations), ReefID)
    load <- rep(0L, folds)
    sizes$fold <- NA_integer_
    for (i in seq_len(nrow(sizes))) {
        chosen <- which.min(load)
        sizes$fold[[i]] <- chosen
        load[[chosen]] <- load[[chosen]] + sizes$observations[[i]]
    }
    sizes$fold[match(rows$ReefID, sizes$ReefID)]
}
data$reef_fold <- assign_reef_folds(data)

prepare_rows <- function(rows, rules) {
    for (i in seq_len(nrow(rules))) {
        predictor <- rules$predictor[[i]]
        rows[[predictor]][!is.finite(rows[[predictor]])] <- rules$median[[i]]
    }
    add_interactions(rows)
}

prepare_fold <- function(analysis, assessment, auxiliary) {
    rules <- tibble(predictor = numeric_predictors, median = NA_real_)
    for (i in seq_along(numeric_predictors)) {
        predictor <- numeric_predictors[[i]]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        if (length(observed) == 0L) stop('No training values for ', predictor)
        rules$median[[i]] <- median(observed)
    }
    list(
        analysis = prepare_rows(analysis, rules),
        assessment = prepare_rows(assessment, rules),
        auxiliary = prepare_rows(auxiliary, rules), rules = rules
    )
}

model_weights <- function(rows) {
    reef_n <- table(rows$ReefID)
    programme_n <- table(rows$programme_key)
    weight <- 1 / as.numeric(reef_n[rows$ReefID]) *
        1 / as.numeric(programme_n[rows$programme_key])
    if (n_distinct(rows$model_source) > 1L) {
        source_n <- table(rows$model_source)
        weight <- weight / as.numeric(source_n[rows$model_source])
    }
    weight / mean(weight)
}

fit_extreme_models <- function(annual, auxiliary, seed) {
    combined <- bind_rows(annual, auxiliary)
    severe <- combined |>
        filter(severe_loss == 1)
    if (nrow(severe) < 30L) stop('Too few severe rows for the tail expert')
    set.seed(seed)
    gate <- gbm(
        as.formula(paste(
            'severe_loss ~', paste(model_predictors, collapse = '+')
        )),
        data = combined, distribution = 'bernoulli',
        weights = model_weights(combined), n.trees = n_trees,
        interaction.depth = 3L, shrinkage = 0.025,
        n.minobsinnode = 10L, bag.fraction = 0.7,
        train.fraction = 1, keep.data = FALSE, verbose = FALSE
    )
    set.seed(seed + 1L)
    magnitude <- gbm(
        as.formula(paste(
            'tail_change_pp ~', paste(model_predictors, collapse = '+')
        )),
        data = severe, distribution = 'laplace',
        weights = model_weights(severe), n.trees = n_trees,
        interaction.depth = 2L, shrinkage = 0.025,
        n.minobsinnode = 5L, bag.fraction = 0.7,
        train.fraction = 1, keep.data = FALSE, verbose = FALSE
    )
    list(gate = gate, magnitude = magnitude, severe_rows = nrow(severe))
}

predict_extreme_models <- function(models, assessment, base_prediction) {
    gate_probability <- pmin(pmax(predict(
        models$gate, assessment, n.trees = n_trees, type = 'response'
    ), 0), 1)
    tail_prediction <- predict(
        models$magnitude, assessment,
        n.trees = n_trees, type = 'response'
    )
    tail_prediction <- pmax(
        pmin(tail_prediction, -10), -100 * assessment$pre_cover
    )
    tibble(
        severe_probability = gate_probability,
        tail_expert_change_pp = tail_prediction,
        base_prediction_pp = base_prediction,
        soft_gate_prediction_pp = base_prediction + gate_probability *
            (tail_prediction - base_prediction),
        aggressive_gate_prediction_pp = base_prediction +
            sqrt(gate_probability) * (tail_prediction - base_prediction)
    )
}

base_predictions <- read_csv(
    'output/decomposed_annual_change/decomposed_change_predictions.csv',
    show_col_types = FALSE
) |>
    filter(learner == 'Mortality-augmented decomposed BRT') |>
    select(scheme, fold, transition_id, predicted_change_pp)

schemes <- c(
    'recent_year_holdout', 'reef_blocked_5fold', 'forward_2024_2025'
)
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
        held_out_reefs <- unique(data$ReefID[assessment_index])
        auxiliary_index <- switch(
            scheme,
            recent_year_holdout =
                mortality_auxiliary$event_year != as.integer(fold),
            reef_blocked_5fold =
                !mortality_auxiliary$ReefID %in% held_out_reefs,
            forward_2024_2025 = mortality_auxiliary$event_year <= 2023L
        )
        assessment <- data[assessment_index, , drop = FALSE]
        prepared <- prepare_fold(
            data[analysis_index, , drop = FALSE], assessment,
            mortality_auxiliary[auxiliary_index, , drop = FALSE]
        )
        base <- base_predictions |>
            filter(scheme == .env$scheme, fold == as.character(.env$fold)) |>
            right_join(
                select(assessment, transition_id),
                by = 'transition_id', relationship = 'one-to-one'
            ) |>
            arrange(match(transition_id, assessment$transition_id))
        if (any(!is.finite(base$predicted_change_pp))) {
            stop('Base predictions do not match extreme assessment rows')
        }
        models <- fit_extreme_models(
            prepared$analysis, prepared$auxiliary,
            202416L + match(scheme, schemes) * 100L + match(fold, folds)
        )
        prediction <- predict_extreme_models(
            models, prepared$assessment, base$predicted_change_pp
        )
        fold_rows <- bind_rows(
            tibble(
                learner = 'Mortality-augmented base',
                predicted_change_pp = prediction$base_prediction_pp
            ),
            tibble(
                learner = 'Soft-gated extreme ensemble',
                predicted_change_pp = prediction$soft_gate_prediction_pp
            ),
            tibble(
                learner = 'Aggressive extreme ensemble',
                predicted_change_pp = prediction$aggressive_gate_prediction_pp
            )
        ) |>
            group_by(learner) |>
            mutate(row_index = row_number()) |>
            ungroup() |>
            left_join(
                assessment |>
                    mutate(row_index = row_number()) |>
                    select(
                        row_index, transition_id, programme_key,
                        ReefID, ReefName, event_year, pre_cover, post_cover,
                        cover_change_pp, sst_maxdhw, prop_acropora_pre,
                        wqc_freqcc12, wqc_prior10_percentile,
                        log_coastal_rain30, era5_wind_calm_fraction,
                        cyc_maxHrs4mw, cot_idwmeanpertow,
                        disturbance_text
                    ),
                by = 'row_index', relationship = 'many-to-one'
            ) |>
            left_join(
                prediction |>
                    mutate(row_index = row_number()) |>
                    select(
                        row_index, severe_probability,
                        tail_expert_change_pp
                    ),
                by = 'row_index', relationship = 'many-to-one'
            ) |>
            mutate(
                predicted_post_cover = pmin(pmax(
                    pre_cover + predicted_change_pp / 100, 0
                ), 1),
                predicted_change_pp = 100 * (
                    predicted_post_cover - pre_cover
                ),
                error_pp = cover_change_pp - predicted_change_pp,
                scheme = scheme, fold = as.character(fold)
            ) |>
            select(-row_index)
        all_predictions <- bind_rows(all_predictions, fold_rows)
        preprocessing <- bind_rows(
            preprocessing,
            prepared$rules |>
                mutate(scheme = scheme, fold = as.character(fold))
        )
        cat('Completed extreme-loss fold:', scheme, fold, '\n')
    }
}

auc_rank <- function(observed, score) {
    positive <- observed == 1
    n_positive <- sum(positive)
    n_negative <- sum(!positive)
    if (n_positive == 0L || n_negative == 0L) return(NA_real_)
    (sum(rank(score, ties.method = 'average')[positive]) -
        n_positive * (n_positive + 1) / 2) /
        (n_positive * n_negative)
}

metrics <- all_predictions |>
    group_by(learner, scheme) |>
    summarise(
        transitions = n(), severe_n = sum(cover_change_pp <= -10),
        rmse_pp = sqrt(mean(error_pp^2)),
        mae_pp = mean(abs(error_pp)),
        predictive_r2 = 1 - sum(error_pp^2) /
            sum((cover_change_pp - mean(cover_change_pp))^2),
        severe_mae_pp = mean(abs(error_pp)[cover_change_pp <= -10]),
        severe_bias_pp = mean(
            predicted_change_pp[cover_change_pp <= -10] -
                cover_change_pp[cover_change_pp <= -10]
        ),
        nonsevere_mae_pp = mean(abs(error_pp)[cover_change_pp > -10]),
        predicted_severe_rate = mean(predicted_change_pp <= -10),
        severe_recall = mean(
            predicted_change_pp[cover_change_pp <= -10] <= -10
        ),
        false_extreme_rate = mean(
            predicted_change_pp[cover_change_pp > -10] <= -10
        ),
        gate_brier = mean(
            ((cover_change_pp <= -10) - severe_probability)^2
        ),
        gate_auc = auc_rank(
            as.numeric(cover_change_pp <= -10), severe_probability
        ),
        .groups = 'drop'
    )

event_metrics <- all_predictions |>
    group_by(learner, scheme, event_year) |>
    summarise(
        transitions = n(), severe_n = sum(cover_change_pp <= -10),
        observed_change_pp = mean(cover_change_pp),
        severe_observed_pp = if_else(
            severe_n > 0, mean(cover_change_pp[cover_change_pp <= -10]),
            NA_real_
        ),
        severe_predicted_pp = if_else(
            severe_n > 0, mean(predicted_change_pp[cover_change_pp <= -10]),
            NA_real_
        ),
        predicted_change_pp = mean(predicted_change_pp),
        rmse_pp = sqrt(mean(error_pp^2)),
        .groups = 'drop'
    )

write_csv(
    all_predictions,
    file.path(output_dir, 'extreme_ensemble_predictions.csv')
)
write_csv(metrics, file.path(output_dir, 'extreme_ensemble_metrics.csv'))
write_csv(
    event_metrics,
    file.path(output_dir, 'extreme_ensemble_event_metrics.csv')
)
write_csv(
    preprocessing,
    file.path(output_dir, 'extreme_ensemble_preprocessing.csv')
)

outliers <- all_predictions |>
    filter(
        scheme == 'reef_blocked_5fold', event_year == 2024L,
        learner == 'Aggressive extreme ensemble'
    ) |>
    mutate(absolute_error_pp = abs(error_pp)) |>
    arrange(desc(absolute_error_pp))
write_csv(
    outliers,
    file.path(output_dir, 'extreme_2024_outliers.csv')
)

print(metrics)
print(filter(event_metrics, event_year == 2024L))
