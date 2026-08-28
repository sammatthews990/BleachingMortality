# Decompose approximately annual coral-cover transitions into an expected
# recovery contribution and an expected disturbance-loss contribution. The
# two processes are statistically identified from net gains and net losses;
# annual cover snapshots cannot reveal simultaneous gross growth and loss.
# This pragmatic decomposition is therefore a predictive transition model,
# not a claim that the latent ecological fluxes have been directly observed.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
})
source('scripts/joint_compound_model_helpers.R')

set.seed(202415L)
output_dir <- 'output/decomposed_annual_change'
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
        available_space = pmax(1 - pre_cover, 0),
        has_gain = as.numeric(cover_change_pp > 0),
        has_loss = as.numeric(cover_change_pp < 0),
        gain_fraction = pmin(pmax(
            (post_cover - pre_cover) / pmax(available_space, 1e-6), 0
        ), 1),
        loss_fraction = pmin(pmax(
            (pre_cover - post_cover) / pmax(pre_cover, 1e-6), 0
        ), 1),
        dhw_excess4 = pmax(sst_maxdhw - 4, 0),
        dhw_excess8 = pmax(sst_maxdhw - 8, 0),
        acropora_cover_pre = pre_cover * prop_acropora_pre,
        previous_change_available = as.numeric(previous_change_available),
        wqc_available = as.numeric(wqc_available),
        cots_available = as.numeric(cots_available),
        mortality_environment_available = as.numeric(
            mortality_environment_available
        ),
        weather_available = as.numeric(weather_available)
        ,model_source = 'annual transition'
    )

# The restricted mortality dataset is used as an auxiliary response for the
# loss process. This is stronger transfer than merely copying predictor names:
# it asks the annual loss learner to retain the bleaching-attributed response
# learned by the relative-mortality analysis. Its rows are excluded under the
# same temporal or reef holdout as the annual assessment rows.
mortality_auxiliary <- load_joint_compound_rows() |>
    transmute(
        programme_key, ReefID, ReefName, event_year, region_block,
        pre_cover = observed_pre_cover,
        available_space = pmax(1 - observed_pre_cover, 0),
        interval_years = interval_days / 365.25,
        prior_change_pp = NA_real_,
        prop_acropora_pre,
        acropora_cover_pre = observed_pre_cover * prop_acropora_pre,
        depth,
        sst_maxdhw = ann_maxdhw,
        dhw_excess4, dhw_excess8,
        histmDHW6, yrsince6,
        histmDHW4 = NA_real_, yrsince4 = NA_real_,
        dhw10_load4, dhw_novelty10, dhw_events_since2016_n6,
        dhw_events_prior8_n6 = NA_real_,
        dhw_years_since_last_n6_capped8, dhw_no_prior_n6,
        ann_maxsst = NA_real_,
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
        has_loss = as.numeric(mortality_prop > 0),
        loss_fraction = mortality_prop,
        model_source = 'restricted bleaching mortality'
    )

raw_numeric_predictors <- c(
    'pre_cover', 'available_space', 'interval_years', 'prior_change_pp',
    'prop_acropora_pre', 'acropora_cover_pre', 'depth', 'sst_maxdhw',
    'dhw_excess4', 'dhw_excess8',
    'histmDHW6', 'yrsince6', 'histmDHW4', 'yrsince4',
    'dhw10_load4', 'dhw_novelty10', 'dhw_events_since2016_n6',
    'dhw_events_prior8_n6', 'dhw_years_since_last_n6_capped8',
    'dhw_no_prior_n6', 'ann_maxsst',
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

growth_predictors <- c(
    'programme_factor', 'region_factor',
    'pre_cover', 'available_space', 'interval_years', 'prior_change_pp',
    'prop_acropora_pre', 'acropora_cover_pre', 'depth',
    'sst_maxdhw', 'dhw10_load4', 'dhw_novelty10',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8',
    'wqc_freqcc12', 'log1p_cyc_maxHrs4mw',
    'log1p_cot_idwmeanpertow',
    'previous_change_available', 'wqc_available', 'cots_available',
    'mortality_environment_available'
)

loss_core_predictors <- c(
    'programme_factor', 'region_factor',
    'pre_cover', 'interval_years', 'prior_change_pp',
    'prop_acropora_pre', 'acropora_cover_pre', 'depth',
    'sst_maxdhw', 'dhw_excess4', 'dhw_excess8',
    'previous_change_available'
)

# The rich loss process preserves the relative-mortality feature contract.
# Rainfall is labelled as an exposure proxy rather than measured salinity;
# cyclone waves and COTS remain separate disturbance mechanisms.
loss_modifier_predictors <- c(
    loss_core_predictors,
    'histmDHW6', 'yrsince6', 'histmDHW4', 'yrsince4',
    'dhw10_load4', 'dhw_novelty10', 'dhw_events_since2016_n6',
    'dhw_years_since_last_n6_capped8', 'dhw_no_prior_n6', 'ann_maxsst',
    'winyear_mean', 'winyear_sd', 'mcur_90', 'dist_to_er_km',
    'secc3m', 'secc3m_p10', 'cloudp_90',
    'wqc_freqcc12', 'wqc_prior10_percentile',
    'log1p_cyc_maxHrs4mw', 'log1p_cot_idwmeanpertow',
    'log_coastal_rain30', 'era5_wind_mean',
    'era5_wind_calm_fraction', 'era5_coastal_distance_km',
    'wqc_available', 'cots_available',
    'mortality_environment_available', 'weather_available',
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
            dhw4_x_cloud = dhw_excess4 * cloudp_90
        )
}

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

apply_preprocessing <- function(rows, rules) {
    for (i in seq_len(nrow(rules))) {
        predictor <- rules$predictor[[i]]
        rows[[predictor]][!is.finite(rows[[predictor]])] <- rules$median[[i]]
    }
    rows <- add_interactions(rows)
    rows$programme_factor <- factor(
        rows$programme_key, levels = c('ltmp', 'manta', 'mmp')
    )
    rows$region_factor <- factor(
        rows$region_block,
        levels = c('Northern GBR', 'Central GBR', 'Southern GBR')
    )
    rows
}

prepare_fold <- function(analysis, assessment, auxiliary = NULL) {
    rules <- tibble(predictor = raw_numeric_predictors, median = NA_real_)
    for (i in seq_along(raw_numeric_predictors)) {
        predictor <- raw_numeric_predictors[[i]]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        if (length(observed) == 0L) stop('No training values for ', predictor)
        value <- median(observed)
        rules$median[[i]] <- value
    }
    list(
        analysis = apply_preprocessing(analysis, rules),
        assessment = apply_preprocessing(assessment, rules),
        auxiliary = if (is.null(auxiliary)) NULL else
            apply_preprocessing(auxiliary, rules),
        rules = rules
    )
}

balanced_weights <- function(rows) {
    reef_n <- table(rows$ReefID)
    programme_n <- table(rows$programme_key)
    weight <- 1 / as.numeric(reef_n[rows$ReefID]) *
        1 / as.numeric(programme_n[rows$programme_key])
    weight / mean(weight)
}

fit_component <- function(rows, response, predictors, distribution, seed,
                          depth = 3L) {
    weights <- balanced_weights(rows)
    if ('model_source' %in% names(rows) && n_distinct(rows$model_source) > 1L) {
        source_n <- table(rows$model_source)
        weights <- weights / as.numeric(source_n[rows$model_source])
        weights <- weights / mean(weights)
    }
    set.seed(seed)
    gbm(
        as.formula(paste(response, '~', paste(predictors, collapse = '+'))),
        data = rows, distribution = distribution,
        weights = weights,
        n.trees = 800L, interaction.depth = depth,
        shrinkage = 0.025, n.minobsinnode = 10L,
        bag.fraction = 0.7, train.fraction = 1,
        keep.data = FALSE, verbose = FALSE
    )
}

fit_process <- function(rows, process, predictors, seed, depth = 3L) {
    occurrence_response <- paste0('has_', process)
    fraction_response <- paste0(process, '_fraction')
    positive <- rows[[occurrence_response]] == 1
    magnitude_rows <- rows[positive, , drop = FALSE]
    magnitude_rows$magnitude_logit <- qlogis(pmin(pmax(
        magnitude_rows[[fraction_response]], 1e-4
    ), 1 - 1e-4))
    occurrence <- fit_component(
        rows, occurrence_response, predictors, 'bernoulli', seed, depth
    )
    magnitude <- fit_component(
        magnitude_rows, 'magnitude_logit', predictors, 'gaussian',
        seed + 1L, depth
    )
    fitted_logit <- predict(
        magnitude, magnitude_rows, n.trees = 800L, type = 'response'
    )
    residual_quantiles <- as.numeric(quantile(
        magnitude_rows$magnitude_logit - fitted_logit,
        probs = seq(0.01, 0.99, length.out = 99L), na.rm = TRUE
    ))
    list(
        occurrence = occurrence, magnitude = magnitude,
        residual_quantiles = residual_quantiles,
        predictors = predictors, process = process
    )
}

predict_process <- function(model, new_data) {
    probability <- pmin(pmax(predict(
        model$occurrence, new_data,
        n.trees = 800L, type = 'response'
    ), 0), 1)
    magnitude_logit <- predict(
        model$magnitude, new_data,
        n.trees = 800L, type = 'response'
    )
    # Duan-style smearing on the logit scale avoids treating inverse-logit of
    # the conditional mean as the expected bounded magnitude.
    magnitude <- rowMeans(plogis(outer(
        magnitude_logit, model$residual_quantiles, '+'
    )))
    tibble(probability = probability, magnitude = magnitude)
}

fit_decomposition <- function(analysis, rich_loss, seed,
                              auxiliary_loss = NULL) {
    loss_predictors <- if (rich_loss) {
        loss_modifier_predictors
    } else {
        loss_core_predictors
    }
    loss_rows <- if (is.null(auxiliary_loss)) {
        analysis
    } else {
        bind_rows(analysis, auxiliary_loss)
    }
    list(
        gain = fit_process(
            analysis, 'gain', growth_predictors, seed, depth = 2L
        ),
        loss = fit_process(
            loss_rows, 'loss', loss_predictors, seed + 10L,
            depth = if (rich_loss) 3L else 2L
        ),
        rich_loss = rich_loss
    )
}

predict_decomposition <- function(model, new_data) {
    gain <- predict_process(model$gain, new_data)
    loss <- predict_process(model$loss, new_data)
    total_probability <- gain$probability + loss$probability
    normaliser <- pmax(total_probability, 1)
    gain_probability <- gain$probability / normaliser
    loss_probability <- loss$probability / normaliser
    expected_gain <- new_data$available_space *
        gain_probability * gain$magnitude
    expected_loss <- new_data$pre_cover *
        loss_probability * loss$magnitude
    predicted_post <- pmin(pmax(
        new_data$pre_cover + expected_gain - expected_loss, 0
    ), 1)
    tibble(
        gain_probability = gain_probability,
        loss_probability = loss_probability,
        conditional_gain_fraction = gain$magnitude,
        conditional_loss_fraction = loss$magnitude,
        expected_gain_pp = 100 * expected_gain,
        expected_loss_pp = 100 * expected_loss,
        predicted_post_cover = predicted_post,
        predicted_change_pp = 100 * (predicted_post - new_data$pre_cover)
    )
}

schemes <- c(
    'recent_year_holdout', 'reef_blocked_5fold', 'forward_2024_2025'
)
all_predictions <- tibble()
all_preprocessing <- tibble()

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
        prepared <- prepare_fold(
            data[analysis_index, , drop = FALSE],
            data[assessment_index, , drop = FALSE],
            mortality_auxiliary[auxiliary_index, , drop = FALSE]
        )
        seed <- 202415L + match(scheme, schemes) * 100L +
            match(fold, folds)

        core <- fit_decomposition(prepared$analysis, FALSE, seed)
        rich <- fit_decomposition(prepared$analysis, TRUE, seed + 30L)
        augmented <- fit_decomposition(
            prepared$analysis, TRUE, seed + 60L,
            auxiliary_loss = prepared$auxiliary
        )
        fold_predictions <- bind_rows(
            predict_decomposition(core, prepared$assessment) |>
                mutate(learner = 'Core decomposed BRT'),
            predict_decomposition(rich, prepared$assessment) |>
                mutate(learner = 'Modifier-rich decomposed BRT'),
            predict_decomposition(augmented, prepared$assessment) |>
                mutate(learner = 'Mortality-augmented decomposed BRT')
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
                        wqc_freqcc12, wqc_prior10_percentile,
                        log_coastal_rain30, era5_wind_calm_fraction,
                        cyc_maxHrs4mw, cot_idwmeanpertow,
                        mortality_environment_available, weather_available,
                        DISTURBANCE_TYPE, disturbance_text
                    ),
                by = 'row_index', relationship = 'many-to-one'
            ) |>
            mutate(
                post_cover_error = post_cover - predicted_post_cover,
                change_error_pp = cover_change_pp - predicted_change_pp,
                scheme = scheme, fold = as.character(fold)
            ) |>
            select(-row_index)

        all_predictions <- bind_rows(all_predictions, fold_predictions)
        all_preprocessing <- bind_rows(
            all_preprocessing,
            prepared$rules |>
                mutate(scheme = scheme, fold = as.character(fold))
        )
        cat('Completed decomposed annual fold:', scheme, fold, '\n')
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
            change_bias_pp = mean(predicted_change_pp - cover_change_pp),
            change_predictive_r2 = 1 -
                sum(change_error_pp^2) /
                sum((cover_change_pp - mean(cover_change_pp))^2),
            direction_accuracy = mean(
                sign(predicted_change_pp) == sign(cover_change_pp)
            ),
            severe_loss_n = sum(cover_change_pp <= -10),
            severe_loss_mae_pp = mean(
                abs(change_error_pp)[cover_change_pp <= -10]
            ),
            severe_loss_bias_pp = mean(
                predicted_change_pp[cover_change_pp <= -10] -
                    cover_change_pp[cover_change_pp <= -10]
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
    file.path(output_dir, 'decomposed_change_predictions.csv')
)
write_csv(metrics, file.path(output_dir, 'decomposed_change_metrics.csv'))
write_csv(
    all_preprocessing,
    file.path(output_dir, 'decomposed_change_preprocessing.csv')
)

# Matched comparison to the previous annual benchmarks.
benchmarks <- read_csv(
    'output/annual_coral_change/annual_change_predictions.csv',
    show_col_types = FALSE
) |>
    select(
        learner, scheme, fold, transition_id, programme_key,
        ReefID, ReefName, event_year, pre_cover, post_cover,
        cover_change_pp, predicted_post_cover, predicted_change_pp,
        post_cover_error, change_error_pp
    )
comparison_rows <- bind_rows(
    benchmarks,
    all_predictions |>
        select(any_of(names(benchmarks)))
)
comparison_metrics <- bind_rows(
    metric_summary(comparison_rows) |> mutate(programme_key = 'all'),
    metric_summary(comparison_rows, by_programme = TRUE)
)
write_csv(
    comparison_metrics,
    file.path(output_dir, 'matched_model_comparison_metrics.csv')
)

# Refit on all transitions for interpretation and subsequent operational use.
prepared_all <- prepare_fold(data, data)
prepared_auxiliary <- apply_preprocessing(
    mortality_auxiliary, prepared_all$rules
)
production <- fit_decomposition(
    prepared_all$analysis, rich_loss = TRUE, seed = 202415L,
    auxiliary_loss = prepared_auxiliary
)
saveRDS(
    list(
        model = production,
        preprocessing = prepared_all$rules,
        growth_predictors = growth_predictors,
        loss_predictors = loss_modifier_predictors,
        n_trees = 800L,
        model_version = 'decomposed_mortality_augmented_brt_v2'
    ),
    file.path(output_dir, 'decomposed_change_production_model.rds')
)

importance <- bind_rows(
    summary(production$gain$occurrence, plotit = FALSE) |>
        transmute(
            component = 'gain occurrence', variable = var,
            influence = rel.inf
        ),
    summary(production$gain$magnitude, plotit = FALSE) |>
        transmute(
            component = 'gain magnitude', variable = var,
            influence = rel.inf
        ),
    summary(production$loss$occurrence, plotit = FALSE) |>
        transmute(
            component = 'loss occurrence', variable = var,
            influence = rel.inf
        ),
    summary(production$loss$magnitude, plotit = FALSE) |>
        transmute(
            component = 'loss magnitude', variable = var,
            influence = rel.inf
        )
)
write_csv(
    importance,
    file.path(output_dir, 'decomposed_component_importance.csv')
)

outliers <- all_predictions |>
    filter(
        scheme == 'forward_2024_2025',
        learner == 'Mortality-augmented decomposed BRT'
    ) |>
    mutate(absolute_change_error_pp = abs(change_error_pp)) |>
    arrange(desc(absolute_change_error_pp))
write_csv(
    outliers,
    file.path(output_dir, 'mortality_augmented_forward_outliers.csv')
)

event_summary <- all_predictions |>
    group_by(learner, scheme, event_year) |>
    summarise(
        transitions = n(),
        observed_change_pp = mean(cover_change_pp),
        severe_loss_n = sum(cover_change_pp <= -10),
        severe_loss_observed_pp = mean(
            cover_change_pp[cover_change_pp <= -10]
        ),
        severe_loss_predicted_pp = mean(
            predicted_change_pp[cover_change_pp <= -10]
        ),
        predicted_change_pp = mean(predicted_change_pp),
        predicted_gain_pp = mean(expected_gain_pp),
        predicted_loss_pp = mean(expected_loss_pp),
        rmse_pp = sqrt(mean(change_error_pp^2)),
        .groups = 'drop'
    )
write_csv(
    event_summary,
    file.path(output_dir, 'decomposed_event_summary.csv')
)

print(filter(comparison_metrics, programme_key == 'all'))
