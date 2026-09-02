# Distributionally appropriate annual transition benchmark. Gain and loss
# occurrence use binomial GAMs; conditional magnitudes use beta GAMs. The loss
# process retains the relative-mortality modifier contract and restricted
# bleaching-mortality observations as a source-balanced auxiliary response.

suppressPackageStartupMessages({
    library(dplyr)
    library(mgcv)
    library(readr)
})
source('src/lib/joint_compound_model_helpers.R')

set.seed(202417L)
output_dir <- 'output/beta_binomial_annual'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data <- readRDS('data/processed/annual_coral_transitions.rds') |>
    mutate(
        available_space = pmax(1 - pre_cover, 0),
        has_gain = as.numeric(cover_change_pp > 0),
        has_loss = as.numeric(cover_change_pp < 0),
        gain_fraction = pmin(pmax(
            (post_cover - pre_cover) / pmax(available_space, 1e-6), 0
        ), 1),
        loss_fraction = pmin(pmax(
            (pre_cover - post_cover) / pmax(pre_cover, 1e-6), 0
        ), 1),
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
    transmute(
        programme_key, ReefID, event_year, region_block,
        pre_cover = observed_pre_cover,
        available_space = pmax(1 - observed_pre_cover, 0),
        interval_years = interval_days / 365.25,
        prior_change_pp = NA_real_, prop_acropora_pre,
        acropora_cover_pre = observed_pre_cover * prop_acropora_pre,
        depth, sst_maxdhw = ann_maxdhw,
        histmDHW6, yrsince6, dhw10_load4, dhw_novelty10,
        dhw_events_since2016_n6, dhw_years_since_last_n6_capped8,
        dhw_no_prior_n6,
        secc3m_p10, cloudp_90, wqc_freqcc12,
        wqc_prior10_percentile, log1p_cyc_maxHrs4mw,
        log1p_cot_idwmeanpertow, log_coastal_rain30,
        era5_wind_mean, era5_wind_calm_fraction,
        era5_coastal_distance_km,
        previous_change_available = 0,
        wqc_available = as.numeric(is.finite(wqc_freqcc12)),
        cots_available = as.numeric(is.finite(log1p_cot_idwmeanpertow)),
        mortality_environment_available = 1,
        weather_available = 1,
        has_loss = as.numeric(mortality_prop > 0),
        loss_fraction = mortality_prop,
        model_source = 'restricted bleaching mortality'
    )

numeric_predictors <- c(
    'pre_cover', 'available_space', 'interval_years', 'prior_change_pp',
    'prop_acropora_pre', 'acropora_cover_pre', 'depth', 'sst_maxdhw',
    'histmDHW6', 'yrsince6', 'dhw10_load4', 'dhw_novelty10',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8',
    'dhw_no_prior_n6',
    'secc3m_p10', 'cloudp_90', 'wqc_freqcc12',
    'wqc_prior10_percentile', 'log1p_cyc_maxHrs4mw',
    'log1p_cot_idwmeanpertow', 'log_coastal_rain30',
    'era5_wind_mean', 'era5_wind_calm_fraction',
    'era5_coastal_distance_km', 'uq_coralsink1',
    'previous_change_available', 'wqc_available', 'cots_available',
    'mortality_environment_available', 'weather_available'
)

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
        if (!predictor %in% names(rows)) rows[[predictor]] <- NA_real_
        rows[[predictor]][!is.finite(rows[[predictor]])] <- rules$median[[i]]
    }
    rows |>
        mutate(
            programme_factor = factor(
                programme_key, levels = c('ltmp', 'manta', 'mmp')
            ),
            region_factor = factor(
                region_block,
                levels = c('Northern GBR', 'Central GBR', 'Southern GBR')
            )
        )
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

growth_rhs <- paste(c(
    'programme_factor', 'region_factor',
    's(pre_cover, k = 5)', 's(interval_years, k = 4)',
    's(prior_change_pp, k = 5)', 's(prop_acropora_pre, k = 5)',
    's(sst_maxdhw, k = 6)', 's(dhw10_load4, k = 5)',
    's(dhw_novelty10, k = 5)', 's(log1p_cyc_maxHrs4mw, k = 5)',
    's(log1p_cot_idwmeanpertow, k = 5)', 's(uq_coralsink1, k = 5)',
    'previous_change_available', 'mortality_environment_available'
), collapse = ' + ')

loss_rhs <- paste(c(
    'programme_factor', 'region_factor',
    's(pre_cover, k = 5)', 's(interval_years, k = 4)',
    's(prior_change_pp, k = 5)', 's(prop_acropora_pre, k = 5)',
    's(sst_maxdhw, k = 6)',
    'ti(sst_maxdhw, prop_acropora_pre, k = c(5, 4))',
    's(histmDHW6, k = 5)', 's(yrsince6, k = 5)',
    's(dhw10_load4, k = 5)', 's(dhw_novelty10, k = 5)',
    's(dhw_events_since2016_n6, k = 4)',
    's(dhw_years_since_last_n6_capped8, k = 5)', 'dhw_no_prior_n6',
    'ti(sst_maxdhw, dhw_events_since2016_n6, k = c(5, 4))',
    'ti(sst_maxdhw, dhw_years_since_last_n6_capped8, k = c(5, 4))',
    's(secc3m_p10, k = 5)', 's(cloudp_90, k = 5)',
    's(wqc_freqcc12, k = 5)',
    'ti(sst_maxdhw, wqc_freqcc12, k = c(5, 4))',
    's(log_coastal_rain30, k = 5)',
    'ti(sst_maxdhw, log_coastal_rain30, k = c(5, 4))',
    's(era5_wind_calm_fraction, k = 5)',
    'ti(sst_maxdhw, era5_wind_calm_fraction, k = c(5, 4))',
    's(log1p_cyc_maxHrs4mw, k = 5)',
    's(log1p_cot_idwmeanpertow, k = 5)',
    'wqc_available', 'cots_available',
    'mortality_environment_available', 'weather_available'
), collapse = ' + ')

fit_process <- function(rows, process, rhs) {
    occurrence_name <- paste0('has_', process)
    fraction_name <- paste0(process, '_fraction')
    occurrence_formula <- as.formula(paste(occurrence_name, '~', rhs))
    occurrence <- bam(
        occurrence_formula, data = rows, family = binomial(),
        weights = model_weights(rows), method = 'fREML', discrete = TRUE
    )
    positive <- rows[[occurrence_name]] == 1
    magnitude_rows <- rows[positive, , drop = FALSE]
    magnitude_rows$bounded_fraction <- pmin(pmax(
        magnitude_rows[[fraction_name]], 1e-4
    ), 1 - 1e-4)
    magnitude <- bam(
        as.formula(paste('bounded_fraction ~', rhs)),
        data = magnitude_rows, family = betar(link = 'logit'),
        weights = model_weights(magnitude_rows),
        method = 'fREML', discrete = TRUE
    )
    list(occurrence = occurrence, magnitude = magnitude)
}

fit_decomposition <- function(annual, auxiliary) {
    list(
        gain = fit_process(annual, 'gain', growth_rhs),
        loss = fit_process(bind_rows(annual, auxiliary), 'loss', loss_rhs)
    )
}

predict_process <- function(model, assessment) {
    tibble(
        probability = pmin(pmax(predict(
            model$occurrence, assessment, type = 'response'
        ), 0), 1),
        magnitude = pmin(pmax(predict(
            model$magnitude, assessment, type = 'response'
        ), 0), 1)
    )
}

predict_decomposition <- function(model, assessment) {
    gain <- predict_process(model$gain, assessment)
    loss <- predict_process(model$loss, assessment)
    total <- gain$probability + loss$probability
    normaliser <- pmax(total, 1)
    gain_probability <- gain$probability / normaliser
    loss_probability <- loss$probability / normaliser
    expected_gain <- assessment$available_space *
        gain_probability * gain$magnitude
    expected_loss <- assessment$pre_cover *
        loss_probability * loss$magnitude
    predicted_post <- pmin(pmax(
        assessment$pre_cover + expected_gain - expected_loss, 0
    ), 1)
    tibble(
        gain_probability, loss_probability,
        conditional_gain_fraction = gain$magnitude,
        conditional_loss_fraction = loss$magnitude,
        expected_gain_pp = 100 * expected_gain,
        expected_loss_pp = 100 * expected_loss,
        predicted_post_cover = predicted_post,
        predicted_change_pp = 100 * (predicted_post - assessment$pre_cover)
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
        model <- fit_decomposition(prepared$analysis, prepared$auxiliary)
        prediction <- predict_decomposition(model, prepared$assessment)
        fold_rows <- prepared$assessment |>
            select(
                transition_id, programme_key, ReefID, ReefName, event_year,
                pre_cover, post_cover, cover_change_pp, sst_maxdhw,
                prop_acropora_pre, wqc_freqcc12, log_coastal_rain30,
                era5_wind_calm_fraction, cyc_maxHrs4mw,
                cot_idwmeanpertow, disturbance_text
            ) |>
            bind_cols(prediction) |>
            mutate(
                error_pp = cover_change_pp - predicted_change_pp,
                learner = 'Binomial + beta hurdle GAM',
                scheme = scheme, fold = as.character(fold)
            )
        all_predictions <- bind_rows(all_predictions, fold_rows)
        all_preprocessing <- bind_rows(
            all_preprocessing,
            prepared$rules |>
                mutate(scheme = scheme, fold = as.character(fold))
        )
        cat('Completed beta-binomial fold:', scheme, fold, '\n')
    }
}

metrics <- all_predictions |>
    group_by(learner, scheme) |>
    summarise(
        transitions = n(), severe_n = sum(cover_change_pp <= -10),
        rmse_pp = sqrt(mean(error_pp^2)), mae_pp = mean(abs(error_pp)),
        predictive_r2 = 1 - sum(error_pp^2) /
            sum((cover_change_pp - mean(cover_change_pp))^2),
        severe_mae_pp = mean(abs(error_pp)[cover_change_pp <= -10]),
        severe_bias_pp = mean(
            predicted_change_pp[cover_change_pp <= -10] -
                cover_change_pp[cover_change_pp <= -10]
        ),
        nonsevere_mae_pp = mean(abs(error_pp)[cover_change_pp > -10]),
        .groups = 'drop'
    )

prepared_all <- prepare_fold(data, data, mortality_auxiliary)
production <- fit_decomposition(
    prepared_all$analysis, prepared_all$auxiliary
)
saveRDS(
    list(
        model = production, preprocessing = prepared_all$rules,
        growth_formula = growth_rhs, loss_formula = loss_rhs,
        model_version = 'beta_binomial_gam_v1'
    ),
    file.path(output_dir, 'beta_binomial_production_model.rds')
)

write_csv(
    all_predictions,
    file.path(output_dir, 'beta_binomial_predictions.csv')
)
write_csv(metrics, file.path(output_dir, 'beta_binomial_metrics.csv'))
write_csv(
    all_preprocessing,
    file.path(output_dir, 'beta_binomial_preprocessing.csv')
)
print(metrics)
