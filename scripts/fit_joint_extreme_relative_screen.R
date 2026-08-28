# Reef-blocked screen of a coherent Acropora parameterisation and an
# extreme-aware relative-mortality ensemble.
#
# Zeros remain in the Bernoulli occurrence model. Positive loss is modelled on
# the logit scale. A separate >=30% mortality gate and tail-magnitude learner
# test whether severe cases can be improved without broadly inflating ordinary
# predictions. All environmental preprocessing is learned inside each analysis
# fold, and every reef is confined to one side of a fold.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
    library(tidyr)
})
source('scripts/joint_compound_model_helpers.R')

set.seed(20260826L)
output_dir <- 'output/joint_extreme_relative'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tail_threshold <- 0.30
n_trees <- 900L

data <- load_joint_compound_rows() |>
    mutate(
        acropora_cover_pre = observed_pre_cover * prop_acropora_pre,
        region_factor = factor(
            region_block,
            levels = c('Northern GBR', 'Central GBR', 'Southern GBR')
        )
    )

raw_predictors <- c(
    'ann_maxdhw', 'dhw_excess4', 'dhw_excess8',
    'prop_acropora_pre', 'acropora_cover_pre', 'observed_pre_cover',
    'histmDHW6', 'dhw10_load4', 'dhw_novelty10',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8',
    'secc3m_p10', 'cloudp_90',
    'log_coastal_rain30', 'wqc_prior10_percentile',
    'log1p_cyc_maxHrs4mw', 'log1p_cot_idwmeanpertow'
)

shared_terms <- c(
    'programme_factor', 'region_factor',
    'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z',
    'observed_pre_cover_z', 'histmDHW6_z', 'dhw10_load4_z',
    'dhw_novelty10_z', 'dhw_events_since2016_n6_z',
    'dhw_years_since_last_n6_capped8_z',
    'secc3m_p10_z', 'cloudp_90_z',
    'log_coastal_rain30_z', 'wqc_prior10_percentile_z',
    'log1p_cyc_maxHrs4mw_z', 'log1p_cot_idwmeanpertow_z',
    'dhw_x_novelty', 'dhw_x_freshwater', 'dhw_x_cloud'
)

composition_variants <- list(
    proportion_only = c(
        'prop_acropora_pre_z', 'dhw_x_acropora_proportion'
    ),
    absolute_cover_only = c(
        'acropora_cover_pre_z', 'dhw_x_acropora_absolute'
    ),
    legacy_dual = c(
        'prop_acropora_pre_z', 'acropora_cover_pre_z',
        'dhw_x_acropora_proportion'
    )
)

add_interactions <- function(rows) {
    rows |>
        mutate(
            freshwater_joint_z = (
                log_coastal_rain30_z + wqc_prior10_percentile_z
            ) / sqrt(2),
            dhw_x_acropora_proportion =
                dhw_excess4_z * prop_acropora_pre_z,
            dhw_x_acropora_absolute =
                dhw_excess4_z * acropora_cover_pre_z,
            dhw_x_novelty = dhw_excess4_z * dhw_novelty10_z,
            dhw_x_freshwater = dhw_excess4_z * freshwater_joint_z,
            dhw_x_cloud = dhw_excess4_z * cloudp_90_z,
            tail_loss = as.numeric(mortality_prop >= tail_threshold)
        )
}

prepare_screen_fold <- function(analysis, assessment) {
    prepared <- prepare_fold_predictors(
        analysis, assessment, raw_predictors
    )
    list(
        analysis = add_interactions(prepared$analysis),
        assessment = add_interactions(prepared$assessment),
        preprocessing = prepared$preprocessing
    )
}

reef_event_weights <- function(rows) {
    event_size <- table(rows$reef_event_effect)
    programme_size <- table(rows$programme_key)
    weight <- 1 / as.numeric(event_size[rows$reef_event_effect])
    weight <- weight / as.numeric(programme_size[rows$programme_key])
    weight / mean(weight)
}

fit_component <- function(rows, response, distribution, predictors,
                          weights, depth, minobs, seed) {
    outcome <- rows[[response]]
    if (length(unique(outcome)) < 2L || sd(outcome) == 0) {
        return(list(type = 'constant', value = weighted.mean(outcome, weights)))
    }
    monotone <- ifelse(
        predictors %in% c(
            'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z',
            'prop_acropora_pre_z', 'acropora_cover_pre_z'
        ),
        1L, 0L
    )
    set.seed(seed)
    model <- do.call(gbm, list(
        formula = as.formula(paste(
            response, '~', paste(predictors, collapse = ' + ')
        )),
        data = rows, distribution = distribution,
        weights = as.numeric(weights), var.monotone = monotone,
        n.trees = n_trees, interaction.depth = depth,
        shrinkage = 0.025,
        n.minobsinnode = min(minobs, max(1L, floor(nrow(rows) / 6L))),
        bag.fraction = 0.70, train.fraction = 1,
        keep.data = FALSE, verbose = FALSE
    ))
    list(type = 'gbm', model = model)
}

predict_component <- function(component, newdata, type = 'response') {
    if (component$type == 'constant') {
        return(rep(component$value, nrow(newdata)))
    }
    predict(
        component$model, newdata = newdata,
        n.trees = n_trees, type = type
    )
}

fit_base_model <- function(rows, composition, seed) {
    predictors <- c(shared_terms, composition_variants[[composition]])
    weights <- reef_event_weights(rows)
    occurrence_rows <- rows |>
        mutate(has_loss = as.numeric(mortality_prop > 0))
    positive <- occurrence_rows |> filter(has_loss == 1)
    positive_weights <- weights[occurrence_rows$has_loss == 1]
    epsilon <- min(0.01, 0.5 / nrow(positive))
    positive$positive_logit <- qlogis(pmin(pmax(
        positive$mortality_prop, epsilon
    ), 1 - epsilon))
    list(
        occurrence = fit_component(
            occurrence_rows, 'has_loss', 'bernoulli', predictors,
            weights, depth = 3L, minobs = 8L, seed = seed
        ),
        magnitude = fit_component(
            positive, 'positive_logit', 'gaussian', predictors,
            positive_weights, depth = 3L, minobs = 5L, seed = seed + 1L
        ),
        predictors = predictors
    )
}

predict_base_model <- function(model, newdata) {
    occurrence <- pmin(pmax(predict_component(
        model$occurrence, newdata, 'response'
    ), 0), 1)
    positive <- plogis(predict_component(
        model$magnitude, newdata, 'response'
    ))
    tibble(
        predicted_occurrence = occurrence,
        predicted_positive_mortality = positive,
        predicted_mortality = occurrence * positive
    )
}

fit_tail_model <- function(rows, seed) {
    predictors <- c(
        shared_terms,
        composition_variants$proportion_only,
        'freshwater_joint_z'
    )
    weights <- reef_event_weights(rows)
    tail_rows <- rows |> filter(tail_loss == 1)
    tail_weights <- weights[rows$tail_loss == 1]
    epsilon <- min(0.01, 0.5 / nrow(tail_rows))
    tail_rows$tail_logit <- qlogis(pmin(pmax(
        tail_rows$mortality_prop, epsilon
    ), 1 - epsilon))
    list(
        gate = fit_component(
            rows, 'tail_loss', 'bernoulli', predictors,
            weights, depth = 3L, minobs = 8L, seed = seed
        ),
        magnitude = fit_component(
            tail_rows, 'tail_logit', 'laplace', predictors,
            tail_weights, depth = 2L, minobs = 3L, seed = seed + 1L
        ),
        predictors = predictors,
        tail_rows = nrow(tail_rows)
    )
}

predict_tail_model <- function(model, newdata) {
    probability <- pmin(pmax(predict_component(
        model$gate, newdata, 'response'
    ), 0), 1)
    magnitude <- plogis(predict_component(
        model$magnitude, newdata, 'response'
    ))
    magnitude <- pmin(pmax(magnitude, tail_threshold), 1)
    tibble(
        tail_probability = probability,
        tail_magnitude = magnitude
    )
}

composition_predictions <- tibble()
ensemble_predictions <- tibble()
runtime <- tibble()

for (fold in sort(unique(data$joint_reef_fold))) {
    assessment_index <- data$joint_reef_fold == fold
    held_out_reefs <- unique(data$ReefID[assessment_index])
    analysis <- data[!data$ReefID %in% held_out_reefs, , drop = FALSE]
    assessment <- data[assessment_index, , drop = FALSE]
    prepared <- prepare_screen_fold(analysis, assessment)

    selected_base <- NULL
    selected_prediction <- NULL
    for (composition in names(composition_variants)) {
        started <- proc.time()[['elapsed']]
        model <- fit_base_model(
            prepared$analysis, composition,
            seed = 20260826L + fold * 100L +
                match(composition, names(composition_variants)) * 10L
        )
        prediction <- predict_base_model(model, prepared$assessment)
        elapsed <- proc.time()[['elapsed']] - started
        composition_predictions <- bind_rows(
            composition_predictions,
            assessment |>
                transmute(
                    programme_key, source_observation_id, ReefID, ReefName,
                    event_year, ann_maxdhw, prop_acropora_pre,
                    acropora_cover_pre, observed_mortality = mortality_prop,
                    fold = as.character(.env$fold),
                    composition = .env$composition
                ) |>
                bind_cols(prediction)
        )
        runtime <- bind_rows(runtime, tibble(
            fold, stage = 'composition', candidate = composition,
            elapsed_seconds = elapsed
        ))
        if (composition == 'proportion_only') {
            selected_base <- model
            selected_prediction <- prediction
        }
    }

    started <- proc.time()[['elapsed']]
    tail_model <- fit_tail_model(
        prepared$analysis, seed = 20261826L + fold * 100L
    )
    tail_prediction <- predict_tail_model(
        tail_model, prepared$assessment
    )
    elapsed <- proc.time()[['elapsed']] - started
    runtime <- bind_rows(runtime, tibble(
        fold, stage = 'tail', candidate = 'tail_gate_and_magnitude',
        elapsed_seconds = elapsed
    ))

    base_rows <- assessment |>
        transmute(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, region_block, ann_maxdhw, prop_acropora_pre,
            observed_pre_cover, dhw_novelty10, cloudp_90,
            log_coastal_rain30, wqc_prior10_percentile,
            observed_mortality = mortality_prop,
            fold = as.character(.env$fold)
        ) |>
        bind_cols(selected_prediction, tail_prediction)

    ensemble_predictions <- bind_rows(
        ensemble_predictions,
        base_rows |>
            mutate(
                learner = 'Standard two-part BRT',
                ensemble_prediction = predicted_mortality
            ),
        base_rows |>
            mutate(
                learner = 'Conservative tail ensemble',
                ensemble_prediction = predicted_mortality +
                    0.5 * tail_probability *
                    (tail_magnitude - predicted_mortality)
            ),
        base_rows |>
            mutate(
                learner = 'Soft-gated tail ensemble',
                ensemble_prediction = predicted_mortality +
                    tail_probability *
                    (tail_magnitude - predicted_mortality)
            ),
        base_rows |>
            mutate(
                learner = 'Aggressive tail ensemble',
                ensemble_prediction = predicted_mortality +
                    sqrt(tail_probability) *
                    (tail_magnitude - predicted_mortality)
            )
    )
    message('Completed joint extreme fold ', fold)
}

metric_summary <- function(rows) {
    rows |>
        summarise(
            n = n(), severe_n = sum(observed_mortality >= 0.50),
            observed_mean = mean(observed_mortality),
            predicted_mean = mean(ensemble_prediction),
            rmse = sqrt(mean((observed_mortality - ensemble_prediction)^2)),
            mae = mean(abs(observed_mortality - ensemble_prediction)),
            predictive_r2 = 1 -
                sum((observed_mortality - ensemble_prediction)^2) /
                sum((observed_mortality - mean(observed_mortality))^2),
            bias = mean(ensemble_prediction - observed_mortality),
            severe_observed = if_else(
                severe_n > 0,
                mean(observed_mortality[observed_mortality >= 0.50]), NA_real_
            ),
            severe_predicted = if_else(
                severe_n > 0,
                mean(ensemble_prediction[observed_mortality >= 0.50]), NA_real_
            ),
            severe_bias = severe_predicted - severe_observed,
            severe_rmse = if_else(
                severe_n > 0,
                sqrt(mean((
                    observed_mortality[observed_mortality >= 0.50] -
                        ensemble_prediction[observed_mortality >= 0.50]
                )^2)), NA_real_
            ),
            false_extreme_rate = mean(
                ensemble_prediction[observed_mortality < tail_threshold] >=
                    tail_threshold
            ),
            .groups = 'drop'
        )
}

composition_metrics <- composition_predictions |>
    group_by(programme_key, composition) |>
    summarise(
        n = n(),
        rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
        mae = mean(abs(observed_mortality - predicted_mortality)),
        bias = mean(predicted_mortality - observed_mortality),
        severe_n = sum(observed_mortality >= 0.50),
        severe_predicted = if_else(
            severe_n > 0,
            mean(predicted_mortality[observed_mortality >= 0.50]), NA_real_
        ),
        severe_observed = if_else(
            severe_n > 0,
            mean(observed_mortality[observed_mortality >= 0.50]), NA_real_
        ),
        .groups = 'drop'
    )

ensemble_metrics <- ensemble_predictions |>
    group_by(programme_key, learner) |>
    metric_summary()
ensemble_metrics_2024 <- ensemble_predictions |>
    filter(event_year == 2024) |>
    group_by(programme_key, learner) |>
    metric_summary()
event_metrics <- ensemble_predictions |>
    group_by(programme_key, learner, event_year) |>
    metric_summary()

# Production fit and joint low/median/high compound scenarios.
production_prepared <- prepare_screen_fold(data, data[0, , drop = FALSE])
production_base <- fit_base_model(
    production_prepared$analysis, 'proportion_only', 20260826L
)
production_tail <- fit_tail_model(
    production_prepared$analysis, 20261826L
)

scenario_levels <- c('Low compound risk', 'Median conditions',
                     'High compound risk')
scenario_probabilities <- c(0.10, 0.50, 0.90)
scenario_grid <- expand_grid(
    ann_maxdhw = seq(0, 16, by = 0.25),
    scenario = factor(scenario_levels, levels = scenario_levels),
    programme_key = joint_programmes
)
raw_grid <- data[rep(1, nrow(scenario_grid)), , drop = FALSE]
for (predictor in raw_predictors) {
    raw_grid[[predictor]] <- ave(
        data[[predictor]], data$programme_key,
        FUN = function(x) median(x, na.rm = TRUE)
    )[match(scenario_grid$programme_key, data$programme_key)]
}
raw_grid$programme_key <- scenario_grid$programme_key
raw_grid$programme_factor <- factor(
    raw_grid$programme_key, levels = joint_programmes
)
modal_region <- data |>
    count(programme_key, region_factor, name = 'n') |>
    group_by(programme_key) |>
    slice_max(n, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(programme_key, region_factor)
raw_grid$region_factor <- modal_region$region_factor[
    match(raw_grid$programme_key, modal_region$programme_key)
]
raw_grid$ann_maxdhw <- scenario_grid$ann_maxdhw
raw_grid$dhw_excess4 <- pmax(raw_grid$ann_maxdhw - 4, 0)
raw_grid$dhw_excess8 <- pmax(raw_grid$ann_maxdhw - 8, 0)

set_scenario_quantile <- function(rows, predictor, reverse = FALSE) {
    for (programme in joint_programmes) {
        source <- data[[predictor]][data$programme_key == programme]
        values <- quantile(
            source[is.finite(source)], scenario_probabilities,
            names = FALSE
        )
        if (reverse) values <- rev(values)
        for (i in seq_along(scenario_levels)) {
            use <- scenario_grid$programme_key == programme &
                scenario_grid$scenario == scenario_levels[[i]]
            rows[[predictor]][use] <- values[[i]]
        }
    }
    rows
}
raw_grid <- set_scenario_quantile(raw_grid, 'prop_acropora_pre')
raw_grid <- set_scenario_quantile(raw_grid, 'dhw_novelty10')
raw_grid <- set_scenario_quantile(raw_grid, 'log_coastal_rain30')
raw_grid <- set_scenario_quantile(raw_grid, 'wqc_prior10_percentile')
raw_grid <- set_scenario_quantile(raw_grid, 'cloudp_90', reverse = TRUE)
raw_grid$acropora_cover_pre <-
    raw_grid$observed_pre_cover * raw_grid$prop_acropora_pre

prepared_grid <- prepare_fold_predictors(data, raw_grid, raw_predictors)
profile_rows <- add_interactions(prepared_grid$assessment)
base_curve <- predict_base_model(production_base, profile_rows)
tail_curve <- predict_tail_model(production_tail, profile_rows)
dhw_support <- data |>
    group_by(programme_key) |>
    summarise(
        observed_dhw_min = min(ann_maxdhw, na.rm = TRUE),
        observed_dhw_max = max(ann_maxdhw, na.rm = TRUE),
        .groups = 'drop'
    )

scenario_curves <- bind_cols(
    scenario_grid,
    profile_rows |>
        select(
            prop_acropora_pre, dhw_novelty10, log_coastal_rain30,
            wqc_prior10_percentile, cloudp_90
        ),
    base_curve, tail_curve
) |>
    mutate(
        conservative_prediction = predicted_mortality +
            0.5 * tail_probability *
            (tail_magnitude - predicted_mortality),
        soft_gate_prediction = predicted_mortality +
            tail_probability *
            (tail_magnitude - predicted_mortality),
        aggressive_prediction = predicted_mortality +
            sqrt(tail_probability) *
            (tail_magnitude - predicted_mortality)
    ) |>
    left_join(dhw_support, by = 'programme_key') |>
    mutate(
        outside_observed_dhw = ann_maxdhw < observed_dhw_min |
            ann_maxdhw > observed_dhw_max
    )

production_influence <- bind_rows(lapply(list(
    base_occurrence = production_base$occurrence,
    base_magnitude = production_base$magnitude,
    tail_gate = production_tail$gate,
    tail_magnitude = production_tail$magnitude
), function(component) {
    if (component$type != 'gbm') return(tibble())
    summary(component$model, plotit = FALSE) |>
        as_tibble() |>
        transmute(predictor = var, relative_influence = rel.inf)
}), .id = 'component')

stopifnot(
    nrow(composition_predictions) == nrow(data) * length(composition_variants),
    nrow(ensemble_predictions) == nrow(data) * 4L,
    all(is.finite(composition_predictions$predicted_mortality)),
    all(composition_predictions$predicted_mortality >= 0 &
            composition_predictions$predicted_mortality <= 1),
    all(is.finite(ensemble_predictions$ensemble_prediction)),
    all(ensemble_predictions$ensemble_prediction >= 0 &
            ensemble_predictions$ensemble_prediction <= 1),
    all((data |>
             distinct(ReefID, joint_reef_fold) |>
             count(ReefID) |>
             pull(n)) == 1L)
)

write_csv(composition_predictions, file.path(output_dir, 'composition_predictions.csv'), na = '')
write_csv(composition_metrics, file.path(output_dir, 'composition_metrics.csv'), na = '')
write_csv(ensemble_predictions, file.path(output_dir, 'ensemble_predictions.csv'), na = '')
write_csv(ensemble_metrics, file.path(output_dir, 'ensemble_metrics.csv'), na = '')
write_csv(ensemble_metrics_2024, file.path(output_dir, 'ensemble_metrics_2024.csv'), na = '')
write_csv(event_metrics, file.path(output_dir, 'event_metrics.csv'), na = '')
write_csv(scenario_curves, file.path(output_dir, 'joint_scenario_curves.csv'), na = '')
write_csv(production_influence, file.path(output_dir, 'production_influence.csv'), na = '')
write_csv(runtime, file.path(output_dir, 'runtime.csv'), na = '')

print(composition_metrics)
print(ensemble_metrics)
print(ensemble_metrics_2024)
