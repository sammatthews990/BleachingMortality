# Severe-tail BRT correction for the selected joint INLA mortality model.
#
# The INLA prediction remains the ecological baseline. A separate Bernoulli
# gate estimates Pr(mortality >= 30%) and a robust magnitude learner estimates
# mortality within that tail. WQC threshold/cumulative features and current
# BoM cyclone-track features are tested as explicit ablations. Cyclone-track
# fields are deliberately isolated here because a better Jasper/Kirrily data
# set is expected and should replace this provisional block without changing
# the rest of the model.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(readr)
    library(tidyr)
})
source('src/lib/joint_compound_model_helpers.R')

set.seed(20260827L)
output_dir <- Sys.getenv('INLA_TAIL_OUTPUT_DIR', 'output/inla_severe_tail')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tail_threshold <- 0.30
severe_threshold <- 0.50
n_trees <- 900L
selected_inla <- Sys.getenv(
    'INLA_TAIL_BASELINE', 'persistent_rw1_cyclone_cots_partial_pool'
)
use_lizard_adjusted_dhw <- grepl('lizard_adjusted', selected_inla)
lizard_cluster_reef_ids <- c(
    '14-114', '14-116a', '14-118', '14-123', '14-126', '14-143'
)
lizard_logger_uplift <- read_csv(
    'data/processed/lizard_cluster_logger_dhw_uplift_2024.csv',
    show_col_types = FALSE
) |>
    filter(grepl('FL', series)) |>
    summarise(uplift = mean(dhw_uplift_vs_official_noaa, na.rm = TRUE)) |>
    pull(uplift)

sst_chla <- read_csv(
    'data/processed/sst_chla_features_validation.csv',
    show_col_types = FALSE
) |>
    select(
        ReefID, year, sst_summer_skewness,
        sst_summer_excess_kurtosis, chla_wetseason_median
    )

cyclone <- read_csv(
    'data/processed/bom_cyclone_reef_year.csv', show_col_types = FALSE
) |>
    select(
        ReefID, event_year, tc_min_distance_km,
        tc_nearest_name, tc_nearest_max_wind_ms,
        tc_max_wind_within300km_ms, tc_storms_within300km
    )

cots_hindcast <- read_csv(
    'data/gbrPredsAdj_20262408.csv', show_col_types = FALSE
) |>
    transmute(
        ReefName = reefName, event_year = as.integer(year),
        cots_outbreak_probability = outbrProb
    )

disease <- read_csv(
    'data/processed/noaa_disease_risk_validation.csv',
    show_col_types = FALSE
) |>
    select(
        ReefID, year, disease_risk_applicable,
        disease_risk_max, disease_risk_days_ge1
    )

data <- load_joint_compound_rows() |>
    left_join(
        sst_chla, by = c('ReefID', 'event_year' = 'year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        cyclone, by = c('ReefID', 'event_year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        cots_hindcast, by = c('ReefName', 'event_year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        disease, by = c('ReefID', 'event_year' = 'year'),
        relationship = 'many-to-one'
    ) |>
    mutate(
        ann_maxdhw_original = ann_maxdhw,
        applied_dhw_uplift = if_else(
            use_lizard_adjusted_dhw & event_year == 2024L &
                ReefID %in% lizard_cluster_reef_ids,
            lizard_logger_uplift, 0
        ),
        ann_maxdhw = ann_maxdhw + applied_dhw_uplift,
        dhw_excess4 = pmax(ann_maxdhw - 4, 0),
        dhw_excess8 = pmax(ann_maxdhw - 8, 0),
        dhw_novelty10 = dhw_novelty10 + applied_dhw_uplift,
        region_factor = factor(
            region_block,
            levels = c('Northern GBR', 'Central GBR', 'Southern GBR')
        ),
        acropora_cover_pre = observed_pre_cover * prop_acropora_pre,
        tc_proximity100 = exp(-pmin(tc_min_distance_km, 1000) / 100),
        tc_wind_distance_index = pmax(
            tc_nearest_max_wind_ms - 17, 0
        ) * tc_proximity100,
        tail_loss = as.numeric(mortality_prop >= tail_threshold)
    )

raw_common <- c(
    'ann_maxdhw', 'dhw_excess4', 'dhw_excess8',
    'prop_acropora_pre', 'acropora_cover_pre', 'observed_pre_cover',
    'dhw10_load4', 'dhw_novelty10',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8',
    'secc3m_p10', 'cloudp_90', 'log_coastal_rain30',
    'wqc_prior10_percentile', 'log1p_cyc_maxHrs4mw',
    'log1p_cot_idwmeanpertow', 'cots_outbreak_probability', 'mcur_90',
    'sst_summer_skewness', 'sst_summer_excess_kurtosis',
    'chla_wetseason_median'
)
raw_wqc <- c('wqc_freqcc12', 'wqc_excess50', 'wqc_excess50_10yr_sum')
raw_track <- c('tc_proximity100', 'tc_wind_distance_index')
raw_disease <- c('disease_risk_max', 'disease_risk_days_ge1')
raw_predictors <- unique(c(raw_common, raw_wqc, raw_track, raw_disease))

common_terms <- c(
    'programme_factor', 'region_factor',
    paste0(raw_common, '_z'),
    'dhw_x_acropora', 'dhw_x_novelty', 'dhw_x_freshwater',
    'dhw_x_cloud', 'dhw_x_current'
)
feature_sets <- list(
    legacy_tail = common_terms,
    wqc_threshold_tail = c(
        common_terms, paste0(raw_wqc, '_z'),
        'dhw_x_wqc_excess50', 'dhw_x_wqc_cumulative'
    ),
    wqc_threshold_track_tail = c(
        common_terms, paste0(raw_wqc, '_z'), paste0(raw_track, '_z'),
        'dhw_x_wqc_excess50', 'dhw_x_wqc_cumulative',
        'dhw_x_tc_track'
    ),
    wqc_threshold_track_disease_tail = c(
        common_terms, paste0(raw_wqc, '_z'), paste0(raw_track, '_z'),
        paste0(raw_disease, '_z'),
        'dhw_x_wqc_excess50', 'dhw_x_wqc_cumulative',
        'dhw_x_tc_track', 'dhw_x_disease'
    )
)

add_derived_terms <- function(rows) {
    rows |>
        mutate(
            freshwater_joint_z = (
                log_coastal_rain30_z + wqc_prior10_percentile_z
            ) / sqrt(2),
            dhw_x_acropora = dhw_excess4_z * prop_acropora_pre_z,
            dhw_x_novelty = dhw_excess4_z * dhw_novelty10_z,
            dhw_x_freshwater = dhw_excess4_z * freshwater_joint_z,
            dhw_x_cloud = dhw_excess4_z * cloudp_90_z,
            dhw_x_current = ann_maxdhw_z * mcur_90_z,
            dhw_x_wqc_excess50 = dhw_excess4_z * wqc_excess50_z,
            dhw_x_wqc_cumulative =
                dhw_excess4_z * wqc_excess50_10yr_sum_z,
            dhw_x_tc_track =
                dhw_excess4_z * tc_wind_distance_index_z,
            dhw_x_disease = dhw_excess4_z * disease_risk_max_z
        )
}

prepare_tail_fold <- function(analysis, assessment) {
    prepared <- prepare_fold_predictors(
        analysis, assessment, raw_predictors
    )
    list(
        analysis = add_derived_terms(prepared$analysis),
        assessment = add_derived_terms(prepared$assessment),
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
        return(list(
            type = 'constant', value = weighted.mean(outcome, weights)
        ))
    }
    monotone_positive <- c(
        'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z',
        'prop_acropora_pre_z', 'acropora_cover_pre_z',
        'wqc_excess50_z', 'wqc_excess50_10yr_sum_z',
        'tc_wind_distance_index_z', 'cots_outbreak_probability_z',
        'disease_risk_max_z', 'disease_risk_days_ge1_z'
    )
    monotone <- ifelse(predictors %in% monotone_positive, 1L, 0L)
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

predict_component <- function(component, newdata) {
    if (component$type == 'constant') {
        return(rep(component$value, nrow(newdata)))
    }
    predict(
        component$model, newdata = newdata,
        n.trees = n_trees, type = 'response'
    )
}

fit_tail_model <- function(rows, predictors, seed) {
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
        )
    )
}

predict_tail_model <- function(model, rows) {
    probability <- pmin(pmax(predict_component(model$gate, rows), 0), 1)
    magnitude <- plogis(predict_component(model$magnitude, rows))
    tibble(
        tail_probability = probability,
        tail_magnitude = pmin(pmax(magnitude, tail_threshold), 1)
    )
}

base_predictions <- read_csv(
    'output/inla_spatiotemporal/cv_predictions.csv',
    show_col_types = FALSE
) |>
    filter(candidate == selected_inla) |>
    mutate(source_observation_id = as.character(source_observation_id)) |>
    rename(inla_prediction = predicted_mortality)

keys <- c(
    'programme_key', 'source_observation_id', 'ReefID', 'event_year'
)
fold_map <- base_predictions |>
    select(all_of(keys), scheme, fold, inla_prediction)

cv_predictions <- tibble()
runtime <- tibble()

for (scheme_name in unique(fold_map$scheme)) {
    scheme_map <- fold_map |> filter(scheme == scheme_name)
    for (fold_name in unique(scheme_map$fold)) {
        assessment_map <- scheme_map |> filter(fold == fold_name)
        assessment <- data |>
            inner_join(assessment_map, by = keys, relationship = 'one-to-one')
        analysis <- data |>
            anti_join(assessment_map, by = keys)
        prepared <- prepare_tail_fold(analysis, assessment)

        if (nrow(assessment) != nrow(assessment_map)) {
            stop('INLA assessment join lost rows for ', scheme_name, ' ', fold_name)
        }

        base_rows <- assessment |>
            transmute(
                programme_key, source_observation_id, ReefID, ReefName,
                event_year, region_block, ann_maxdhw, prop_acropora_pre,
                observed_pre_cover, wqc_freqcc12, wqc_excess50,
                wqc_excess50_10yr_sum, wqc_prior10_percentile,
                log1p_cyc_maxHrs4mw, tc_min_distance_km,
                tc_nearest_name, tc_nearest_max_wind_ms,
                tc_wind_distance_index,
                observed_mortality = mortality_prop,
                inla_prediction, scheme = scheme_name,
                fold = as.character(fold_name)
            )

        cv_predictions <- bind_rows(
            cv_predictions,
            base_rows |>
                mutate(
                    tail_features = 'none', learner = 'INLA only',
                    tail_probability = NA_real_, tail_magnitude = NA_real_,
                    predicted_mortality = inla_prediction
                )
        )

        for (feature_name in names(feature_sets)) {
            started <- proc.time()[['elapsed']]
            model <- fit_tail_model(
                prepared$analysis, feature_sets[[feature_name]],
                seed = 20260827L + match(
                    feature_name, names(feature_sets)
                ) * 100L + as.integer(as.factor(fold_name))
            )
            tail <- predict_tail_model(model, prepared$assessment)
            elapsed <- proc.time()[['elapsed']] - started
            runtime <- bind_rows(runtime, tibble(
                scheme = scheme_name, fold = as.character(fold_name),
                tail_features = feature_name, elapsed_seconds = elapsed
            ))

            candidate_rows <- bind_cols(base_rows, tail) |>
                mutate(
                    correction = pmax(
                        tail_magnitude - inla_prediction, 0
                    )
                )
            cv_predictions <- bind_rows(
                cv_predictions,
                candidate_rows |>
                    mutate(
                        tail_features = feature_name,
                        learner = 'Conservative tail (0.5 gate)',
                        predicted_mortality = pmin(
                            inla_prediction +
                                0.5 * tail_probability * correction, 1
                        )
                    ),
                candidate_rows |>
                    mutate(
                        tail_features = feature_name,
                        learner = 'Soft-gated tail',
                        predicted_mortality = pmin(
                            inla_prediction +
                                tail_probability * correction, 1
                        )
                    )
            )
        }
        message('Completed ', scheme_name, ' fold ', fold_name)
    }
}

metric_summary <- function(rows) {
    rows |>
        summarise(
            n = n(), severe_n = sum(observed_mortality >= severe_threshold),
            observed_mean = mean(observed_mortality),
            predicted_mean = mean(predicted_mortality),
            rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
            mae = mean(abs(observed_mortality - predicted_mortality)),
            predictive_r2 = 1 -
                sum((observed_mortality - predicted_mortality)^2) /
                sum((observed_mortality - mean(observed_mortality))^2),
            bias = mean(predicted_mortality - observed_mortality),
            severe_observed = if_else(
                severe_n > 0,
                mean(observed_mortality[observed_mortality >= severe_threshold]),
                NA_real_
            ),
            severe_predicted = if_else(
                severe_n > 0,
                mean(predicted_mortality[observed_mortality >= severe_threshold]),
                NA_real_
            ),
            severe_rmse = if_else(
                severe_n > 0,
                sqrt(mean((
                    observed_mortality[observed_mortality >= severe_threshold] -
                        predicted_mortality[observed_mortality >= severe_threshold]
                )^2)), NA_real_
            ),
            false_extreme_rate = mean(
                predicted_mortality[observed_mortality < tail_threshold] >=
                    tail_threshold
            ),
            .groups = 'drop'
        )
}

metrics_pooled <- cv_predictions |>
    group_by(scheme, tail_features, learner) |>
    metric_summary() |>
    mutate(programme_key = 'pooled', .before = n)
metrics_programme <- cv_predictions |>
    group_by(scheme, programme_key, tail_features, learner) |>
    metric_summary()
metrics <- bind_rows(metrics_pooled, metrics_programme)

event_metrics <- cv_predictions |>
    group_by(
        scheme, event_year, programme_key, tail_features, learner
    ) |>
    metric_summary()

tail_rows <- cv_predictions |>
    filter(observed_mortality >= severe_threshold) |>
    mutate(error = predicted_mortality - observed_mortality) |>
    arrange(error)

# A reef-event view prevents repeat observations at one reef (notably Snapper
# in 2024) from occupying several places in a diagnostic ranking. It is a
# diagnostic aggregation only; the fitting data and inclusion logic are
# unchanged.
reef_event_predictions <- cv_predictions |>
    group_by(
        scheme, tail_features, learner, programme_key,
        ReefID, ReefName, event_year, region_block
    ) |>
    summarise(
        observations = n(),
        observed_mortality = mean(observed_mortality),
        predicted_mortality = mean(predicted_mortality),
        ann_maxdhw = mean(ann_maxdhw),
        wqc_freqcc12 = mean(wqc_freqcc12),
        tc_min_distance_km = mean(tc_min_distance_km),
        .groups = 'drop'
    ) |>
    mutate(
        error = predicted_mortality - observed_mortality,
        underprediction = observed_mortality - predicted_mortality
    )

# Production influence is descriptive only; validation above determines value.
production_prepared <- prepare_tail_fold(data, data[0, , drop = FALSE])
production_models <- lapply(seq_along(feature_sets), function(i) {
    fit_tail_model(
        production_prepared$analysis, feature_sets[[i]],
        20261827L + i * 100L
    )
})
names(production_models) <- names(feature_sets)
production_influence <- bind_rows(lapply(names(production_models), function(name) {
    model <- production_models[[name]]
    bind_rows(lapply(c('gate', 'magnitude'), function(component_name) {
        component <- model[[component_name]]
        if (component$type != 'gbm') return(tibble())
        summary(component$model, plotit = FALSE) |>
            as_tibble() |>
            transmute(
                tail_features = name, component = component_name,
                predictor = var, relative_influence = rel.inf
            )
    }))
}))

stopifnot(
    all(is.finite(cv_predictions$predicted_mortality)),
    all(cv_predictions$predicted_mortality >= 0 &
            cv_predictions$predicted_mortality <= 1),
    nrow(cv_predictions |> filter(learner == 'INLA only')) ==
        nrow(data) * n_distinct(fold_map$scheme)
)

write_csv(cv_predictions, file.path(output_dir, 'cv_predictions.csv'), na = '')
write_csv(metrics, file.path(output_dir, 'metrics.csv'), na = '')
write_csv(event_metrics, file.path(output_dir, 'event_metrics.csv'), na = '')
write_csv(tail_rows, file.path(output_dir, 'severe_rows.csv'), na = '')
write_csv(
    reef_event_predictions,
    file.path(output_dir, 'reef_event_predictions.csv'), na = ''
)
write_csv(
    production_influence,
    file.path(output_dir, 'production_influence.csv'), na = ''
)
write_csv(runtime, file.path(output_dir, 'runtime.csv'), na = '')

print(metrics |> filter(programme_key == 'pooled'))
