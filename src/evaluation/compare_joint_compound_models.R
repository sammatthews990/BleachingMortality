# Compare the joint compound-stress candidates with the existing routed and
# reef-blocked benchmarks. The BRT core ablation isolates weather value from
# the effect of pooling observation programmes.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(posterior)
    library(readr)
    library(tidyr)
})
source('src/lib/joint_compound_model_helpers.R')

output_dir <- 'output/joint_compound_models'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

joint_brms <- read_csv(
    file.path(output_dir, 'joint_brms_predictions.csv'),
    show_col_types = FALSE
) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )
joint_brt <- read_csv(
    file.path(output_dir, 'joint_brt_predictions.csv'),
    show_col_types = FALSE
) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )
joint_core_brt <- read_csv(
    file.path(output_dir, 'joint_brt_core_predictions.csv'),
    show_col_types = FALSE
) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )

event_benchmark <- read_csv(
    'output/shape_aware_models/comparison_rows.csv', show_col_types = FALSE
) |>
    filter(learner == 'recommended_route') |>
    transmute(
        programme_key,
        source_observation_id = as.character(source_observation_id),
        ReefID, ReefName, event_year, region_block,
        joint_reef_fold = NA_integer_, ann_maxdhw, prop_acropora_pre,
        observed_pre_cover, log_coastal_rain30 = NA_real_,
        era5_wind_mean = NA_real_, era5_wind_calm_fraction = NA_real_,
        era5_coastal_distance_km = NA_real_, DISTURBANCE_TYPE = NA_character_,
        description = NA_character_, observed_mortality,
        predicted_occurrence, predicted_positive_loss, predicted_mortality,
        prediction_q05 = NA_real_, prediction_q50 = NA_real_,
        prediction_q95 = NA_real_, learner = 'routed_benchmark',
        scheme = 'leave_one_event_out', fold = as.character(fold),
        residual = observed_mortality - predicted_mortality,
        absolute_error = abs(residual)
    )

spatial_benchmark <- read_csv(
    'output/formal_models/model_comparison_predictions.csv',
    show_col_types = FALSE
) |>
    filter(
        scheme == 'reef_blocked_5fold', event_year == 2024,
        learner == 'formal_ensemble'
    ) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold),
        learner = 'formal_spatial_benchmark',
        scheme = 'reef_blocked_2024',
        residual = observed_mortality - predicted_mortality,
        absolute_error = abs(residual)
    ) |>
    left_join(
        joint_brt |>
            filter(scheme == 'reef_blocked_2024') |>
            select(
                programme_key, source_observation_id,
                joint_reef_fold,
                ann_maxdhw, prop_acropora_pre, observed_pre_cover,
                log_coastal_rain30, era5_wind_mean,
                era5_wind_calm_fraction, era5_coastal_distance_km,
                DISTURBANCE_TYPE, description
            ),
        by = c('programme_key', 'source_observation_id'),
        relationship = 'one-to-one'
    )

common_columns <- intersect(
    names(joint_brms),
    Reduce(intersect, list(
        names(joint_brt), names(joint_core_brt),
        names(event_benchmark), names(spatial_benchmark)
    ))
)
comparison <- bind_rows(
    joint_brms |> select(all_of(common_columns)),
    joint_brt |> select(all_of(common_columns)),
    joint_core_brt |> select(all_of(common_columns)),
    event_benchmark |> select(all_of(common_columns)),
    spatial_benchmark |> select(all_of(common_columns))
) |>
    mutate(
        severe = observed_mortality >= 0.2,
        nonsevere = !severe
    )

metrics <- comparison |>
    group_by(scheme, programme_key, learner) |>
    summarise(
        n = n(), severe_n = sum(severe),
        rmse = sqrt(mean(residual^2)),
        mae = mean(absolute_error),
        predictive_r2 = 1 - sum(residual^2) /
            sum((observed_mortality - mean(observed_mortality))^2),
        bias_observed_minus_predicted = mean(residual),
        severe_mae = mean(absolute_error[severe]),
        severe_bias = mean(residual[severe]),
        nonsevere_mae = mean(absolute_error[nonsevere]),
        .groups = 'drop'
    )

event_metrics <- comparison |>
    filter(scheme == 'leave_one_event_out') |>
    group_by(programme_key, learner, event_year) |>
    summarise(
        n = n(), severe_n = sum(severe),
        observed_mean = mean(observed_mortality),
        predicted_mean = mean(predicted_mortality),
        rmse = sqrt(mean(residual^2)),
        severe_mae = mean(absolute_error[severe]),
        severe_bias = mean(residual[severe]),
        .groups = 'drop'
    )

ablation <- metrics |>
    filter(learner %in% c('joint_core_brt', 'joint_compound_brt')) |>
    select(
        scheme, programme_key, learner,
        rmse, severe_mae, severe_bias, nonsevere_mae
    ) |>
    pivot_wider(
        names_from = learner,
        values_from = c(rmse, severe_mae, severe_bias, nonsevere_mae)
    ) |>
    mutate(
        delta_rmse_compound_minus_core =
            rmse_joint_compound_brt - rmse_joint_core_brt,
        delta_severe_mae_compound_minus_core =
            severe_mae_joint_compound_brt - severe_mae_joint_core_brt,
        delta_nonsevere_mae_compound_minus_core =
            nonsevere_mae_joint_compound_brt - nonsevere_mae_joint_core_brt
    )

production <- readRDS(
    'output/models/joint_compound/joint_brms_v2_production.rds'
)
fit <- production$fit
draws <- posterior::as_draws_df(fit)
effects <- tibble(
    effect = c(
        'Rainfall main effect on positive mortality',
        'DHW x rainfall on positive mortality',
        'Calm main effect on positive mortality',
        'DHW x calm on positive mortality',
        'Rainfall main effect on zero inflation',
        'DHW x rainfall on zero inflation',
        'Calm main effect on zero inflation',
        'DHW x calm on zero inflation'
    ),
    parameter = c(
        'b_log_coastal_rain30_z',
        'b_ann_maxdhw_z:log_coastal_rain30_z',
        'b_era5_wind_calm_fraction_z',
        'b_ann_maxdhw_z:era5_wind_calm_fraction_z',
        'b_zoi_log_coastal_rain30_z',
        'b_zoi_ann_maxdhw_z:log_coastal_rain30_z',
        'b_zoi_era5_wind_calm_fraction_z',
        'b_zoi_ann_maxdhw_z:era5_wind_calm_fraction_z'
    ),
    ecological_direction = c(
        'positive', 'positive', 'positive', 'positive',
        'negative', 'negative', 'negative', 'negative'
    )
) |>
    rowwise() |>
    mutate(
        estimate = mean(draws[[parameter]]),
        q025 = quantile(draws[[parameter]], 0.025),
        q50 = quantile(draws[[parameter]], 0.5),
        q975 = quantile(draws[[parameter]], 0.975),
        probability_ecological_direction = if_else(
            ecological_direction == 'positive',
            mean(draws[[parameter]] > 0),
            mean(draws[[parameter]] < 0)
        )
    ) |>
    ungroup()

manta_spatial <- comparison |>
    filter(
        scheme == 'reef_blocked_2024', programme_key == 'manta', severe
    ) |>
    select(
        source_observation_id, ReefID, ReefName, region_block,
        ann_maxdhw, prop_acropora_pre, log_coastal_rain30,
        era5_wind_calm_fraction, DISTURBANCE_TYPE, description,
        observed_mortality, learner, predicted_mortality, absolute_error
    ) |>
    pivot_wider(
        names_from = learner,
        values_from = c(predicted_mortality, absolute_error)
    ) |>
    mutate(
        coastal_rain30_mm = expm1(log_coastal_rain30),
        weather_gain_vs_core =
            absolute_error_joint_core_brt -
            absolute_error_joint_compound_brt,
        weather_gain_vs_formal =
            absolute_error_formal_spatial_benchmark -
            absolute_error_joint_compound_brt
    ) |>
    arrange(desc(weather_gain_vs_core))

manta_regions <- comparison |>
    filter(
        scheme == 'reef_blocked_2024', programme_key == 'manta', severe
    ) |>
    group_by(learner, region_block) |>
    summarise(
        n = n(), observed_mean = mean(observed_mortality),
        predicted_mean = mean(predicted_mortality),
        severe_mae = mean(absolute_error),
        severe_bias = mean(residual), .groups = 'drop'
    )

decision <- tribble(
    ~candidate, ~prediction_target, ~status, ~reason,
    'Joint compound BRMS', 'Future unseen event', 'Reject',
    'Worse RMSE and severe MAE for all programmes; positive DHW-rainfall synergy not supported.',
    'Joint compound BRT', 'Future unseen event', 'Reject',
    'Weather improves some average errors but worsens severe-tail error and does not beat the routed benchmark.',
    'Joint compound BRMS', '2024 spatial reconstruction', 'Reject',
    'Does not improve manta or MMP severe spatial error over the formal benchmark.',
    'Joint compound BRT', '2024 spatial reconstruction', 'Sensitivity only',
    'Weather ablation improves manta spatial RMSE and severe MAE, but the gain is framework- and programme-specific.'
)

write_csv(comparison, file.path(output_dir, 'comparison_rows.csv'))
write_csv(metrics, file.path(output_dir, 'comparison_metrics.csv'))
write_csv(event_metrics, file.path(output_dir, 'event_metrics.csv'))
write_csv(ablation, file.path(output_dir, 'weather_ablation.csv'))
write_csv(effects, file.path(output_dir, 'brms_weather_effects.csv'))
write_csv(manta_spatial, file.path(output_dir, 'manta_2024_severe_reefs.csv'))
write_csv(manta_regions, file.path(output_dir, 'manta_2024_regions.csv'))
write_csv(decision, file.path(output_dir, 'model_decision.csv'))

print(metrics)
print(ablation)
print(effects)
