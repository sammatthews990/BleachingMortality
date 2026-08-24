# Post-fit diagnostics for the decomposed annual transition candidates. Kept
# separate from model fitting so residual and modifier summaries can be rebuilt
# without refitting the BRTs.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

output_dir <- 'output/decomposed_annual_change'
predictions <- read_csv(
    file.path(output_dir, 'decomposed_change_predictions.csv'),
    show_col_types = FALSE
)

forward_outliers <- predictions |>
    filter(
        scheme == 'forward_2024_2025',
        learner == 'Mortality-augmented decomposed BRT'
    ) |>
    mutate(absolute_change_error_pp = abs(change_error_pp)) |>
    arrange(desc(absolute_change_error_pp))
write_csv(
    forward_outliers,
    file.path(output_dir, 'mortality_augmented_forward_outliers.csv')
)

event_summary <- predictions |>
    group_by(learner, scheme, event_year) |>
    summarise(
        transitions = n(),
        observed_change_pp = mean(cover_change_pp),
        severe_loss_n = sum(cover_change_pp <= -10),
        severe_loss_observed_pp = if_else(
            severe_loss_n > 0,
            mean(cover_change_pp[cover_change_pp <= -10]), NA_real_
        ),
        severe_loss_predicted_pp = if_else(
            severe_loss_n > 0,
            mean(predicted_change_pp[cover_change_pp <= -10]), NA_real_
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

candidate_predictions <- predictions |>
    select(
        scheme, fold, transition_id, programme_key, ReefID, ReefName,
        event_year, pre_cover, post_cover, cover_change_pp,
        sst_maxdhw, prop_acropora_pre,
        wqc_freqcc12, wqc_prior10_percentile,
        log_coastal_rain30, era5_wind_calm_fraction,
        cyc_maxHrs4mw, cot_idwmeanpertow,
        DISTURBANCE_TYPE, disturbance_text,
        learner, predicted_change_pp, expected_gain_pp, expected_loss_pp,
        change_error_pp
    ) |>
    pivot_wider(
        names_from = learner,
        values_from = c(
            predicted_change_pp, expected_gain_pp, expected_loss_pp,
            change_error_pp
        ),
        names_glue = '{.value}__{learner}'
    ) |>
    rename_with(~ gsub('[^A-Za-z0-9]+', '_', .x)) |>
    mutate(
        rich_minus_core_change_pp =
            predicted_change_pp_Modifier_rich_decomposed_BRT -
            predicted_change_pp_Core_decomposed_BRT,
        augmented_minus_rich_change_pp =
            predicted_change_pp_Mortality_augmented_decomposed_BRT -
            predicted_change_pp_Modifier_rich_decomposed_BRT,
        rich_improves_absolute_error =
            abs(change_error_pp_Modifier_rich_decomposed_BRT) <
            abs(change_error_pp_Core_decomposed_BRT),
        auxiliary_improves_absolute_error =
            abs(change_error_pp_Mortality_augmented_decomposed_BRT) <
            abs(change_error_pp_Modifier_rich_decomposed_BRT),
        severe_loss = cover_change_pp <= -10
    )
write_csv(
    candidate_predictions,
    file.path(output_dir, 'modifier_prediction_contrasts.csv')
)

contrast_summary <- candidate_predictions |>
    group_by(scheme, event_year, severe_loss) |>
    summarise(
        transitions = n(),
        mean_rich_minus_core_change_pp = mean(rich_minus_core_change_pp),
        rich_improvement_fraction = mean(rich_improves_absolute_error),
        mean_augmented_minus_rich_change_pp =
            mean(augmented_minus_rich_change_pp),
        auxiliary_improvement_fraction =
            mean(auxiliary_improves_absolute_error),
        .groups = 'drop'
    )
write_csv(
    contrast_summary,
    file.path(output_dir, 'modifier_contrast_summary.csv')
)

severe_2024 <- candidate_predictions |>
    filter(scheme == 'reef_blocked_5fold', event_year == 2024L, severe_loss) |>
    mutate(
        augmented_absolute_error_pp = abs(
            change_error_pp_Mortality_augmented_decomposed_BRT
        )
    ) |>
    arrange(desc(augmented_absolute_error_pp))
write_csv(
    severe_2024,
    file.path(output_dir, 'severe_2024_reef_diagnostics.csv')
)

numeric_diagnostics <- c(
    'sst_maxdhw', 'prop_acropora_pre', 'wqc_freqcc12',
    'wqc_prior10_percentile', 'log_coastal_rain30',
    'era5_wind_calm_fraction', 'cyc_maxHrs4mw', 'cot_idwmeanpertow'
)
residual_correlations <- bind_rows(lapply(numeric_diagnostics, function(x) {
    complete <- is.finite(severe_2024[[x]]) & is.finite(
        severe_2024$change_error_pp_Mortality_augmented_decomposed_BRT
    )
    tibble(
        predictor = x, n = sum(complete),
        spearman_with_residual = if (sum(complete) >= 5L) cor(
            severe_2024[[x]][complete],
            severe_2024$change_error_pp_Mortality_augmented_decomposed_BRT[complete],
            method = 'spearman'
        ) else NA_real_
    )
}))
write_csv(
    residual_correlations,
    file.path(output_dir, 'severe_2024_residual_correlations.csv')
)

# A fixed equal blend is reported as a sensitivity only. It is not tuned on
# these assessment outcomes and is not promoted without cross-fitting.
blend_diagnostics <- candidate_predictions |>
    mutate(
        predicted_change_pp = 0.5 * (
            predicted_change_pp_Modifier_rich_decomposed_BRT +
            predicted_change_pp_Mortality_augmented_decomposed_BRT
        ),
        error_pp = cover_change_pp - predicted_change_pp
    ) |>
    group_by(scheme) |>
    summarise(
        transitions = n(),
        rmse_pp = sqrt(mean(error_pp^2)),
        predictive_r2 = 1 - sum(error_pp^2) /
            sum((cover_change_pp - mean(cover_change_pp))^2),
        severe_mae_pp = mean(abs(error_pp)[cover_change_pp <= -10]),
        severe_bias_pp = mean(
            predicted_change_pp[cover_change_pp <= -10] -
                cover_change_pp[cover_change_pp <= -10]
        ),
        .groups = 'drop'
    )
write_csv(
    blend_diagnostics,
    file.path(output_dir, 'equal_blend_diagnostics.csv')
)

print(filter(
    contrast_summary,
    scheme == 'reef_blocked_5fold', event_year == 2024L
))
