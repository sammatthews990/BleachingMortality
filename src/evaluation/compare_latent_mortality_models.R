# Compare the shared latent response against existing pooled BRMS and BRT
# predictions on the same outer validation observations.

suppressPackageStartupMessages({
    library(dplyr)
    library(posterior)
    library(readr)
    library(tidyr)
})
source('src/lib/latent_mortality_model_helpers.R')

output_dir <- 'output/latent_mortality'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

latent <- read_csv(
    file.path(output_dir, 'latent_predictions.csv'), show_col_types = FALSE
) |>
    transmute(
        programme_key, source_observation_id = as.character(source_observation_id),
        ReefID, ReefName, event_year, scheme, fold = as.character(fold),
        observed_mortality,
        predicted_mortality,
        prediction_q05, prediction_q95,
        learner = 'Shared latent BRMS'
    )

pooled <- read_csv(
    'output/joint_compound_models/joint_brms_predictions.csv',
    show_col_types = FALSE
) |>
    filter(scheme %in% unique(latent$scheme)) |>
    transmute(
        programme_key, source_observation_id = as.character(source_observation_id),
        ReefID, ReefName, event_year, scheme, fold = as.character(fold),
        observed_mortality,
        predicted_mortality,
        prediction_q05, prediction_q95,
        learner = 'Pooled-response BRMS'
    )

brt <- read_csv(
    'output/rrn_pressure_assessment/mortality_brt_paired_rows.csv',
    show_col_types = FALSE
) |>
    filter(scheme %in% unique(latent$scheme)) |>
    transmute(
        programme_key, source_observation_id = as.character(source_observation_id),
        ReefID, ReefName, event_year, scheme, fold = as.character(fold),
        observed_mortality,
        predicted_mortality = predicted_mortality_weather_rrn,
        prediction_q05 = NA_real_, prediction_q95 = NA_real_,
        learner = 'Weather + RRN BRT'
    )

predictions <- bind_rows(latent, pooled, brt)

metric_summary <- function(data, programme = FALSE) {
    groups <- c('learner', 'scheme')
    if (programme) groups <- c(groups, 'programme_key')
    data |>
        group_by(across(all_of(groups))) |>
        summarise(
            observations = n(),
            rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
            mae = mean(abs(observed_mortality - predicted_mortality)),
            severe_mae = mean(
                abs(observed_mortality - predicted_mortality)[
                    observed_mortality >= 0.2
                ]
            ),
            predictive_r2 = 1 -
                sum((observed_mortality - predicted_mortality)^2) /
                sum((observed_mortality - mean(observed_mortality))^2),
            bias_observed_minus_predicted = mean(
                observed_mortality - predicted_mortality
            ),
            interval90_coverage = if (all(is.na(prediction_q05))) {
                NA_real_
            } else {
                mean(
                    observed_mortality >= prediction_q05 &
                        observed_mortality <= prediction_q95
                )
            },
            .groups = 'drop'
        )
}

metrics <- bind_rows(
    metric_summary(predictions) |> mutate(programme_key = 'all'),
    metric_summary(predictions, programme = TRUE)
)
write_csv(predictions, file.path(output_dir, 'model_comparison_predictions.csv'))
write_csv(metrics, file.path(output_dir, 'model_comparison_metrics.csv'))

data <- load_joint_compound_rows()
overlap <- data |>
    distinct(ReefID, event_year, programme_key) |>
    mutate(value = 1L) |>
    pivot_wider(
        names_from = programme_key, values_from = value, values_fill = 0L
    ) |>
    count(ltmp, manta, mmp, name = 'reef_events')
write_csv(overlap, file.path(output_dir, 'programme_overlap.csv'))

paired <- data |>
    select(ReefID, event_year, programme_key, mortality_prop) |>
    group_by(ReefID, event_year, programme_key) |>
    summarise(mortality_prop = median(mortality_prop), .groups = 'drop') |>
    pivot_wider(names_from = programme_key, values_from = mortality_prop) |>
    filter(is.finite(ltmp), is.finite(manta))
paired_summary <- tibble(
    paired_reef_events = nrow(paired),
    ltmp_manta_correlation = cor(paired$ltmp, paired$manta),
    mean_ltmp_minus_manta = mean(paired$ltmp - paired$manta),
    median_absolute_difference = median(abs(paired$ltmp - paired$manta))
)
write_csv(paired_summary, file.path(output_dir, 'ltmp_manta_agreement.csv'))

production <- readRDS(
    'output/models/latent_mortality/latent_production.rds'
)
draws <- posterior::as_draws_df(production$fit)
effect_variables <- grep(
    '^b_(latenteta|bias)|^b_(phi|zoi|coi)_programme|^sd_',
    names(draws), value = TRUE
)
effects <- bind_rows(lapply(effect_variables, function(variable) {
    values <- as.numeric(draws[[variable]])
    tibble(
        variable,
        mean = mean(values), sd = sd(values),
        q025 = quantile(values, 0.025),
        q50 = quantile(values, 0.5),
        q975 = quantile(values, 0.975),
        probability_positive = mean(values > 0)
    )
}))
write_csv(effects, file.path(output_dir, 'latent_production_effects.csv'))

loo_result <- production$loo
r2 <- as.numeric(production$r2[, 'Estimate'])
fit_summary <- tibble(
    elpd_loo = loo_result$estimates['elpd_loo', 'Estimate'],
    elpd_loo_se = loo_result$estimates['elpd_loo', 'SE'],
    looic = loo_result$estimates['looic', 'Estimate'],
    looic_se = loo_result$estimates['looic', 'SE'],
    pareto_k_over_07 = sum(loo_result$diagnostics$pareto_k > 0.7),
    bayes_r2 = r2
)
write_csv(fit_summary, file.path(output_dir, 'latent_fit_summary.csv'))

# Smoke-test the programme-free operational prediction path. New identifiers
# ensure that no fitted reef-event effect can leak into the predictions.
event_example <- make_latent_event_rows(data) |>
    slice_head(n = 12) |>
    mutate(ReefID = paste0('new_reef_', row_number()))
operational <- predict_operational_latent_mortality(
    production$fit, event_example, production$preprocessing
)
write_csv(
    operational,
    file.path(output_dir, 'operational_prediction_smoke_test.csv')
)

print(filter(metrics, programme_key == 'all'))
print(paired_summary)
print(fit_summary)
