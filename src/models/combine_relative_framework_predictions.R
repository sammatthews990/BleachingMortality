# Combine untouched reef-blocked predictions from the INLA beta--Bernoulli
# screen and the joint-modifier BRT screen. Equal weights are deliberately
# fixed in advance; optimising them here would contaminate validation.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

output_dir <- 'output/joint_extreme_relative'

brt <- read_csv(
    file.path(output_dir, 'ensemble_predictions.csv'),
    show_col_types = FALSE
) |>
    filter(learner %in% c(
        'Standard two-part BRT',
        'Conservative tail ensemble',
        'Soft-gated tail ensemble'
    )) |>
    select(
        programme_key, source_observation_id, ReefID, ReefName, event_year,
        ann_maxdhw, observed_mortality, learner, ensemble_prediction
    )

# Selected in the existing reef-blocked INLA history screen.
inla_choice <- tribble(
    ~programme_key, ~candidate,
    'ltmp', 'legacy_history',
    'manta', 'paper_history',
    'mmp', 'legacy_history'
)

inla <- read_csv(
    'output/inla_history_screen/predictions.csv',
    show_col_types = FALSE
) |>
    filter(scheme == 'reef_blocked_5fold') |>
    inner_join(inla_choice, by = c('programme_key', 'candidate')) |>
    select(
        programme_key, source_observation_id,
        inla_candidate = candidate,
        inla_prediction = predicted_mortality
    )

stopifnot(
    nrow(inla) == n_distinct(
        paste(inla$programme_key, inla$source_observation_id)
    )
)

joined <- brt |>
    left_join(inla, by = c('programme_key', 'source_observation_id'))

if (anyNA(joined$inla_prediction)) {
    stop('INLA and BRT held-out observations do not match')
}

identity_columns <- c(
    'programme_key', 'source_observation_id', 'ReefID', 'ReefName',
    'event_year', 'ann_maxdhw', 'observed_mortality'
)

framework_predictions <- bind_rows(
    joined |>
        transmute(
            across(all_of(identity_columns)), learner,
            prediction = ensemble_prediction
        ),
    joined |>
        distinct(
            across(all_of(identity_columns)),
            inla_candidate, inla_prediction
        ) |>
        transmute(
            across(all_of(identity_columns)),
            learner = 'INLA beta + Bernoulli',
            prediction = inla_prediction
        ),
    joined |>
        mutate(
            learner = paste0('Equal INLA + ', learner),
            prediction = 0.5 * inla_prediction +
                0.5 * ensemble_prediction
        ) |>
        select(all_of(identity_columns), learner, prediction)
)

metric_summary <- function(rows) {
    rows |>
        summarise(
            n = n(),
            severe_n = sum(observed_mortality >= 0.50),
            rmse = sqrt(mean((observed_mortality - prediction)^2)),
            mae = mean(abs(observed_mortality - prediction)),
            predictive_r2 = 1 -
                sum((observed_mortality - prediction)^2) /
                sum((observed_mortality - mean(observed_mortality))^2),
            bias = mean(prediction - observed_mortality),
            severe_observed = if_else(
                severe_n > 0,
                mean(observed_mortality[observed_mortality >= 0.50]),
                NA_real_
            ),
            severe_predicted = if_else(
                severe_n > 0,
                mean(prediction[observed_mortality >= 0.50]),
                NA_real_
            ),
            severe_rmse = if_else(
                severe_n > 0,
                sqrt(mean((observed_mortality[observed_mortality >= 0.50] -
                               prediction[observed_mortality >= 0.50])^2)),
                NA_real_
            ),
            false_extreme_rate = mean(
                prediction[observed_mortality < 0.30] >= 0.30
            ),
            .groups = 'drop'
        )
}

framework_metrics <- framework_predictions |>
    group_by(programme_key, learner) |>
    metric_summary()

framework_metrics_2024 <- framework_predictions |>
    filter(event_year == 2024) |>
    group_by(programme_key, learner) |>
    metric_summary()

pooled_metrics <- framework_predictions |>
    group_by(learner) |>
    metric_summary()

pooled_metrics_2024 <- framework_predictions |>
    filter(event_year == 2024) |>
    group_by(learner) |>
    metric_summary()

comparison_wide <- framework_predictions |>
    filter(learner %in% c(
        'Standard two-part BRT',
        'Conservative tail ensemble'
    )) |>
    select(
        all_of(identity_columns), learner, prediction
    ) |>
    pivot_wider(names_from = learner, values_from = prediction) |>
    mutate(
        standard_absolute_error = abs(
            observed_mortality - `Standard two-part BRT`
        ),
        conservative_absolute_error = abs(
            observed_mortality - `Conservative tail ensemble`
        ),
        error_reduction = standard_absolute_error -
            conservative_absolute_error
    ) |>
    arrange(desc(error_reduction))

stopifnot(
    all(is.finite(framework_predictions$prediction)),
    all(framework_predictions$prediction >= 0 &
            framework_predictions$prediction <= 1)
)

write_csv(
    framework_predictions,
    file.path(output_dir, 'framework_predictions.csv'), na = ''
)
write_csv(
    framework_metrics,
    file.path(output_dir, 'framework_metrics.csv'), na = ''
)
write_csv(
    framework_metrics_2024,
    file.path(output_dir, 'framework_metrics_2024.csv'), na = ''
)
write_csv(
    pooled_metrics,
    file.path(output_dir, 'framework_metrics_pooled.csv'), na = ''
)
write_csv(
    pooled_metrics_2024,
    file.path(output_dir, 'framework_metrics_pooled_2024.csv'), na = ''
)
write_csv(
    comparison_wide,
    file.path(output_dir, 'tail_error_changes.csv'), na = ''
)

print(framework_metrics)
