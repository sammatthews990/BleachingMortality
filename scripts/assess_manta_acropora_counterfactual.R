# Quantify whether the five severe manta misses could be explained solely by
# inaccurate preceding Acropora proportions. Predictions use each reef's
# original reef-blocked 2024 core-BRT fold, not the production fit.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})
source('scripts/joint_compound_model_helpers.R')

output_dir <- 'output/rrn_target_reefs'
target_ids <- c('16-015', '14-126', '14-116a', '14-118', '15-030')
scenario_acropora <- c(0, 0.15, 0.30, 0.50, 0.70, 0.90, 1)

predict_component <- function(component, new_data, n_trees, type) {
    if (component$type == 'constant') {
        return(rep(component$value, nrow(new_data)))
    }
    predict(
        component$model, newdata = new_data,
        n.trees = n_trees, type = type
    )
}

predict_core_brt <- function(model, new_data, n_trees) {
    occurrence <- predict_component(
        model$occurrence, new_data, n_trees, 'response'
    )
    positive <- predict_component(
        model$magnitude, new_data, n_trees, 'response'
    )
    if (model$magnitude_mode == 'logit_gaussian') positive <- plogis(positive)
    pmin(pmax(occurrence, 0), 1) * pmin(pmax(positive, 0), 1)
}

fold_rows <- read_csv(
    'output/rrn_pressure_assessment/mortality_brt_paired_rows.csv',
    show_col_types = FALSE
) |>
    filter(
        scheme == 'reef_blocked_2024',
        programme_key == 'manta',
        event_year == 2024L,
        ReefID %in% target_ids
    ) |>
    distinct(ReefID, source_observation_id, .keep_all = TRUE) |>
    mutate(source_observation_id = as.character(source_observation_id)) |>
    select(
        ReefID, source_observation_id, fold,
        observed_mortality, predicted_mortality_core
    )

data <- load_joint_compound_rows() |>
    filter(
        programme_key == 'manta',
        event_year == 2024L,
        ReefID %in% target_ids
    ) |>
    inner_join(
        fold_rows,
        by = c('ReefID', 'source_observation_id'),
        relationship = 'one-to-one'
    )

predictions <- tibble()
for (row_index in seq_len(nrow(data))) {
    reef_row <- data[row_index, , drop = FALSE]
    fold <- as.character(reef_row$fold[[1]])
    reef_scenarios <- sort(unique(c(
        scenario_acropora, reef_row$prop_acropora_pre[[1]]
    )))
    cache_file <- file.path(
        'output/models/joint_compound',
        paste0('joint_brt_core_reef_blocked_2024_', fold, '.rds')
    )
    cached <- readRDS(cache_file)
    new_data <- reef_row[rep(1, length(reef_scenarios)), , drop = FALSE]
    new_data$prop_acropora_pre <- reef_scenarios
    n_trees <- as.integer(cached$tuning$best$n_trees[[1]])
    predictions <- bind_rows(
        predictions,
        tibble(
            ReefID = reef_row$ReefID[[1]],
            ReefName = reef_row$ReefName[[1]],
            fold = fold,
            observed_mortality = reef_row$observed_mortality[[1]],
            assigned_prop_acropora_pre = reef_row$prop_acropora_pre[[1]],
            acropora_source = reef_row$acropora_source[[1]],
            acropora_interpolated = reef_row$acropora_interpolated[[1]],
            scenario_prop_acropora_pre = reef_scenarios,
            predicted_mortality = predict_core_brt(
                cached$model, new_data, n_trees
            )
        )
    )
}

summary <- predictions |>
    group_by(ReefID, ReefName) |>
    summarise(
        observed_mortality = first(observed_mortality),
        assigned_prop_acropora_pre = first(assigned_prop_acropora_pre),
        acropora_source = first(acropora_source),
        acropora_interpolated = first(acropora_interpolated),
        prediction_at_assigned = predicted_mortality[
            abs(
                scenario_prop_acropora_pre -
                    first(assigned_prop_acropora_pre)
            ) < 1e-12
        ][1],
        prediction_at_50pct_acropora =
            predicted_mortality[scenario_prop_acropora_pre == 0.5],
        prediction_at_90pct_acropora =
            predicted_mortality[scenario_prop_acropora_pre == 0.9],
        maximum_prediction =
            max(predicted_mortality[scenario_prop_acropora_pre == 1]),
        gain_assigned_to_90pct =
            prediction_at_90pct_acropora - prediction_at_assigned,
        remaining_underprediction_at_90pct =
            first(observed_mortality) - prediction_at_90pct_acropora,
        .groups = 'drop'
    ) |>
    arrange(desc(observed_mortality))

write_csv(
    predictions,
    file.path(output_dir, 'manta_acropora_counterfactual_predictions.csv')
)
write_csv(
    summary,
    file.path(output_dir, 'manta_acropora_counterfactual_summary.csv')
)
print(summary, n = Inf, width = Inf)
