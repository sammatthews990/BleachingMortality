# Full-data explanatory comparison: event-specific compound INLA response
# versus a compound model containing only local-first DHW spline terms and
# programme observation layers. These are apparent/full-data predictive
# metrics, not future-event validation.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tibble)
})
source('scripts/fit_inla_explanatory_event_dhw.R')

comparison_dir <- file.path(root, 'output', 'explanatory_event_dhw_comparison')
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

event_candidate <- candidate |>
    mutate(candidate = 'explanatory_event_dhw_full_data')
dhw_candidate <- tibble(
    candidate = 'dhw_only_full_data',
    spatial_structure = 'none', temporal_dhw_slope = FALSE,
    feature_set = 'core', partial_pool = FALSE,
    lizard_dhw_adjustment = 'local_first'
)

event_make_formula <- make_formula
make_formula <- function(candidate_row) {
    if (identical(as.character(candidate_row$candidate), 'dhw_only_full_data')) {
        return(as.formula(paste(
            'response ~ -1 +',
            paste(
                c(
                    paste0('layer_', seq_along(families)),
                    'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z'
                ),
                collapse = ' + '
            )
        )))
    }
    event_make_formula(candidate_row)
}

fit_or_load <- function(candidate_row, filename, seed) {
    path <- file.path(comparison_dir, filename)
    if (file.exists(path)) return(readRDS(path))
    fitted <- fit_candidate(
        data, data, candidate_row,
        compute_criteria = TRUE, seed = seed
    )
    saveRDS(fitted, path)
    fitted
}
event_result <- fit_or_load(
    event_candidate, 'event_dhw_full_data_predictions.rds', 20260911L
)
dhw_result <- fit_or_load(
    dhw_candidate, 'dhw_only_full_data_predictions.rds', 20260912L
)

metric_row <- function(result, model_id) {
    predictions <- result$predictions
    tibble(
        model = model_id,
        n = nrow(predictions),
        rmse = sqrt(mean(
            (predictions$observed_mortality -
                 predictions$predicted_mortality)^2
        )),
        mae = mean(abs(
            predictions$observed_mortality -
                predictions$predicted_mortality
        )),
        predictive_r2 = 1 - sum(
            (predictions$observed_mortality -
                 predictions$predicted_mortality)^2
        ) / sum(
            (predictions$observed_mortality -
                 mean(predictions$observed_mortality))^2
        ),
        correlation_squared = cor(
            predictions$observed_mortality,
            predictions$predicted_mortality
        )^2,
        waic = result$fit$waic$waic,
        dic = result$fit$dic$dic
    )
}
metrics <- bind_rows(
    metric_row(event_result, 'Full event-DHW ecological model'),
    metric_row(dhw_result, 'DHW-only compound model')
)
write_csv(metrics, file.path(comparison_dir, 'model_comparison.csv'))

plot_data <- bind_rows(
    event_result$predictions |>
        transmute(
            model = 'Full event-DHW ecological model',
            observed = observed_mortality,
            predicted = predicted_mortality,
            programme = programme_key, event_year
        ),
    dhw_result$predictions |>
        transmute(
            model = 'DHW-only compound model',
            observed = observed_mortality,
            predicted = predicted_mortality,
            programme = programme_key, event_year
        )
) |>
    left_join(
        metrics |> select(model, predictive_r2),
        by = 'model'
    ) |>
    mutate(
        panel = paste0(
            model, '\nFull-data predictive R-squared = ',
            sprintf('%.3f', predictive_r2)
        )
    )
write_csv(plot_data, file.path(comparison_dir, 'predictions.csv'))
plot <- ggplot(plot_data, aes(predicted, observed, colour = programme)) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey40') +
    geom_point(alpha = 0.55, size = 1.7) +
    facet_wrap(~ panel, ncol = 2) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
        title = 'Full-data explanatory fit versus DHW-only reference',
        subtitle = 'Known events and reefs; this is not leave-one-event-out performance',
        x = 'Full-data predicted relative mortality',
        y = 'Observed relative mortality', colour = 'Programme'
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = 'bottom')
save_figure_bundle(
    plot, 'Fig-INLA-09_explanatory_vs_dhw_only', plot_data,
    'Observed versus full-data fitted mortality for the explanatory event-DHW ecological model and a programme-layer compound model using only the local-first DHW spline.',
    'Quantifies the additional apparent variation explained by event-specific ecological, disturbance and spatial structure relative to thermal exposure alone.',
    'These predictive R-squared values reuse the fitting data and known event/reef effects. They are explanatory apparent fit, not future-event forecast skill.',
    'explanatory_event_dhw_comparison', 'INLA', 'model_comparison',
    'full_data_explanatory', root, TRUE,
    code_source = 'scripts/compare_explanatory_event_dhw_r2.R'
)
write_figure_readme(root)
message('Wrote explanatory versus DHW-only comparison')
