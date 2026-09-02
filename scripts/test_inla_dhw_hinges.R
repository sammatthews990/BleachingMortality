# Sensitivity test for the first DHW hinge in the selected operational INLA
# model. Compare 4+8 DHW with 6+8 DHW using identical data, folds, priors,
# covariates, interactions and latent structure.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(purrr)
    library(readr)
    library(tidyr)
})

Sys.setenv(INLA_ST_RUN = '0')
source('scripts/fit_inla_spatiotemporal_screen.R')
source('scripts/model_registry.R')
source('scripts/model_diagnostics.R')

root <- project_root()
hinge_dir <- file.path(root, 'output', 'inla_hinge_sensitivity')
dir.create(hinge_dir, recursive = TRUE, showWarnings = FALSE)

registered_best <- read_model_registry(root)$best_model
selected_id <- if (!is.null(registered_best$thermal_core_id)) {
    registered_best$thermal_core_id
} else registered_best$id
base_candidate <- candidates |>
    filter(candidate == selected_id)
stopifnot(nrow(base_candidate) == 1L)

hinge_candidates <- bind_rows(
    base_candidate |>
        mutate(candidate = 'selected_hinge4_plus8', first_hinge_dhw = 4),
    base_candidate |>
        mutate(candidate = 'test_hinge6_plus8', first_hinge_dhw = 6)
)

# Keep the established column name so all interactions and model-building
# machinery are unchanged; only its ecological definition changes.
base_apply_candidate_dhw <- apply_candidate_dhw
apply_candidate_dhw <- function(rows, candidate) {
    first_hinge <- if ('first_hinge_dhw' %in% names(candidate)) {
        as.numeric(candidate$first_hinge_dhw[[1]])
    } else 4
    base_apply_candidate_dhw(rows, candidate) |>
        mutate(dhw_excess4 = pmax(ann_maxdhw - first_hinge, 0))
}

full_criteria <- tibble()
full_fixed <- tibble()
predictions <- tibble()

for (candidate_index in seq_len(nrow(hinge_candidates))) {
    candidate <- hinge_candidates[candidate_index, ]
    first_hinge <- candidate$first_hinge_dhw[[1]]
    full_path <- file.path(
        hinge_dir, paste0('full_hinge', first_hinge, '_plus8.rds')
    )
    if (file.exists(full_path)) {
        full_result <- readRDS(full_path)
    } else {
        full_result <- fit_candidate(
            data, data[1, , drop = FALSE], candidate,
            compute_criteria = TRUE,
            seed = 20260831L + candidate_index
        )
        saveRDS(full_result, full_path)
    }
    full_summary <- extract_fit_summary(full_result, candidate)
    full_criteria <- bind_rows(
        full_criteria,
        full_summary$criteria |> mutate(first_hinge_dhw = first_hinge)
    )
    full_fixed <- bind_rows(
        full_fixed,
        full_summary$fixed |> mutate(first_hinge_dhw = first_hinge)
    )

    for (scheme in c('leave_one_event_out', 'reef_blocked_5fold')) {
        folds <- if (scheme == 'leave_one_event_out') {
            event_years
        } else {
            sort(unique(data$joint_reef_fold))
        }
        for (fold in folds) {
            if (scheme == 'leave_one_event_out') {
                assessment <- data |> filter(event_year == fold)
                analysis <- data |> filter(event_year != fold)
            } else {
                held_reefs <- unique(data$ReefID[data$joint_reef_fold == fold])
                assessment <- data |> filter(ReefID %in% held_reefs)
                analysis <- data |> filter(!ReefID %in% held_reefs)
            }
            cache_path <- file.path(
                hinge_dir,
                paste0(
                    'cv_hinge', first_hinge, '_plus8_', scheme, '_',
                    fold, '.rds'
                )
            )
            if (file.exists(cache_path)) {
                fold_result <- readRDS(cache_path)
            } else {
                fitted <- fit_candidate(
                    analysis, assessment, candidate,
                    compute_criteria = FALSE,
                    seed = 20260831L + candidate_index * 100L + as.integer(fold)
                )
                fold_result <- list(
                    predictions = fitted$predictions,
                    elapsed_seconds = fitted$elapsed_seconds
                )
                saveRDS(fold_result, cache_path)
            }
            predictions <- bind_rows(
                predictions,
                fold_result$predictions |>
                    mutate(
                        candidate = candidate$candidate,
                        first_hinge_dhw = first_hinge,
                        scheme = scheme,
                        fold = as.character(fold)
                    )
            )
        }
    }
}

metrics <- bind_rows(
    predictions |>
        group_by(candidate, first_hinge_dhw, scheme) |>
        metric_summary() |>
        mutate(programme_key = 'all'),
    predictions |>
        group_by(candidate, first_hinge_dhw, scheme, programme_key) |>
        metric_summary()
)
comparison <- metrics |>
    filter(programme_key == 'all') |>
    left_join(full_criteria, by = c('candidate', 'first_hinge_dhw')) |>
    arrange(scheme, rmse)

write_csv(predictions, file.path(hinge_dir, 'cv_predictions.csv'))
write_csv(metrics, file.path(hinge_dir, 'validation_metrics.csv'))
write_csv(comparison, file.path(hinge_dir, 'hinge_comparison.csv'))
write_csv(full_fixed, file.path(hinge_dir, 'full_fixed_effects.csv'))
write_csv(full_criteria, file.path(hinge_dir, 'full_criteria.csv'))

metric_plot_data <- comparison |>
    select(candidate, first_hinge_dhw, scheme, rmse, severe_rmse,
           predictive_r2, false_extreme_rate) |>
    pivot_longer(
        c(rmse, severe_rmse, predictive_r2, false_extreme_rate),
        names_to = 'metric', values_to = 'value'
    ) |>
    mutate(
        hinge = paste0(first_hinge_dhw, '+8 DHW'),
        metric = factor(
            metric,
            levels = c('rmse', 'severe_rmse', 'predictive_r2',
                       'false_extreme_rate'),
            labels = c('RMSE', 'Severe-event RMSE', 'Predictive R-squared',
                       'False-extreme rate')
        )
    )
metric_plot <- ggplot(
    metric_plot_data,
    aes(hinge, value, fill = hinge)
) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(aes(label = sprintf('%.3f', value)), vjust = -0.35, size = 3.2) +
    facet_grid(metric ~ scheme, scales = 'free_y') +
    scale_fill_manual(values = c('4+8 DHW' = '#0072B2', '6+8 DHW' = '#D55E00')) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
        title = 'First DHW hinge sensitivity: held-out performance',
        subtitle = 'All covariates, interactions, priors and latent effects are held constant',
        x = 'Thermal spline hinges', y = NULL
    ) +
    theme_bw(base_size = 11)
save_figure_bundle(
    metric_plot, 'Fig-INLA-05_hinge_validation', metric_plot_data,
    'Held-out model performance comparing thermal spline hinges at 4+8 and 6+8 DHW.',
    'Tests whether moving the first change in thermal slope from 4 to 6 DHW improves transfer to unseen events and reefs.',
    'The alternative also moves the onset of DHW interactions with Acropora, novelty, freshwater and cloud; it is therefore a whole-response sensitivity test.',
    'dhw_hinge_sensitivity', 'INLA', 'model_comparison',
    'operational_candidate_test', root, TRUE,
    code_source = 'scripts/test_inla_dhw_hinges.R'
)

response_data <- predictions |>
    filter(scheme == 'leave_one_event_out') |>
    mutate(
        dhw_bin = cut(
            ann_maxdhw, breaks = seq(0, 22, by = 1),
            include.lowest = TRUE, right = FALSE
        )
    ) |>
    group_by(candidate, first_hinge_dhw, dhw_bin) |>
    summarise(
        dhw = mean(ann_maxdhw),
        observed = mean(observed_mortality),
        predicted = mean(predicted_mortality),
        n = n(), .groups = 'drop'
    ) |>
    filter(n >= 3) |>
    pivot_longer(c(observed, predicted), names_to = 'series', values_to = 'mortality') |>
    mutate(hinge = paste0(first_hinge_dhw, '+8 DHW'))
response_plot <- ggplot(
    response_data,
    aes(dhw, mortality, colour = series, linetype = series)
) +
    geom_vline(
        data = distinct(response_data, hinge, first_hinge_dhw),
        aes(xintercept = first_hinge_dhw), colour = 'grey55', linetype = 3,
        inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 8, colour = 'grey55', linetype = 3) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(size = n), alpha = 0.65) +
    facet_wrap(~ hinge, ncol = 2) +
    scale_colour_manual(values = c(observed = 'black', predicted = '#D55E00')) +
    scale_linetype_manual(values = c(observed = 2, predicted = 1)) +
    coord_cartesian(ylim = c(0, 0.8)) +
    labs(
        title = 'Held-out mortality response under alternative DHW hinges',
        subtitle = 'One-DHW bins; points are shown where at least three observations are available',
        x = 'Local-first annual maximum DHW', y = 'Mean relative mortality',
        colour = NULL, linetype = NULL, size = 'Observations'
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = 'bottom')
save_figure_bundle(
    response_plot, 'Fig-INLA-06_hinge_response', response_data,
    'Binned held-out observed and predicted mortality for the 4+8 and 6+8 DHW spline candidates.',
    'Shows where moving the first hinge changes practical predictions across the observed DHW range.',
    'This is a binned validation diagnostic rather than a covariate-adjusted posterior dose-response curve.',
    'dhw_hinge_sensitivity', 'INLA', 'hinge_response',
    'operational_candidate_test', root, TRUE,
    code_source = 'scripts/test_inla_dhw_hinges.R'
)

write_figure_readme(root)
message('Wrote hinge sensitivity outputs to ', hinge_dir)
