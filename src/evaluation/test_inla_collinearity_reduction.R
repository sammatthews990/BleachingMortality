# Controlled validation of a collinearity-reduced version of the selected
# operational INLA model. Structural spline bases/interactions are retained;
# redundant ecological summaries are removed as complete mechanism pairs.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(tidyr)
})

Sys.setenv(INLA_ST_RUN = '0')
source('src/models/fit_inla_spatiotemporal_screen.R')
source('src/lib/model_registry.R')
source('src/lib/model_diagnostics.R')

root <- project_root()
audit_dir <- file.path(root, 'output', 'inla_collinearity_sensitivity')
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
registered_best <- read_model_registry(root)$best_model
selected_id <- if (!is.null(registered_best$thermal_core_id)) {
    registered_best$thermal_core_id
} else registered_best$id
candidate <- candidates |>
    filter(candidate == selected_id) |>
    mutate(
        candidate = 'test_collinearity_reduced',
        collinearity_variant = 'reduced'
    )
stopifnot(nrow(candidate) == 1L)

removed_terms <- c(
    'dhw10_load4_z',
    'wqc_10yr_sum_z',
    'dhw_x_wqc_cumulative'
)

decisions <- tribble(
    ~term_a, ~term_b, ~correlation, ~decision, ~reason,
    'ann_maxdhw_z', 'dhw_excess4_z', 0.966, 'retain_both_structural',
    'Required spline basis; assess the thermal mechanism jointly.',
    'dhw_excess4_z', 'dhw_excess8_z', 0.887, 'retain_both_structural',
    'Required spline basis; assess the thermal mechanism jointly.',
    'dhw_events_since2016_n6_z', 'dhw10_load4_z', 0.935,
    'remove_dhw10_load4_z',
    'Event count matches the pre-specified repeat-exposure mechanism and literature formulation.',
    'dhw_events_since2016_n6_z', 'dhw_years_since_last_n6_capped8_z', -0.815,
    'retain_frequency_and_recency',
    'Frequency and recovery interval are distinct pre-specified mechanisms; require ablation review if instability persists.',
    'wqc_freqcc12_z', 'wqc_10yr_sum_z', 0.950,
    'remove_inclusive_wqc_10yr_sum',
    'The rolling sum includes current WQC, mechanically duplicating the current exposure signal.',
    'dhw_x_wqc_current', 'dhw_x_wqc_cumulative', 0.973,
    'remove_cumulative_interaction',
    'Remove the interaction paired with the redundant inclusive cumulative main effect.'
)
write_csv(decisions, file.path(audit_dir, 'collinearity_decisions.csv'))

base_candidate_shared_terms <- candidate_shared_terms
candidate_shared_terms <- function(candidate_row) {
    terms <- base_candidate_shared_terms(candidate_row)
    if ('collinearity_variant' %in% names(candidate_row) &&
            identical(as.character(candidate_row$collinearity_variant), 'reduced')) {
        setdiff(terms, removed_terms)
    } else terms
}

# The expanded programme-specific freshwater deviation must follow the same
# reduction and therefore excludes the cumulative WQC interaction.
active_reduced <- FALSE
base_make_fixed_rows <- make_fixed_rows
make_fixed_rows <- function(rows, component, observed_event_keys,
                            observed_event_indices) {
    if (isTRUE(active_reduced)) {
        rows <- rows |>
            mutate(
                freshwater_amplification_weight = (
                    dhw_x_rainfall + dhw_x_wqc_current +
                        dhw_x_wqc_relative
                ) / sqrt(3)
            )
    }
    base_make_fixed_rows(
        rows, component, observed_event_keys, observed_event_indices
    )
}
base_fit_candidate <- fit_candidate
fit_candidate <- function(analysis, assessment, candidate,
                          compute_criteria = FALSE, seed = 1L) {
    old_active <- active_reduced
    active_reduced <<- 'collinearity_variant' %in% names(candidate) &&
        identical(as.character(candidate$collinearity_variant), 'reduced')
    on.exit(active_reduced <<- old_active, add = TRUE)
    base_fit_candidate(
        analysis, assessment, candidate,
        compute_criteria = compute_criteria, seed = seed
    )
}

full_path <- file.path(audit_dir, 'full_collinearity_reduced.rds')
if (file.exists(full_path)) {
    full_result <- readRDS(full_path)
} else {
    full_result <- fit_candidate(
        data, data[1, , drop = FALSE], candidate,
        compute_criteria = TRUE, seed = 20260901L
    )
    saveRDS(full_result, full_path)
}
full_summary <- extract_fit_summary(full_result, candidate)
write_csv(full_summary$criteria, file.path(audit_dir, 'full_criteria.csv'))
write_csv(full_summary$fixed, file.path(audit_dir, 'full_fixed_effects.csv'))

reduced_predictions <- tibble()
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
            audit_dir, paste0('cv_reduced_', scheme, '_', fold, '.rds')
        )
        if (file.exists(cache_path)) {
            fold_result <- readRDS(cache_path)
        } else {
            fitted <- fit_candidate(
                analysis, assessment, candidate,
                compute_criteria = FALSE,
                seed = 20260901L + as.integer(fold)
            )
            fold_result <- list(predictions = fitted$predictions)
            saveRDS(fold_result, cache_path)
        }
        reduced_predictions <- bind_rows(
            reduced_predictions,
            fold_result$predictions |>
                mutate(
                    candidate = 'collinearity_reduced',
                    scheme = scheme, fold = as.character(fold)
                )
        )
    }
}

baseline_predictions <- read_csv(
    file.path(root, 'output', 'inla_hinge_sensitivity', 'cv_predictions.csv'),
    show_col_types = FALSE
) |>
    filter(candidate == 'selected_hinge4_plus8') |>
    mutate(
        candidate = 'selected_full',
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )
reduced_predictions <- reduced_predictions |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )
predictions <- bind_rows(baseline_predictions, reduced_predictions)
metrics <- predictions |>
    group_by(candidate, scheme) |>
    metric_summary()
baseline_criteria <- read_csv(
    file.path(root, 'output', 'inla_hinge_sensitivity', 'full_criteria.csv'),
    show_col_types = FALSE
) |>
    filter(candidate == 'selected_hinge4_plus8') |>
    transmute(candidate = 'selected_full', waic, dic, log_marginal_likelihood,
              mean_neg_log_cpo)
criteria <- bind_rows(
    baseline_criteria,
    full_summary$criteria |>
        transmute(
            candidate = 'collinearity_reduced', waic, dic,
            log_marginal_likelihood, mean_neg_log_cpo
        )
)
comparison <- metrics |>
    left_join(criteria, by = 'candidate') |>
    arrange(scheme, rmse)
write_csv(predictions, file.path(audit_dir, 'cv_predictions.csv'))
write_csv(comparison, file.path(audit_dir, 'model_comparison.csv'))

plot_data <- comparison |>
    select(candidate, scheme, rmse, severe_rmse, predictive_r2,
           false_extreme_rate) |>
    pivot_longer(
        c(rmse, severe_rmse, predictive_r2, false_extreme_rate),
        names_to = 'metric', values_to = 'value'
    ) |>
    mutate(
        model = recode(
            candidate,
            selected_full = 'Selected full',
            collinearity_reduced = 'Reduced'
        ),
        metric = factor(
            metric,
            levels = c('rmse', 'severe_rmse', 'predictive_r2',
                       'false_extreme_rate'),
            labels = c('RMSE', 'Severe-event RMSE', 'Predictive R-squared',
                       'False-extreme rate')
        )
    )
plot <- ggplot(plot_data, aes(model, value, fill = model)) +
    geom_col(width = 0.68, show.legend = FALSE) +
    geom_text(aes(label = sprintf('%.3f', value)), vjust = -0.35, size = 3.1) +
    facet_grid(metric ~ scheme, scales = 'free_y') +
    scale_fill_manual(values = c('Selected full' = '#0072B2', 'Reduced' = '#009E73')) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
        title = 'Collinearity-reduced candidate: held-out performance',
        subtitle = 'Removed ten-year DHW load and inclusive cumulative-WQC main/interaction terms',
        x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11)
save_figure_bundle(
    plot, 'Fig-INLA-08_collinearity_reduction', plot_data,
    'Held-out performance of the selected full model and a candidate with redundant heat-history and inclusive cumulative-WQC terms removed.',
    'Tests whether reducing genuine ecological redundancy preserves or improves transfer to unseen events and reefs.',
    'Structural DHW spline bases and their modifier interactions remain correlated by construction and are assessed jointly rather than pruned.',
    'collinearity_reduced_candidate', 'INLA', 'model_comparison',
    'operational_candidate_test', root, TRUE,
    code_source = 'src/evaluation/test_inla_collinearity_reduction.R'
)
write_figure_readme(root)
message('Wrote collinearity-reduction sensitivity outputs')
