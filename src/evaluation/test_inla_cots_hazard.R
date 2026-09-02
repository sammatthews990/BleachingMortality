# Controlled test of a complementary COTS hazard representation.
#
# The hindcast estimates Pr(COTS density > 0.22 per tow), whereas the RRN
# interval maximum describes realised pressure.  The candidate retains the
# probability term for outbreak occurrence and replaces log intensity with a
# probability-weighted excess above the hindcast outbreak threshold
# (0.22 COTS/tow).

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
audit_dir <- file.path(root, 'output', 'inla_cots_hazard_sensitivity')
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

cots_outbreak_threshold <- 0.22
data <- data |>
    mutate(
        cots_severe_excess = log1p(pmax(
            cot_interval_idw_max - cots_outbreak_threshold, 0
        )),
        cots_expected_severe_burden =
            cots_outbreak_probability * cots_severe_excess
    )

raw_predictors <- unique(c(raw_predictors, 'cots_expected_severe_burden'))
shared_terms <- unique(c(shared_terms, 'cots_expected_severe_burden_z'))

registered_best <- read_model_registry(root)$best_model
selected_id <- if (!is.null(registered_best$thermal_core_id)) {
    registered_best$thermal_core_id
} else registered_best$id
candidate <- candidates |>
    filter(candidate == selected_id) |>
    mutate(candidate = 'test_cots_probability_severity')
stopifnot(nrow(candidate) == 1L)

base_candidate_shared_terms <- candidate_shared_terms
candidate_shared_terms <- function(candidate_row) {
    terms <- base_candidate_shared_terms(candidate_row)
    if (identical(
        as.character(candidate_row$candidate),
        'test_cots_probability_severity'
    )) {
        c(
            setdiff(terms, 'log1p_cot_interval_idw_max_z'),
            'cots_expected_severe_burden_z'
        )
    } else terms
}

# The programme-specific COTS deviation uses the same occurrence + severity
# representation as the shared ecological effect.
active_combined_cots <- FALSE
base_make_fixed_rows <- make_fixed_rows
make_fixed_rows <- function(rows, component, observed_event_keys,
                            observed_event_indices) {
    if (isTRUE(active_combined_cots)) {
        rows <- rows |>
            mutate(
                cots_hazard_weight = (
                    cots_outbreak_probability_z +
                        cots_expected_severe_burden_z
                ) / sqrt(2)
            )
    }
    base_make_fixed_rows(
        rows, component, observed_event_keys, observed_event_indices
    )
}

base_fit_candidate <- fit_candidate
fit_candidate <- function(analysis, assessment, candidate,
                          compute_criteria = FALSE, seed = 1L) {
    old_active <- active_combined_cots
    active_combined_cots <<- identical(
        as.character(candidate$candidate),
        'test_cots_probability_severity'
    )
    on.exit(active_combined_cots <<- old_active, add = TRUE)
    base_fit_candidate(
        analysis, assessment, candidate,
        compute_criteria = compute_criteria, seed = seed
    )
}

full_path <- file.path(audit_dir, 'full_cots_probability_severity.rds')
if (file.exists(full_path)) {
    full_result <- readRDS(full_path)
} else {
    full_result <- fit_candidate(
        data, data[1, , drop = FALSE], candidate,
        compute_criteria = TRUE, seed = 20260902L
    )
    saveRDS(full_result, full_path)
}
full_summary <- extract_fit_summary(full_result, candidate)
write_csv(full_summary$criteria, file.path(audit_dir, 'full_criteria.csv'))
write_csv(full_summary$fixed, file.path(audit_dir, 'full_fixed_effects.csv'))

candidate_predictions <- tibble()
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
            audit_dir, paste0('cv_cots_', scheme, '_', fold, '.rds')
        )
        if (file.exists(cache_path)) {
            fold_result <- readRDS(cache_path)
        } else {
            fitted <- fit_candidate(
                analysis, assessment, candidate,
                compute_criteria = FALSE,
                seed = 20260902L + as.integer(fold)
            )
            fold_result <- list(predictions = fitted$predictions)
            saveRDS(fold_result, cache_path)
        }
        candidate_predictions <- bind_rows(
            candidate_predictions,
            fold_result$predictions |>
                mutate(
                    candidate = 'cots_probability_severity',
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
        candidate = 'selected_additive_cots',
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )
candidate_predictions <- candidate_predictions |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )
predictions <- bind_rows(baseline_predictions, candidate_predictions)

metrics <- predictions |>
    group_by(candidate, scheme) |>
    metric_summary()
baseline_criteria <- read_csv(
    file.path(root, 'output', 'inla_hinge_sensitivity', 'full_criteria.csv'),
    show_col_types = FALSE
) |>
    filter(candidate == 'selected_hinge4_plus8') |>
    transmute(
        candidate = 'selected_additive_cots', waic, dic,
        log_marginal_likelihood, mean_neg_log_cpo
    )
criteria <- bind_rows(
    baseline_criteria,
    full_summary$criteria |>
        transmute(
            candidate = 'cots_probability_severity', waic, dic,
            log_marginal_likelihood, mean_neg_log_cpo
        )
)
comparison <- metrics |>
    left_join(criteria, by = 'candidate') |>
    arrange(scheme, rmse)
write_csv(predictions, file.path(audit_dir, 'cv_predictions.csv'))
write_csv(comparison, file.path(audit_dir, 'model_comparison.csv'))

gannett_ids <- data |>
    filter(grepl('Gannet', ReefName, ignore.case = TRUE)) |>
    distinct(source_observation_id) |>
    pull(source_observation_id) |>
    as.character()
gannett_audit <- predictions |>
    filter(
        scheme == 'leave_one_event_out',
        as.character(source_observation_id) %in% gannett_ids
    ) |>
    select(candidate, source_observation_id, ReefID, ReefName, event_year,
           observed_mortality, predicted_mortality) |>
    mutate(residual = observed_mortality - predicted_mortality) |>
    left_join(
        data |>
            transmute(
                source_observation_id = as.character(source_observation_id),
                cot_interval_idw_max, cots_outbreak_probability,
                cots_expected_severe_burden
            ),
        by = 'source_observation_id'
    ) |>
    distinct()
write_csv(gannett_audit, file.path(audit_dir, 'gannett_cay_audit.csv'))

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
            selected_additive_cots = 'Current additive COTS',
            cots_probability_severity = 'Probability + weighted intensity'
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
    geom_text(aes(label = sprintf('%.3f', value)),
              vjust = -0.35, size = 3.1) +
    facet_grid(metric ~ scheme, scales = 'free_y') +
    scale_fill_manual(values = c(
        'Current additive COTS' = '#0072B2',
        'Probability + weighted intensity' = '#D55E00'
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
        title = 'COTS outbreak probability and observed-intensity candidate',
        subtitle = paste0(
            'Intensity is log excess above ', cots_outbreak_threshold,
            ' COTS/tow, weighted by hindcast outbreak probability'
        ),
        x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11)
save_figure_bundle(
    plot, 'Fig-INLA-10_cots_probability_severity', plot_data,
    'Held-out comparison of the current additive COTS predictors and a complementary outbreak-probability plus probability-weighted-intensity representation.',
    'Tests whether realised COTS intensity above the 0.22 COTS/tow outbreak threshold adds complementary information to the hindcast probability signal.',
    'Both terms use the hindcast outbreak definition: probability of density above 0.22 COTS/tow, and probability multiplied by log intensity excess above 0.22 COTS/tow.',
    'cots_probability_severity_candidate', 'INLA', 'model_comparison',
    'operational_candidate_test', root, TRUE,
    code_source = 'src/evaluation/test_inla_cots_hazard.R'
)
write_figure_readme(root)
message('Wrote COTS hazard sensitivity outputs')
