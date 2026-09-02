# Operationally available COTS timing/cover support and cyclone activation
# sensitivity. Historical disturbance labels supervise cause-model training;
# they are never used as predictors for held-out or future reef-events.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(stringr)
    library(tidyr)
})

source('src/models/fit_cause_aware_competing_hazards.R')

assessment_dir <- file.path(
    root, 'output', 'cots_timing_cover_cyclone_soft_gate'
)
dir.create(assessment_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Leakage-safe COTS history available at the start of each event year
# -------------------------------------------------------------------------

cots_history <- read_csv(
    file.path(root, 'data', 'processed', 'rrn_pressure_reef_year.csv'),
    show_col_types = FALSE
) |>
    transmute(
        ReefID = LABEL_ID,
        event_year = as.integer(event_year),
        cots_density_year = as.numeric(cot_idwmeanpertow),
        cots_burden_year = log1p(pmax(cots_density_year - 0.22, 0))
    ) |>
    arrange(ReefID, event_year) |>
    group_by(ReefID) |>
    mutate(
        cots_prior3_burden =
            coalesce(lag(cots_burden_year, 1), 0) +
            coalesce(lag(cots_burden_year, 2), 0) +
            coalesce(lag(cots_burden_year, 3), 0),
        latest_outbreak_year = cummax(if_else(
            coalesce(cots_density_year, 0) > 0.22,
            event_year, -Inf
        )),
        prior_outbreak_year = lag(latest_outbreak_year),
        cots_years_since_prior_outbreak = pmin(
            event_year - prior_outbreak_year, 10
        ),
        cots_years_since_prior_outbreak = if_else(
            is.finite(cots_years_since_prior_outbreak),
            cots_years_since_prior_outbreak, 10
        )
    ) |>
    ungroup() |>
    select(
        ReefID, event_year, cots_density_year, cots_prior3_burden,
        cots_years_since_prior_outbreak
    )

add_cots_support <- function(rows) {
    rows |>
        left_join(
            cots_history,
            by = c('ReefID', 'event_year'), relationship = 'many-to-one'
        ) |>
        mutate(
            cots_prior3_burden = coalesce(cots_prior3_burden, 0),
            cots_years_since_prior_outbreak = coalesce(
                cots_years_since_prior_outbreak, 10
            ),
            cots_pressure_x_cover = cots_weighted_excess *
                sqrt(pmax(pre_cover, 0)),
            cots_prior_burden_x_cover = cots_prior3_burden *
                sqrt(pmax(pre_cover, 0))
        )
}

annual_supported <- add_cots_support(annual)

support_features <- c(
    'cots_probability', 'cots_weighted_excess', 'pre_cover',
    'cots_pressure_x_cover', 'cots_prior3_burden',
    'cots_prior_burden_x_cover', 'cots_years_since_prior_outbreak',
    'acropora', 'lon', 'lat', 'is_manta', 'is_mmp'
)

# Correlation is diagnostic only; hierarchy is retained for the two pressure x
# cover interactions and candidate promotion is based on held-out prediction.
correlation_rows <- annual_supported |>
    select(all_of(support_features)) |>
    mutate(across(everything(), as.numeric))
support_correlation <- cor(
    correlation_rows, use = 'pairwise.complete.obs', method = 'spearman'
)
write_csv(
    as.data.frame(as.table(support_correlation), stringsAsFactors = FALSE) |>
        as_tibble() |>
        rename(feature_1 = Var1, feature_2 = Var2,
               spearman_correlation = Freq),
    file.path(assessment_dir, 'cots_support_feature_correlations.csv')
)

# -------------------------------------------------------------------------
# 2. Cross-validated COTS timing/cover predictions
# -------------------------------------------------------------------------

supported_cots_predictions <- tibble()

for (scheme in c('leave_one_event_out', 'reef_blocked_5fold')) {
    folds <- if (scheme == 'leave_one_event_out') {
        event_years
    } else sort(unique(data$joint_reef_fold))

    for (fold in folds) {
        if (scheme == 'leave_one_event_out') {
            assessment <- data |> filter(event_year == fold)
            annual_training <- annual_supported |>
                filter(!(baseline_report_year <= fold & report_year >= fold))
        } else {
            held_reefs <- unique(data$ReefID[data$joint_reef_fold == fold])
            assessment <- data |> filter(ReefID %in% held_reefs)
            annual_training <- annual_supported |>
                filter(!ReefID %in% held_reefs)
        }

        assessment_features <- make_assessment_features(assessment) |>
            add_cots_support()
        fitted <- fit_cause_hurdle(
            annual_training, assessment_features, 'cots', support_features
        )
        supported_cots_predictions <- bind_rows(
            supported_cots_predictions,
            assessment_features |>
                transmute(
                    source_observation_id,
                    supported_cots_prediction = fitted$prediction,
                    supported_cots_occurrence = fitted$occurrence,
                    cots_density_year, cots_prior3_burden,
                    cots_years_since_prior_outbreak,
                    cots_pressure_x_cover, cots_prior_burden_x_cover,
                    scheme, fold = as.character(fold)
                )
        )
    }
}

write_csv(
    supported_cots_predictions,
    file.path(assessment_dir, 'supported_cots_predictions.csv')
)

# -------------------------------------------------------------------------
# 3. Factorial COTS support x cyclone activation comparison
# -------------------------------------------------------------------------

activation_functions <- function(wave_hours) {
    wave_hours <- pmax(coalesce(wave_hours, 0), 0)
    tibble(
        hard20 = as.numeric(wave_hours > 20),
        ramp5_20 = pmin(pmax((wave_hours - 5) / 15, 0), 1),
        logistic20_5 = plogis((wave_hours - 20) / 5)
    )
}

base_rows <- gated_predictions |>
    left_join(
        supported_cots_predictions,
        by = c('source_observation_id', 'scheme', 'fold'),
        relationship = 'one-to-one'
    )

make_combined <- function(rows, candidate_name, use_supported_cots,
                          cyclone_activation) {
    activations <- activation_functions(rows$cyclone_wave_hours)
    cyclone_gate_value <- activations[[cyclone_activation]]
    cots_prediction_value <- if (use_supported_cots) {
        rows$supported_cots_prediction
    } else rows$cots_prediction
    cots_occurrence_value <- if (use_supported_cots) {
        rows$supported_cots_occurrence
    } else rows$cots_occurrence
    rows |>
        mutate(
            cots_prediction_test = cots_prediction_value,
            cots_occurrence_test = cots_occurrence_value,
            cots_activation = as.numeric(coalesce(cots_interval_max, 0) > 0.22),
            cyclone_activation = cyclone_gate_value,
            predicted_mortality = 1 -
                (1 - thermal_prediction) *
                (1 - cots_prediction_test * cots_activation) *
                (1 - cyclone_prediction * cyclone_activation),
            predicted_occurrence = 1 -
                (1 - thermal_occurrence) *
                (1 - cots_occurrence_test * cots_activation) *
                (1 - cyclone_occurrence * cyclone_activation),
            residual = observed_mortality - predicted_mortality,
            candidate = candidate_name
        )
}

factorial_predictions <- bind_rows(
    make_combined(
        base_rows, 'cause_aware_exposure_gated_competing_hazards',
        FALSE, 'hard20'
    ),
    make_combined(
        base_rows, 'cots_timing_cover_hard20_cyclone', TRUE, 'hard20'
    ),
    make_combined(
        base_rows, 'baseline_cots_ramp5_20_cyclone', FALSE, 'ramp5_20'
    ),
    make_combined(
        base_rows, 'baseline_cots_logistic20_5_cyclone',
        FALSE, 'logistic20_5'
    ),
    make_combined(
        base_rows, 'cots_timing_cover_ramp5_20_cyclone', TRUE, 'ramp5_20'
    ),
    make_combined(
        base_rows, 'cots_timing_cover_logistic20_5_cyclone',
        TRUE, 'logistic20_5'
    )
)

comparison <- factorial_predictions |>
    group_by(candidate, scheme) |>
    metric_summary() |>
    arrange(scheme, rmse)
event_metrics <- factorial_predictions |>
    filter(scheme == 'leave_one_event_out') |>
    group_by(candidate, event_year) |>
    metric_summary()

write_csv(
    factorial_predictions,
    file.path(assessment_dir, 'cv_predictions.csv')
)
write_csv(comparison, file.path(assessment_dir, 'model_comparison.csv'))
write_csv(event_metrics, file.path(assessment_dir, 'event_metrics.csv'))

# Full-data supported-COTS coefficient summary.
full_assessment_supported <- make_assessment_features(data) |>
    add_cots_support()
full_supported_cots <- fit_cause_hurdle(
    annual_supported, full_assessment_supported, 'cots', support_features
)
write_csv(
    full_supported_cots$fixed,
    file.path(assessment_dir, 'supported_cots_fixed_effects.csv')
)

# -------------------------------------------------------------------------
# 4. Severe and La Nina event residual audits
# -------------------------------------------------------------------------

best_candidate <- comparison |>
    filter(scheme == 'leave_one_event_out') |>
    arrange(rmse, severe_rmse, false_extreme_rate) |>
    slice(1) |>
    pull(candidate)

selected_rows <- factorial_predictions |>
    filter(
        candidate == best_candidate,
        scheme == 'leave_one_event_out'
    )

model_context <- read_csv(
    file.path(root, 'output', 'explanatory_event_dhw',
              'event_dhw_brt_data.csv'),
    show_col_types = FALSE
) |>
    mutate(source_observation_id = as.character(source_observation_id)) |>
    distinct(source_observation_id, .keep_all = TRUE) |>
    select(
        source_observation_id, observed_pre_cover, prop_acropora_pre,
        ann_maxdhw, applied_dhw_uplift, log_coastal_rain30,
        wqc_freqcc12, wqc_prior10_percentile, cloudp_90, mcur_90,
        cyc_interval_maxHrs4mw, tc_interval_wind_distance_index,
        cot_interval_idw_max, cots_outbreak_probability,
        dhw10_load4, dhw_novelty10, dhw_events_since2016_n6,
        dhw_years_since_last_n6_capped8, disturbance_text
    )

audit_rows <- selected_rows |>
    mutate(source_observation_id = as.character(source_observation_id)) |>
    select(-any_of(names(model_context)[-1])) |>
    left_join(model_context, by = 'source_observation_id') |>
    mutate(
        squared_error = residual^2,
        direction = if_else(residual >= 0, 'Underpredicted', 'Overpredicted')
    )

severe_audit <- audit_rows |>
    filter(observed_mortality >= 0.5 | predicted_mortality >= 0.5) |>
    arrange(desc(abs(residual)))
write_csv(severe_audit, file.path(assessment_dir, 'severe_residual_audit.csv'))

lanina_audit <- audit_rows |>
    filter(event_year %in% c(2020, 2022)) |>
    group_by(event_year) |>
    mutate(
        event_sse = sum(squared_error),
        sse_share = squared_error / event_sse,
        observed_event_mean = mean(observed_mortality),
        observed_event_sd = sd(observed_mortality)
    ) |>
    ungroup() |>
    arrange(event_year, desc(squared_error))
write_csv(lanina_audit, file.path(assessment_dir, 'lanina_residual_audit.csv'))

lanina_reef_event <- lanina_audit |>
    group_by(ReefID, ReefName, event_year) |>
    summarise(
        programmes = paste(sort(unique(programme_key)), collapse = ' + '),
        observations = n(),
        observed_mortality = mean(observed_mortality),
        predicted_mortality = mean(predicted_mortality),
        residual = mean(residual),
        squared_error = sum(squared_error),
        sse_share = sum(sse_share),
        ann_maxdhw = mean(ann_maxdhw),
        pre_cover = mean(observed_pre_cover),
        acropora = mean(prop_acropora_pre, na.rm = TRUE),
        rainfall = mean(log_coastal_rain30, na.rm = TRUE),
        wqc = mean(wqc_freqcc12, na.rm = TRUE),
        cloud = mean(cloudp_90, na.rm = TRUE),
        current = mean(mcur_90, na.rm = TRUE),
        cots = max(cot_interval_idw_max, na.rm = TRUE),
        cyclone_wave_hours = max(cyc_interval_maxHrs4mw, na.rm = TRUE),
        disturbance_text = paste(unique(na.omit(disturbance_text)),
                                 collapse = ' | '),
        .groups = 'drop'
    ) |>
    mutate(
        direction = if_else(residual >= 0, 'Underpredicted', 'Overpredicted')
    ) |>
    arrange(event_year, desc(squared_error))
write_csv(
    lanina_reef_event,
    file.path(assessment_dir, 'lanina_reef_event_contributors.csv')
)

event_distribution <- audit_rows |>
    group_by(event_year) |>
    summarise(
        observations = n(),
        observed_mean = mean(observed_mortality),
        observed_sd = sd(observed_mortality),
        observed_max = max(observed_mortality),
        predicted_mean = mean(predicted_mortality),
        rmse = sqrt(mean(squared_error)),
        mean_benchmark_rmse = sqrt(mean(
            (observed_mortality - mean(observed_mortality))^2
        )),
        predictive_r2 = 1 - sum(squared_error) / sum(
            (observed_mortality - mean(observed_mortality))^2
        ),
        .groups = 'drop'
    )
write_csv(
    event_distribution,
    file.path(assessment_dir, 'event_distribution_diagnostics.csv')
)

# -------------------------------------------------------------------------
# 5. Report/paper figures
# -------------------------------------------------------------------------

comparison_plot_data <- comparison |>
    filter(scheme == 'leave_one_event_out') |>
    select(candidate, rmse, severe_rmse, predictive_r2, false_extreme_rate) |>
    pivot_longer(-candidate, names_to = 'metric', values_to = 'value') |>
    mutate(
        candidate = recode(
            candidate,
            cause_aware_exposure_gated_competing_hazards = 'Current hard gates',
            cots_timing_cover_hard20_cyclone = 'COTS timing/cover',
            baseline_cots_ramp5_20_cyclone = 'Cyclone ramp',
            baseline_cots_logistic20_5_cyclone = 'Cyclone logistic',
            cots_timing_cover_ramp5_20_cyclone = 'Timing/cover + ramp',
            cots_timing_cover_logistic20_5_cyclone = 'Timing/cover + logistic'
        ),
        metric = recode(
            metric, rmse = 'RMSE', severe_rmse = 'Severe RMSE',
            predictive_r2 = 'Predictive R-squared',
            false_extreme_rate = 'False-extreme rate'
        )
    )
comparison_plot <- ggplot(
    comparison_plot_data, aes(candidate, value, fill = candidate)
) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf('%.3f', value)), vjust = -0.2, size = 3) +
    facet_wrap(~ metric, scales = 'free_y', ncol = 2) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = 'COTS timing/cover and cyclone activation sensitivity',
        subtitle = 'Leave-one-event-out validation', x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_figure_bundle(
    comparison_plot, 'Fig-INLA-13_cots_timing_cover_validation',
    comparison_plot_data,
    'Leave-one-event-out performance for COTS timing/starting-cover support and hard, ramped or logistic cyclone activation candidates.',
    'Tests whether operational COTS history and a continuous cyclone activation improve event transfer without using future disturbance labels.',
    'The current cyclone exposure layer is provisional; the activation shape must be retested with the forthcoming track/intensity product.',
    best_candidate, 'INLA composite', 'model_comparison',
    'operational_candidate_test', root, TRUE,
    code_source = 'src/evaluation/test_cots_timing_cover_and_cyclone_soft_gate.R',
    width = 11, height = 8
)

activation_grid <- tibble(wave_hours = seq(0, 50, by = 0.25))
activation_plot_data <- bind_cols(
    activation_grid, activation_functions(activation_grid$wave_hours)
) |>
    pivot_longer(-wave_hours, names_to = 'activation', values_to = 'weight') |>
    mutate(activation = recode(
        activation, hard20 = 'Hard >20 h',
        ramp5_20 = 'Ramp 5-20 h', logistic20_5 = 'Logistic centre 20 h'
    ))
activation_plot <- ggplot(
    activation_plot_data, aes(wave_hours, weight, colour = activation)
) +
    geom_line(linewidth = 1.1) +
    geom_vline(xintercept = 20, linetype = 2, colour = 'grey45') +
    labs(
        title = 'Candidate cyclone-hazard activation functions',
        x = 'Hours with waves above 4 m',
        y = 'Fraction of predicted cyclone hazard activated', colour = NULL
    ) +
    theme_bw(base_size = 11) + theme(legend.position = 'bottom')
save_figure_bundle(
    activation_plot, 'Fig-DATA-03_cyclone_activation_functions',
    activation_plot_data,
    'Hard and continuous cyclone-hazard activation functions compared around the provisional 20-hour damaging-wave reference.',
    'Makes the operational consequence of the cyclone gate explicit before the improved exposure layer is available.',
    'These are deterministic sensitivity functions and do not yet propagate cyclone-exposure uncertainty.',
    best_candidate, 'data', 'activation_function',
    'operational_candidate_test', root, TRUE,
    code_source = 'src/evaluation/test_cots_timing_cover_and_cyclone_soft_gate.R'
)

lanina_plot_data <- lanina_reef_event |>
    group_by(event_year, direction) |>
    slice_max(squared_error, n = 8, with_ties = FALSE) |>
    ungroup() |>
    mutate(
        label = str_remove(ReefName, ' \\([0-9].*$'),
        signed_sse = if_else(direction == 'Underpredicted',
                             squared_error, -squared_error)
    )
lanina_plot <- ggplot(
    lanina_plot_data,
    aes(reorder(label, signed_sse), signed_sse, fill = direction)
) +
    geom_col() + coord_flip() + facet_wrap(~ event_year, scales = 'free_y') +
    scale_fill_manual(values = c(
        Underpredicted = '#B2182B', Overpredicted = '#2166AC'
    )) +
    labs(
        title = 'Largest reef-event contributors to 2020/2022 prediction error',
        subtitle = 'Bar magnitude is summed squared error; sign indicates residual direction',
        x = NULL, y = 'Signed squared-error contribution', fill = NULL
    ) + theme_bw(base_size = 11) + theme(legend.position = 'bottom')
save_figure_bundle(
    lanina_plot, 'Fig-INLA-14_lanina_residual_contributors',
    lanina_plot_data,
    'Largest positive and negative reef-event contributors to squared prediction error in the 2020 and 2022 La Nina bleaching events.',
    'Identifies the small number of reefs driving negative event-specific predictive R-squared and separates underprediction from overprediction.',
    'Residual attribution is hypothesis-generating; correlated environmental covariates and observation error prevent causal assignment.',
    best_candidate, 'INLA composite', 'residual_contributors',
    'operational_diagnostic', root, TRUE,
    code_source = 'src/evaluation/test_cots_timing_cover_and_cyclone_soft_gate.R',
    width = 11, height = 8
)

severe_plot_data <- severe_audit |>
    slice_head(n = 20) |>
    mutate(label = paste0(str_remove(ReefName, ' \\([0-9].*$'), ' ', event_year))
severe_plot <- ggplot(
    severe_plot_data,
    aes(predicted_mortality, observed_mortality, colour = factor(event_year))
) +
    geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey45') +
    geom_segment(
        aes(xend = predicted_mortality, yend = predicted_mortality),
        colour = 'grey70'
    ) +
    geom_point(size = 2.8) +
    ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
        title = 'Largest residuals involving observed or predicted severe loss',
        x = 'Held-out predicted mortality', y = 'Observed mortality',
        colour = 'Event'
    ) + theme_bw(base_size = 11)
save_figure_bundle(
    severe_plot, 'Fig-INLA-15_severe_residual_audit', severe_plot_data,
    'Held-out observed and predicted mortality for the largest residuals among observations with observed or predicted mortality of at least 50%.',
    'Shows which severe misses remain after COTS/cyclone cause separation and timing/activation refinement.',
    'Repeated programme/depth rows can represent the same reef-event and should be interpreted with the reef-event audit table.',
    best_candidate, 'INLA composite', 'severe_residual_audit',
    'operational_diagnostic', root, TRUE,
    code_source = 'src/evaluation/test_cots_timing_cover_and_cyclone_soft_gate.R',
    width = 9, height = 8
)

writeLines(best_candidate, file.path(assessment_dir, 'best_candidate.txt'))
write_figure_readme(root)
message('Wrote COTS timing/cover, cyclone activation and La Nina audits')
