# Build the current major-miss register from the selected operational model.
# Repeated programme rows are summarised to reef-event level, while retaining
# their count and maximum absolute observation-level residual.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(stringr)
})
source('scripts/model_registry.R')

root <- project_root()
output_dir <- file.path(root, 'output', 'model_reporting')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
selected_id <- read_model_registry(root)$best_model$id
max_or_na <- function(x) {
    if (all(!is.finite(x))) NA_real_ else max(x[is.finite(x)])
}

prediction_path <- if (selected_id %in% c(
    'prospective_cots_log_peak_recency',
    'operational_rrn_raw_plus_manta_state'
)) file.path(
    root, 'output', 'prospective_cots_nowcast', 'cv_predictions.csv'
) else if (identical(
    selected_id, 'cots_raw_interval_logistic20_5_cyclone'
)) file.path(
    root, 'output', 'cots_raw_enso_occurrence',
    'cots_cv_predictions.csv'
) else if (identical(
    selected_id, 'baseline_cots_logistic20_5_cyclone'
)) file.path(
    root, 'output', 'cots_timing_cover_cyclone_soft_gate',
    'cv_predictions.csv'
) else if (identical(
    selected_id, 'cause_aware_exposure_gated_competing_hazards'
)) file.path(
    root, 'output', 'cause_aware_competing_hazards',
    'gated_cv_predictions.csv'
) else file.path(root, 'output', 'inla_spatiotemporal', 'cv_predictions.csv')

prediction_context <- read_csv(
    file.path(root, 'output', 'explanatory_event_dhw',
              'event_dhw_brt_data.csv'),
    show_col_types = FALSE
) |>
    mutate(source_observation_id = as.character(source_observation_id)) |>
    distinct(source_observation_id, .keep_all = TRUE) |>
    select(
        source_observation_id, ann_maxdhw, cyc_interval_maxHrs4mw,
        cot_interval_idw_max, cots_outbreak_probability,
        wqc_freqcc12, wqc_prior10_percentile
    )

prediction_id <- if (identical(
    selected_id, 'cots_raw_interval_logistic20_5_cyclone'
)) 'cots_raw_interval_relative' else selected_id
predictions <- read_csv(prediction_path, show_col_types = FALSE) |>
    filter(
        candidate == .env$prediction_id,
        scheme == 'leave_one_event_out'
    ) |>
    mutate(source_observation_id = as.character(source_observation_id)) |>
    select(-any_of(names(prediction_context)[-1])) |>
    left_join(prediction_context, by = 'source_observation_id') |>
    mutate(residual = observed_mortality - predicted_mortality)

reef_events <- predictions |>
    group_by(ReefID, ReefName, event_year) |>
    summarise(
        programmes = paste(sort(unique(programme_key)), collapse = ' + '),
        observations = n(),
        observed_mortality = mean(observed_mortality),
        predicted_mortality = mean(predicted_mortality),
        mean_residual = mean(residual),
        max_absolute_residual = max(abs(residual)),
        max_dhw = mean(ann_maxdhw),
        cyclone_wave_hours = max_or_na(cyc_interval_maxHrs4mw),
        cots_interval_max = max_or_na(cot_interval_idw_max),
        cots_outbreak_probability = max_or_na(cots_outbreak_probability),
        current_wqc = max_or_na(wqc_freqcc12),
        relative_wqc = max_or_na(wqc_prior10_percentile),
        .groups = 'drop'
    ) |>
    arrange(desc(max_absolute_residual)) |>
    mutate(current_rank = row_number())

curated <- read_csv(
    file.path(root, 'data', 'curated', 'extreme_miss_explanations.csv'),
    show_col_types = FALSE
) |>
    mutate(event_year = as.integer(event_year))

context <- read_csv(
    file.path(root, 'output', 'explanatory_event_dhw',
              'event_dhw_brt_data.csv'),
    show_col_types = FALSE
) |>
    group_by(ReefID, event_year) |>
    summarise(
        survey_context = paste(
            unique(na.omit(str_squish(
                disturbance_text[!is.na(disturbance_text) &
                    str_squish(disturbance_text) != '']
            ))),
            collapse = ' | '
        ),
        .groups = 'drop'
    )

major_misses <- reef_events |>
    slice_head(n = 20) |>
    left_join(
        curated,
        by = c('ReefID' = 'reef_id', 'event_year')
    ) |>
    left_join(context, by = c('ReefID', 'event_year')) |>
    mutate(
        direction = if_else(
            mean_residual >= 0, 'Underpredicted', 'Overpredicted'
        ),
        note_has_other_disturbance = str_detect(
            str_to_lower(coalesce(survey_context, '')),
            'cyclone|jasper|flood|storm|crown-of-thorns|cots|disease|physical damage'
        ),
        explanation_status = case_when(
            str_detect(
                str_to_lower(coalesce(evidence_status, '')),
                'high confidence|high ecological confidence|high contextual confidence'
            ) ~ 'Mechanism supported; magnitude unresolved',
            !is.na(primary_explanation) ~ 'Partially explained',
            note_has_other_disturbance ~
                'Potential explanation; needs validation',
            TRUE ~ 'Unexplained'
        ),
        current_explanation = case_when(
            !is.na(primary_explanation) ~ primary_explanation,
            note_has_other_disturbance ~ paste0(
                'Other disturbance appears in survey record: ',
                str_trunc(survey_context, 125)
            ),
            TRUE ~ 'No supported explanation currently registered'
        ),
        still_unexplained = explanation_status %in% c(
            'Potential explanation; needs validation', 'Unexplained'
        )
    ) |>
    select(
        current_rank, ReefID, ReefName, event_year, programmes, observations,
        observed_mortality, predicted_mortality, mean_residual,
        max_absolute_residual, direction, explanation_status,
        still_unexplained, current_explanation, evidence_status,
        current_numeric_signal, pending_check, survey_context,
        max_dhw, cyclone_wave_hours, cots_interval_max,
        cots_outbreak_probability, current_wqc, relative_wqc
    )

registered <- reef_events |>
    inner_join(
        curated,
        by = c('ReefID' = 'reef_id', 'event_year')
    ) |>
    arrange(current_rank)

write_csv(
    major_misses,
    file.path(output_dir, 'major_misses_current_selected.csv')
)
write_csv(
    registered,
    file.path(output_dir, 'registered_miss_explanations_current.csv')
)
write_csv(
    major_misses |>
        count(explanation_status, still_unexplained, name = 'reef_events'),
    file.path(output_dir, 'major_miss_status_summary.csv')
)

message('Wrote current selected-model major-miss tables')
