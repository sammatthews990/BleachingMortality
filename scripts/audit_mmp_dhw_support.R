# Explain why the reconstructed MMP modelling field has no values below 4 DHW.
# The legacy survey-linked DHW and the independently reconstructed NOAA field
# are retained side by side; no threshold is used to filter rows.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

output_dir <- 'output/zero_mortality_sensitivity'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outcomes <- readRDS('data/processed/mortality_outcomes_mmp.rds')
validation <- readRDS('data/processed/validation_rows_mmp.rds')

if (nrow(outcomes) != nrow(validation)) {
    stop('MMP row count changed between outcome and validation tables.')
}

row_audit <- validation |>
    transmute(
        source_observation_id, ReefID, ReefName, event_year,
        survey_date, depth, reef_zone, reefpage_category,
        mortality_prop,
        legacy_survey_max_dhw,
        reconstructed_ann_maxdhw = ann_maxdhw,
        dhw_difference = reconstructed_ann_maxdhw - legacy_survey_max_dhw,
        legacy_below4 = legacy_survey_max_dhw < 4,
        reconstructed_below4 = reconstructed_ann_maxdhw < 4,
        legacy_below4_reclassified = legacy_below4 & !reconstructed_below4,
        environmental_feature_count, environment_match_method,
        lon, lat
    ) |>
    arrange(reconstructed_ann_maxdhw, ReefID, source_observation_id)

stage_summary <- bind_rows(
    outcomes |>
        transmute(event_year, dhw = MaxDHW.mean) |>
        mutate(field = 'Legacy survey-linked MaxDHW.mean'),
    validation |>
        transmute(event_year, dhw = ann_maxdhw) |>
        mutate(field = 'Reconstructed NOAA ann_maxdhw')
) |>
    group_by(field, event_year) |>
    summarise(
        rows = n(),
        minimum = min(dhw),
        q10 = quantile(dhw, 0.10),
        median = median(dhw),
        maximum = max(dhw),
        rows_below4 = sum(dhw < 4),
        rows_below4_5 = sum(dhw < 4.5),
        .groups = 'drop'
    )

near_threshold <- row_audit |>
    filter(
        legacy_survey_max_dhw < 4.5 |
            reconstructed_ann_maxdhw < 4.5
    )

write_csv(row_audit, file.path(output_dir, 'mmp_dhw_row_audit.csv'), na = '')
write_csv(stage_summary, file.path(output_dir, 'mmp_dhw_stage_summary.csv'), na = '')
write_csv(near_threshold, file.path(output_dir, 'mmp_dhw_near4_rows.csv'), na = '')

print(stage_summary)
print(near_threshold)
