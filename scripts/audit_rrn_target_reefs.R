# Extract every RRN pressure and metric record for five severe 2024 manta-tow
# underpredictions, then benchmark valid 2023-24 pressures against regional,
# GBR-wide, and within-reef historical distributions.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(readxl)
    library(tidyr)
})

input_file <- 'data/GBRMPA_RRN_2025_AIMS.xlsx'
output_dir <- 'output/rrn_target_reefs'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

target_reefs <- tibble(
    LABEL_ID = c('16-015', '14-126', '14-116a', '14-118', '15-030'),
    reef_short = c('Mackay', 'Linnet', 'Lizard NW', 'Eyrie', 'Swinger')
)

reef_reference <- read_csv(
    'data/AIMS-Reef_Reference.csv', show_col_types = FALSE
) |>
    transmute(
        LABEL_ID = ReefID, ReefName, SECT_NAME, SECTOR,
        region_block = case_when(
            SECTOR %in% c('CG', 'CL', 'PB') ~ 'Northern GBR',
            SECTOR %in% c('CA', 'CU', 'IN', 'TO', 'WH') ~ 'Central GBR',
            SECTOR %in% c('CB', 'PO', 'SW') ~ 'Southern GBR',
            TRUE ~ NA_character_
        )
    ) |>
    distinct(LABEL_ID, .keep_all = TRUE)

pressure <- read_excel(
    input_file, sheet = '4.IndividualPressureData',
    col_types = c('text', 'text', 'numeric', 'text', 'numeric')
) |>
    mutate(
        LABEL_ID = trimws(LABEL_ID),
        latitude_band = suppressWarnings(as.integer(substr(LABEL_ID, 1, 2))),
        source_summer = as.integer(Year),
        event_year = if_else(
            is.na(Year), NA_integer_, as.integer(floor(Year / 100) + 1L)
        )
    ) |>
    left_join(reef_reference, by = 'LABEL_ID', relationship = 'many-to-one')

metrics <- read_excel(
    input_file, sheet = '5.MetricsData',
    col_types = c('text', 'text', 'text', 'numeric', 'text', 'numeric')
) |>
    mutate(
        LABEL_ID = trimws(LABEL_ID),
        latitude_band = suppressWarnings(as.integer(substr(LABEL_ID, 1, 2)))
    ) |>
    left_join(reef_reference, by = 'LABEL_ID', relationship = 'many-to-one') |>
    mutate(
        includes_summer_202425 = case_when(
            Type == 'cmex' & (
                (variable == 'cmex_1y' & Year == 202424) |
                (variable == 'cmex_2y' & Year == 202324) |
                (variable == 'cmex_10y' & Year == 201524) |
                (variable == 'cmex_20y' & Year == 200524)
            ) ~ TRUE,
            Type %in% c('cmfi', 'quad') & (
                (grepl('_1y$', variable) & Year == 202324) |
                (grepl('_2y$', variable) & Year == 202224) |
                (grepl('_10y$', variable) & Year == 201424) |
                (grepl('_20y$', variable) & Year == 200424)
            ) ~ TRUE,
            TRUE ~ FALSE
        )
    )

target_pressure_all <- pressure |>
    semi_join(target_reefs, by = 'LABEL_ID') |>
    left_join(target_reefs, by = 'LABEL_ID') |>
    select(
        reef_short, LABEL_ID, ReefName, SECT_NAME, SECTOR, region_block,
        latitude_band, Type, source_summer, event_year, variable, value
    ) |>
    arrange(reef_short, Type, source_summer)
write_csv(
    target_pressure_all,
    file.path(output_dir, 'target_reef_all_individual_pressures.csv'),
    na = ''
)

target_metrics_all <- metrics |>
    semi_join(target_reefs, by = 'LABEL_ID') |>
    left_join(target_reefs, by = 'LABEL_ID') |>
    select(
        reef_short, LABEL_ID, ReefName, SECT_NAME, SECTOR, region_block,
        latitude_band, Type, Period, Year, variable, value,
        includes_summer_202425
    ) |>
    arrange(reef_short, Type, variable, Year)
write_csv(
    target_metrics_all,
    file.path(output_dir, 'target_reef_all_metrics.csv'),
    na = ''
)

midrank_percent <- function(x) {
    n <- sum(!is.na(x))
    if (n == 0L) return(rep(NA_real_, length(x)))
    100 * (rank(x, ties.method = 'average', na.last = 'keep') - 0.5) / n
}

peer_pressure <- pressure |>
    filter(
        (source_summer == 202324 & Type != 'uq_coralsink1') |
            Type == 'uq_coralsink1'
    ) |>
    group_by(Type, source_summer) |>
    mutate(gbr_percentile = midrank_percent(value)) |>
    ungroup() |>
    group_by(Type, source_summer, latitude_band) |>
    mutate(latitude_band_percentile = midrank_percent(value)) |>
    ungroup() |>
    group_by(Type, source_summer, region_block) |>
    mutate(region_percentile = midrank_percent(value)) |>
    ungroup()

history_10y <- pressure |>
    filter(
        LABEL_ID %in% target_reefs$LABEL_ID,
        Type != 'uq_coralsink1',
        source_summer >= 201415,
        source_summer <= 202324
    ) |>
    group_by(LABEL_ID, Type) |>
    mutate(within_reef_10y_percentile = midrank_percent(value)) |>
    summarise(
        target_202324 = value[source_summer == 202324][1],
        within_reef_10y_percentile =
            within_reef_10y_percentile[source_summer == 202324][1],
        previous_10y_median = median(value[source_summer < 202324], na.rm = TRUE),
        previous_10y_q90 = quantile(
            value[source_summer < 202324], 0.9, na.rm = TRUE
        ),
        .groups = 'drop'
    )

pressure_profile <- peer_pressure |>
    semi_join(target_reefs, by = 'LABEL_ID') |>
    left_join(target_reefs, by = 'LABEL_ID') |>
    left_join(history_10y, by = c('LABEL_ID', 'Type')) |>
    mutate(
        vulnerability_percentile = if_else(
            Type == 'uq_coralsink1', 100 - latitude_band_percentile,
            latitude_band_percentile
        ),
        interpretation_direction = case_when(
            Type == 'uq_coralsink1' ~
                'Lower coral-sink connectivity may reduce recovery/rescue potential',
            Type == 'wqc_freqcc12' ~
                'Higher coloured-water frequency may indicate plume/freshwater exposure',
            Type == 'cyc_maxHrs4mw' ~
                'Higher values indicate mechanical cyclone-wave exposure',
            Type %in% c('cot_meanpertow', 'cot_idwmeanpertow') ~
                'Higher values indicate greater COTS pressure',
            Type == 'sst_maxdhw' ~
                'Higher values indicate greater thermal exposure',
            TRUE ~ NA_character_
        )
    ) |>
    select(
        reef_short, LABEL_ID, ReefName, region_block, SECT_NAME,
        latitude_band, Type, source_summer, value,
        latitude_band_percentile, region_percentile, gbr_percentile,
        within_reef_10y_percentile, previous_10y_median, previous_10y_q90,
        vulnerability_percentile, interpretation_direction
    ) |>
    arrange(reef_short, desc(vulnerability_percentile))
write_csv(
    pressure_profile,
    file.path(output_dir, 'target_reef_2024_pressure_percentiles.csv'),
    na = ''
)

metric_peers <- metrics |>
    group_by(Type, variable, Year) |>
    mutate(gbr_percentile = midrank_percent(value)) |>
    ungroup() |>
    group_by(Type, variable, Year, latitude_band) |>
    mutate(latitude_band_percentile = midrank_percent(value)) |>
    ungroup() |>
    group_by(Type, variable, Year, region_block) |>
    mutate(region_percentile = midrank_percent(value)) |>
    ungroup()

target_metric_percentiles <- metric_peers |>
    semi_join(target_reefs, by = 'LABEL_ID') |>
    left_join(target_reefs, by = 'LABEL_ID') |>
    select(
        reef_short, LABEL_ID, ReefName, region_block, latitude_band,
        Type, Period, Year, variable, value, latitude_band_percentile,
        region_percentile, gbr_percentile,
        includes_summer_202425
    ) |>
    arrange(reef_short, Type, variable, Year)
write_csv(
    target_metric_percentiles,
    file.path(output_dir, 'target_reef_metric_percentiles.csv'),
    na = ''
)

mortality <- read_csv(
    'output/rrn_pressure_assessment/severe_2024_reef_diagnostics.csv',
    show_col_types = FALSE
) |>
    filter(
        scheme == 'reef_blocked_2024',
        programme_key == 'manta',
        ReefID %in% target_reefs$LABEL_ID
    ) |>
    distinct(ReefID, .keep_all = TRUE) |>
    transmute(
        LABEL_ID = ReefID, observed_mortality, ann_maxdhw,
        prop_acropora_pre, core_underprediction, weather_underprediction,
        wq_underprediction, cyclone_underprediction,
        DISTURBANCE_TYPE, storm_name, description, tooltip, disturbance_text
    )

valid_2024_wide <- pressure_profile |>
    select(LABEL_ID, Type, value) |>
    pivot_wider(names_from = Type, values_from = value)

valid_2024_percentiles <- pressure_profile |>
    select(LABEL_ID, Type, vulnerability_percentile) |>
    pivot_wider(
        names_from = Type, values_from = vulnerability_percentile,
        names_prefix = 'vulnerability_pct_'
    )

reef_context <- target_reefs |>
    left_join(reef_reference, by = 'LABEL_ID') |>
    left_join(mortality, by = 'LABEL_ID') |>
    left_join(valid_2024_wide, by = 'LABEL_ID') |>
    left_join(valid_2024_percentiles, by = 'LABEL_ID') |>
    arrange(desc(observed_mortality))
write_csv(
    reef_context,
    file.path(output_dir, 'target_reef_2024_context.csv'),
    na = ''
)

manta_2024 <- read_csv(
    'output/rrn_pressure_assessment/mortality_brt_paired_rows.csv',
    show_col_types = FALSE
) |>
    filter(
        scheme == 'reef_blocked_2024',
        programme_key == 'manta',
        event_year == 2024L
    ) |>
    distinct(ReefID, source_observation_id, .keep_all = TRUE) |>
    mutate(core_residual = observed_mortality - predicted_mortality_core) |>
    select(-any_of(sort(unique(peer_pressure$Type))))

all_2024_pressures <- peer_pressure |>
    select(LABEL_ID, Type, value) |>
    pivot_wider(names_from = Type, values_from = value)

manta_pressure_rows <- manta_2024 |>
    left_join(
        all_2024_pressures,
        by = c('ReefID' = 'LABEL_ID'),
        relationship = 'many-to-one'
    )

manta_associations <- manta_pressure_rows |>
    select(
        ReefID, observed_mortality, core_residual,
        all_of(sort(unique(peer_pressure$Type)))
    ) |>
    pivot_longer(
        all_of(sort(unique(peer_pressure$Type))),
        names_to = 'Type', values_to = 'value'
    ) |>
    filter(!is.na(value)) |>
    group_by(Type) |>
    summarise(
        observations = n(),
        positive_values = sum(value > 0),
        spearman_observed_mortality = cor(
            value, observed_mortality, method = 'spearman'
        ),
        spearman_core_residual = cor(
            value, core_residual, method = 'spearman'
        ),
        upper_quartile_threshold = quantile(value, 0.75),
        mean_residual_upper_quartile = mean(
            core_residual[value >= quantile(value, 0.75)]
        ),
        mean_residual_other = mean(
            core_residual[value < quantile(value, 0.75)]
        ),
        residual_difference = mean_residual_upper_quartile -
            mean_residual_other,
        .groups = 'drop'
    )
write_csv(
    manta_associations,
    file.path(output_dir, 'manta_2024_rrn_residual_associations.csv'),
    na = ''
)

coverage <- tibble(
    item = c(
        'RRN reefs in individual pressure data',
        'RRN reefs matched to AIMS reference',
        'Target reefs requested',
        'Target reefs present in pressure data',
        'Target reefs present in metrics data'
    ),
    value = c(
        n_distinct(pressure$LABEL_ID),
        n_distinct(pressure$LABEL_ID[!is.na(pressure$region_block)]),
        nrow(target_reefs),
        n_distinct(target_pressure_all$LABEL_ID),
        n_distinct(target_metrics_all$LABEL_ID)
    )
)
write_csv(coverage, file.path(output_dir, 'extraction_coverage.csv'))

print(coverage)
print(
    pressure_profile |>
        select(
            reef_short, Type, value, region_percentile,
            latitude_band_percentile, within_reef_10y_percentile,
            vulnerability_percentile
        ),
    n = Inf
)
print(manta_associations)
