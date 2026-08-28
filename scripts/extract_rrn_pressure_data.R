# Extract the Reef 2050 Integrated Monitoring and Reporting Program pressure
# layers needed for the mortality analysis. The source workbook is untouched.
# Austral summer 202324 is assigned to event year 2024. Reef-relative water
# colour percentiles use only the preceding ten event years.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(readxl)
    library(tidyr)
})

input_file <- 'data/GBRMPA_RRN_2025_AIMS.xlsx'
output_file <- 'data/processed/rrn_pressure_reef_year.csv'
metadata_file <- 'data/processed/rrn_pressure_metadata.csv'
dynamic_types <- c(
    'cyc_maxHrs4mw', 'wqc_freqcc12',
    'cot_meanpertow', 'cot_idwmeanpertow', 'sst_maxdhw'
)

if (!file.exists(input_file)) stop('Missing RRN workbook: ', input_file)

pressure <- read_excel(
    input_file,
    sheet = '4.IndividualPressureData',
    col_types = c('text', 'text', 'numeric', 'text', 'numeric')
) |>
    filter(Type %in% dynamic_types) |>
    transmute(
        LABEL_ID = trimws(LABEL_ID),
        pressure_type = Type,
        source_summer = as.integer(Year),
        event_year = as.integer(floor(Year / 100) + 1L),
        value = as.numeric(value)
    )

expected_summer <- pressure |>
    distinct(source_summer, event_year) |>
    mutate(expected_summer = (event_year - 1L) * 100L + event_year %% 100L)
if (any(expected_summer$source_summer != expected_summer$expected_summer)) {
    stop('At least one RRN summer code could not be mapped to its ending year')
}

duplicates <- pressure |>
    count(LABEL_ID, event_year, pressure_type) |>
    filter(n != 1L)
if (nrow(duplicates) > 0L) {
    stop('RRN pressure data contain duplicate reef-year-variable records')
}

wide <- pressure |>
    pivot_wider(names_from = pressure_type, values_from = value) |>
    arrange(event_year, LABEL_ID)

coral_sink <- read_excel(
    input_file,
    sheet = '4.IndividualPressureData',
    col_types = c('text', 'text', 'numeric', 'text', 'numeric')
) |>
    filter(Type == 'uq_coralsink1') |>
    transmute(LABEL_ID = trimws(LABEL_ID), uq_coralsink1 = as.numeric(value))
if (anyDuplicated(coral_sink$LABEL_ID)) {
    stop('RRN coral-sink data contain duplicate reef records')
}
wide <- wide |>
    left_join(coral_sink, by = 'LABEL_ID', relationship = 'many-to-one')

add_prior_wq_context <- function(reef_rows) {
    reef_rows <- arrange(reef_rows, event_year)
    reef_rows$wqc_prior10_n <- 0L
    reef_rows$wqc_prior10_percentile <- NA_real_
    reef_rows$wqc_prior10_delta <- NA_real_
    reef_rows$wqc_excess50_10yr_sum <- NA_real_
    for (i in seq_len(nrow(reef_rows))) {
        year <- reef_rows$event_year[[i]]
        current <- reef_rows$wqc_freqcc12[[i]]
        history <- reef_rows$wqc_freqcc12[
            reef_rows$event_year >= year - 10L &
                reef_rows$event_year < year
        ]
        history <- history[is.finite(history)]
        reef_rows$wqc_prior10_n[[i]] <- length(history)
        if (is.finite(current) && length(history) > 0L) {
            reef_rows$wqc_prior10_percentile[[i]] <- (
                sum(history < current) + 0.5 * sum(history == current)
            ) / length(history)
            reef_rows$wqc_prior10_delta[[i]] <- current - median(history)
        }
        inclusive_window <- reef_rows$wqc_freqcc12[
            reef_rows$event_year >= year - 9L &
                reef_rows$event_year <= year
        ]
        inclusive_window <- inclusive_window[is.finite(inclusive_window)]
        if (length(inclusive_window) > 0L) {
            reef_rows$wqc_excess50_10yr_sum[[i]] <- sum(
                pmax(inclusive_window - 0.50, 0)
            )
        }
    }
    reef_rows
}

wide <- wide |>
    group_by(LABEL_ID) |>
    group_modify(~ add_prior_wq_context(.x)) |>
    ungroup() |>
    mutate(
        wqc_excess50 = pmax(wqc_freqcc12 - 0.50, 0),
        log1p_cyc_maxHrs4mw = log1p(pmax(cyc_maxHrs4mw, 0)),
        log1p_cot_idwmeanpertow = log1p(pmax(cot_idwmeanpertow, 0))
    ) |>
    arrange(event_year, LABEL_ID)

metadata <- tibble(
    variable = c(dynamic_types, 'uq_coralsink1'),
    role = c(
        'mechanical cyclone-wave exposure', 'freshwater/plume proxy',
        'observed COTS pressure', 'spatially complete modelled COTS pressure',
        'annual thermal exposure', 'static coral-larval sink connectivity'
    ),
    units = c(
        'maximum hours exposed to 4 m waves during any cyclone that summer',
        'frequency of colour classes 1-5 on cloud-free wet-season days',
        'mean COTS per manta tow', 'IDW-modelled mean COTS per manta tow'
        , 'maximum degree heating weeks', 'normalised index (0-1)'
    ),
    source_period = c(
        'Nov 1 to Apr 30', 'Dec 1 to Apr 30',
        'Jul 1 to Jun 30', 'Jul 1 to Jun 30', 'Nov 1 to Apr 30',
        'static layer from Hock et al. 2017'
    ),
    operational_note = c(
        'Available through summer 202425 in the supplied workbook',
        'Usually released one wet season later; supplied through summer 202324',
        'Missing where no manta observation was available',
        'Continuous reef-scale surface suitable for prediction',
        'Available through summer 202425 in the supplied workbook',
        'Candidate recovery predictor; not an acute mortality mechanism'
    )
)

metadata <- bind_rows(
    metadata,
    tibble(
        variable = c(
            'wqc_prior10_percentile', 'wqc_prior10_delta',
            'wqc_excess50', 'wqc_excess50_10yr_sum',
            'log1p_cyc_maxHrs4mw', 'log1p_cot_idwmeanpertow'
        ),
        role = c(
            'reef-relative freshwater/plume anomaly',
            'reef-relative freshwater/plume anomaly',
            'high coloured-water threshold exposure',
            'cumulative high coloured-water exposure',
            'right-skewed cyclone-wave predictor',
            'right-skewed spatially complete COTS predictor'
        ),
        units = c(
            'mid-rank percentile against preceding ten event years',
            'current frequency minus preceding-ten-year median',
            'frequency above 0.50',
            'sum of annual frequency excess above 0.50',
            'log(1 + hours)', 'log(1 + modelled mean COTS per tow)'
        ),
        source_period = c(
            'preceding ten event years only', 'preceding ten event years only',
            'current wet season', 'current and preceding nine event years',
            'Nov 1 to Apr 30', 'Jul 1 to Jun 30'
        ),
        operational_note = c(
            'Current event is excluded from the reference distribution',
            'Current event is excluded from the reference distribution',
            'Hinge at 0.50 on the native 0-1 WQC scale',
            'Inclusive rolling window; operationally available only when WQC is released',
            'Retains zero exposure and reduces leverage of rare extremes',
            'Preferred prediction feature; observed COTS retained for audit'
        )
    )
)

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(wide, output_file, na = '')
write_csv(metadata, metadata_file)

coverage <- pressure |>
    group_by(pressure_type) |>
    summarise(
        first_summer = min(source_summer),
        last_summer = max(source_summer),
        first_event_year = min(event_year),
        last_event_year = max(event_year),
        reef_year_rows = n(),
        nonmissing_rows = sum(!is.na(value)),
        reefs = n_distinct(LABEL_ID),
        .groups = 'drop'
    )
print(coverage)
