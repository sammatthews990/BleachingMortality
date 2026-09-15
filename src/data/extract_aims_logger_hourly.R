# Preserve short freshwater pulses from the AIMS logger downloads. Daily means
# can conceal acute low-salinity exposure, so this script summarises the GUI's
# hourly product without writing a very large processed hourly table.

suppressPackageStartupMessages({
    library(dplyr)
    library(jsonlite)
    library(lubridate)
    library(readr)
    library(tidyr)
})

raw_dir <- 'data/raw/aims_water_quality'
processed_dir <- 'data/processed'
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

sites_file <- file.path(raw_dir, 'aims_logger_sites.json')
if (!file.exists(sites_file)) {
    download.file(
        paste0(
            'https://data.aims.gov.au/wq-downloads/api/data/',
            'timeseries/sites'
        ),
        sites_file, mode = 'wb', quiet = TRUE
    )
}
sites <- fromJSON(sites_file) |>
    as_tibble() |>
    rename(site_latitude = latitude, site_longitude = longitude)

finite_min <- function(x) if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_
finite_max <- function(x) if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_
finite_quantile <- function(x, probability) {
    if (!any(is.finite(x))) return(NA_real_)
    quantile(x, probability, na.rm = TRUE, names = FALSE)
}
longest_hourly_spell <- function(time, condition) {
    selected <- which(!is.na(condition) & condition)
    if (length(selected) == 0L) return(0L)
    selected_time <- as.numeric(time[selected]) / 3600
    new_spell <- c(TRUE, diff(selected_time) > 1.5)
    max(tabulate(cumsum(new_spell)))
}

event_parts <- vector('list', nrow(sites))
coverage_parts <- vector('list', nrow(sites))
temperature_parts <- vector('list', nrow(sites))

for (site_index in seq_len(nrow(sites))) {
    site_code <- sites$site[[site_index]]
    source_file <- file.path(
        raw_dir, paste0('AIMS_TimeSeries_Data_', site_code, '_Hours.csv')
    )
    if (!file.exists(source_file)) {
        api_url <- paste0(
            'https://data.aims.gov.au/wq-downloads/api/data/timeseries/',
            site_code, '/Hours'
        )
        response <- fromJSON(api_url)
        if (is.null(response$downloadFileUrl)) {
            stop('No AIMS hourly logger URL returned for ', site_code)
        }
        download.file(
            response$downloadFileUrl, source_file, mode = 'wb', quiet = TRUE
        )
    }

    source_lines <- readLines(source_file, warn = FALSE)
    header_line <- which(grepl(
        '^"SITE","SAMPLE_TIME","PARAMETER","VALUE"', source_lines
    ))
    if (length(header_line) != 1L) {
        stop('Could not locate hourly header in ', source_file)
    }
    hourly <- read_csv(
        source_file, skip = header_line - 1L, show_col_types = FALSE
    ) |>
        transmute(
            site = SITE,
            sample_time = dmy_hms(SAMPLE_TIME, tz = 'UTC'),
            parameter = recode(
                PARAMETER,
                'Salinity (PSU)' = 'salinity_psu',
                'Temperature (C)' = 'temperature_c',
                'Conductivity (S/m)' = 'conductivity_sm',
                'Depth (m)' = 'depth_m'
            ),
            value = as.numeric(VALUE)
        ) |>
        pivot_wider(names_from = parameter, values_from = value) |>
        mutate(
            # The CSV clock time is local at the GBR sites even though the
            # download supplies no timezone field. Preserve that clock time
            # when assigning day phases; converting it as UTC shifts noon to
            # night and reverses the nocturnal contrast.
            local_time = force_tz(sample_time, 'Australia/Brisbane'),
            local_date = as.Date(local_time),
            event_year = if_else(
                month(local_time) >= 11L,
                year(local_time) + 1L, year(local_time)
            ),
            in_bleaching_summer = month(local_time) %in% c(11L, 12L, 1:4)
        ) |>
        arrange(sample_time)

    daily_temperature <- hourly |>
        filter(in_bleaching_summer, is.finite(temperature_c)) |>
        group_by(site, event_year, local_date) |>
        summarise(
            observed_temperature_hours = n(),
            daily_mean_c = mean(temperature_c),
            daily_min_c = min(temperature_c),
            daily_max_c = max(temperature_c),
            daily_range_c = daily_max_c - daily_min_c,
            afternoon_max_c = if (
                any(hour(local_time) %in% 12:17)
            ) max(temperature_c[hour(local_time) %in% 12:17]) else NA_real_,
            predawn_min_c = if (
                any(hour(local_time) %in% 0:5)
            ) min(temperature_c[hour(local_time) %in% 0:5]) else NA_real_,
            .groups = 'drop'
        ) |>
        group_by(site, event_year) |>
        arrange(local_date, .by_group = TRUE) |>
        mutate(
            following_predawn_min_c = lead(predawn_min_c),
            following_date = lead(local_date),
            nocturnal_drop_c = if_else(
                following_date == local_date + 1L,
                afternoon_max_c - following_predawn_min_c,
                NA_real_
            )
        ) |>
        ungroup()

    temperature_parts[[site_index]] <- daily_temperature |>
        group_by(site, event_year) |>
        group_modify(~ {
            peak_date <- .x$local_date[which.max(.x$daily_mean_c)]
            peak_window <- .x |>
                filter(
                    local_date >= peak_date - 15L,
                    local_date <= peak_date + 14L
                )
            tibble(
                logger_peak_date = peak_date,
                logger_peak30_dtr_days = sum(is.finite(
                    peak_window$daily_range_c
                )),
                logger_peak30_mean_dtr_c = mean(
                    peak_window$daily_range_c, na.rm = TRUE
                ),
                logger_peak30_median_dtr_c = median(
                    peak_window$daily_range_c, na.rm = TRUE
                ),
                logger_peak30_nocturnal_days = sum(is.finite(
                    peak_window$nocturnal_drop_c
                )),
                logger_peak30_median_nocturnal_drop_c = median(
                    peak_window$nocturnal_drop_c, na.rm = TRUE
                )
            )
        }) |>
        ungroup()

    coverage_parts[[site_index]] <- hourly |>
        summarise(
            site = site_code,
            site_latitude = sites$site_latitude[[site_index]],
            site_longitude = sites$site_longitude[[site_index]],
            first_time = min(sample_time),
            last_time = max(sample_time),
            hourly_rows = n(),
            salinity_hours = sum(is.finite(salinity_psu)),
            temperature_hours = sum(is.finite(temperature_c))
        )

    event_parts[[site_index]] <- hourly |>
        filter(in_bleaching_summer) |>
        group_by(site, event_year) |>
        summarise(
            site_latitude = sites$site_latitude[[site_index]],
            site_longitude = sites$site_longitude[[site_index]],
            first_time = min(sample_time),
            last_time = max(sample_time),
            observed_hours = n(),
            salinity_observed_hours = sum(is.finite(salinity_psu)),
            salinity_min_psu = finite_min(salinity_psu),
            salinity_p05_psu = finite_quantile(salinity_psu, 0.05),
            salinity_p10_psu = finite_quantile(salinity_psu, 0.10),
            salinity_hours_below34 = sum(salinity_psu < 34, na.rm = TRUE),
            salinity_hours_below32 = sum(salinity_psu < 32, na.rm = TRUE),
            salinity_hours_below30 = sum(salinity_psu < 30, na.rm = TRUE),
            salinity_hours_below25 = sum(salinity_psu < 25, na.rm = TRUE),
            salinity_longest_spell_below30_hours = longest_hourly_spell(
                sample_time, salinity_psu < 30
            ),
            freshwater_deficit30_psu_hours = sum(
                pmax(30 - salinity_psu, 0), na.rm = TRUE
            ),
            temperature_min_c = finite_min(temperature_c),
            temperature_max_c = finite_max(temperature_c),
            .groups = 'drop'
        )
}

hourly_event_summary <- bind_rows(event_parts) |>
    left_join(bind_rows(temperature_parts), by = c('site', 'event_year')) |>
    mutate(
        across(
            c(logger_peak30_mean_dtr_c, logger_peak30_median_dtr_c),
            ~ if_else(logger_peak30_dtr_days >= 20L, .x, NA_real_)
        ),
        logger_peak30_median_nocturnal_drop_c = if_else(
            logger_peak30_nocturnal_days >= 15L,
            logger_peak30_median_nocturnal_drop_c, NA_real_
        )
    )
hourly_coverage <- bind_rows(coverage_parts)

write_csv(
    hourly_event_summary,
    file.path(processed_dir, 'aims_wq_logger_hourly_event_summary.csv')
)
write_csv(
    hourly_coverage,
    file.path(processed_dir, 'aims_wq_logger_hourly_coverage.csv')
)

cat('AIMS hourly logger sites:', nrow(hourly_coverage), '\n')
cat('Hourly site-event summaries:', nrow(hourly_event_summary), '\n')
cat(
    'Latest hourly observation:',
    format(max(hourly_coverage$last_time), tz = 'UTC'), '\n'
)
