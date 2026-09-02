# Authenticated AIMS temperature extraction with deployment depth retained.
#
# Required environment variable:
#   AIMS_DATAPLATFORM_API_KEY
# Optional controls:
#   AIMS_EVENT_YEARS=2016,2017,2020,2022,2024
#   AIMS_REFRESH_TEMPERATURE=1
#   AIMS_PULL_AUTOMATED=1
#   AIMS_AUTOMATED_SITES=Lizard Island
# Use AIMS_AUTOMATED_SITES=all to pull all GBR automated-station water
# temperatures. The default is deliberately limited to Lizard Island because
# the high-frequency all-GBR weather product is very large.

suppressPackageStartupMessages({
    library(dplyr)
    library(httr2)
    library(jsonlite)
    library(lubridate)
    library(readr)
})

api_key <- Sys.getenv('AIMS_DATAPLATFORM_API_KEY')
if (!nzchar(api_key)) {
    stop(
        'AIMS_DATAPLATFORM_API_KEY is not set. ',
        'Store it in the user environment; never put it in this script.'
    )
}

temp_base <- paste0(
    'https://api.aims.gov.au/data-v2.0/',
    '10.25845/5b4eb0f9bb848'
)
weather_base <- paste0(
    'https://api.aims.gov.au/data-v2.0/',
    '10.25845/5c09bf93f315d'
)
event_years <- as.integer(strsplit(
    Sys.getenv('AIMS_EVENT_YEARS', '2016,2017,2020,2022,2024'),
    ',', fixed = TRUE
)[[1]])
refresh <- identical(Sys.getenv('AIMS_REFRESH_TEMPERATURE'), '1')
pull_automated <- !identical(Sys.getenv('AIMS_PULL_AUTOMATED', '1'), '0')
automated_sites <- trimws(strsplit(
    Sys.getenv('AIMS_AUTOMATED_SITES', 'Lizard Island'),
    ',', fixed = TRUE
)[[1]])

raw_dir <- 'data/raw/aims_temperature/api'
processed_dir <- 'data/processed'
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

make_url <- function(base, parameters = list()) {
    if (length(parameters) == 0L) return(base)
    query <- paste0(
        names(parameters), '=',
        vapply(
            parameters,
            function(value) URLencode(as.character(value), reserved = TRUE),
            character(1)
        ),
        collapse = '&'
    )
    paste0(base, '?', query)
}

get_json <- function(url) {
    request(url) |>
        req_headers(`x-api-key` = api_key) |>
        req_timeout(180) |>
        req_retry(max_tries = 5, backoff = function(tries) 2 ^ tries) |>
        req_perform() |>
        resp_body_json(simplifyVector = TRUE)
}

fetch_pages <- function(url) {
    rows <- list()
    page <- 1L
    seen <- character()
    while (nzchar(url) && !url %in% seen) {
        seen <- c(seen, url)
        payload <- get_json(url)
        result <- payload$results
        if (is.null(result) || NROW(result) == 0L) break
        rows[[page]] <- as_tibble(result)
        next_url <- payload$links[['next']]
        if (is.null(next_url) || !nzchar(next_url)) break
        url <- next_url
        page <- page + 1L
    }
    bind_rows(rows)
}

event_window <- function(event_year) {
    list(
        from_date = paste0(event_year - 1L, '-11-01'),
        thru_date = paste0(event_year, '-05-31')
    )
}

deployment_file <- file.path(raw_dir, 'temperature_logger_deployments.csv.gz')
if (refresh || !file.exists(deployment_file)) {
    deployments <- get_json(paste0(
        temp_base, '/data/summary-by-deployment?size=10000'
    )) |>
        as_tibble()
    write_csv(deployments, deployment_file)
} else {
    deployments <- read_csv(deployment_file, show_col_types = FALSE)
}
deployments <- deployments |>
    mutate(
        deployment_id = as.character(deployment_id),
        depth_m = as.numeric(depth),
        depth_class = case_when(
            is.na(depth_m) ~ 'unknown',
            depth_m <= 5 ~ 'shallow_0_5m',
            TRUE ~ 'deep_gt5m'
        )
    )

logger_parts <- vector('list', length(event_years))
for (i in seq_along(event_years)) {
    event_year <- event_years[[i]]
    raw_file <- file.path(
        raw_dir, paste0('temperature_logger_daily_', event_year, '.csv.gz')
    )
    if (refresh || !file.exists(raw_file)) {
        dates <- event_window(event_year)
        url <- make_url(
            paste0(temp_base, '/data/daily'),
            c(dates, list(size = 10000))
        )
        event_rows <- fetch_pages(url)
        write_csv(event_rows, raw_file)
    } else {
        event_rows <- read_csv(raw_file, show_col_types = FALSE)
    }
    logger_parts[[i]] <- event_rows |>
        mutate(
            event_year = event_year,
            deployment_id = as.character(deployment_id),
            date = as.Date(time),
            temperature_c = coalesce(qc_val, cal_val)
        ) |>
        left_join(
            deployments |>
                select(deployment_id, depth_m, depth_class),
            by = 'deployment_id', relationship = 'many-to-one'
        ) |>
        transmute(
            source = 'temperature_logger', event_year, date,
            site, subsite, series, series_id, deployment_id,
            lat, lon, depth_m, depth_class,
            temperature_c, daily_observations = qc_count
        )
    message('Temperature logger event ', event_year, ': ', nrow(event_rows))
}

automated_parts <- list()
if (pull_automated) {
    weather_series <- get_json(paste0(weather_base, '/series')) |>
        as_tibble() |>
        filter(grepl('Water Temperature', series, fixed = TRUE))
    write_csv(
        weather_series,
        file.path(raw_dir, 'automated_water_temperature_series.csv.gz')
    )

    part <- 1L
    for (event_year in event_years) {
        dates <- event_window(event_year)
        requested_sites <- if (
            length(automated_sites) == 1L &&
                tolower(automated_sites[[1]]) == 'all'
        ) NA_character_ else automated_sites
        if (all(is.na(requested_sites))) requested_sites <- NA_character_

        for (site_filter in requested_sites) {
            site_label <- if_else(
                is.na(site_filter), 'all_gbr',
                gsub('[^A-Za-z0-9]+', '_', tolower(site_filter))
            )
            raw_file <- file.path(
                raw_dir,
                paste0('automated_temperature_', site_label, '_', event_year, '.csv.gz')
            )
            if (refresh || !file.exists(raw_file)) {
                parameters <- c(
                    dates,
                    list(
                        parameter = 'Water Temperature',
                        min_lon = 142, max_lon = 154,
                        min_lat = -25, max_lat = -10,
                        size = 10000
                    )
                )
                if (!is.na(site_filter)) parameters$site <- site_filter
                event_rows <- fetch_pages(make_url(
                    paste0(weather_base, '/data'), parameters
                ))
                write_csv(event_rows, raw_file)
            } else {
                event_rows <- read_csv(raw_file, show_col_types = FALSE)
            }
            if (nrow(event_rows) == 0L) next

            automated_parts[[part]] <- event_rows |>
                mutate(
                    event_year = event_year,
                    # The endpoint mixes ISO datetimes and already-aggregated
                    # ISO dates. The first 10 characters are the common,
                    # timezone-independent calendar date in both cases.
                    date = as.Date(substr(as.character(time), 1, 10)),
                    temperature_c = coalesce(qc_val, raw_val),
                    depth_m = as.numeric(depth),
                    depth_class = case_when(
                        is.na(depth_m) ~ 'unknown',
                        depth_m <= 5 ~ 'shallow_0_5m',
                        TRUE ~ 'deep_gt5m'
                    )
                ) |>
                group_by(
                    event_year, date, site, subsite, series, series_id,
                    lat, lon, depth_m, depth_class
                ) |>
                summarise(
                    daily_observations = sum(is.finite(temperature_c)),
                    temperature_c = mean(temperature_c, na.rm = TRUE),
                    .groups = 'drop'
                ) |>
                mutate(
                    source = 'automated_weather',
                    deployment_id = NA_character_,
                    .before = 1
                )
            part <- part + 1L
            message(
                'Automated temperature ', site_label, ' event ',
                event_year, ': ', nrow(event_rows)
            )
        }
    }
}

combined <- bind_rows(bind_rows(logger_parts), bind_rows(automated_parts)) |>
    filter(
        between(lon, 142, 154), between(lat, -25, -10),
        is.finite(temperature_c)
    ) |>
    arrange(source, event_year, site, series, date)

write_csv(
    combined,
    file.path(processed_dir, 'aims_temperature_daily_with_depth.csv.gz')
)
write_csv(
    combined |>
        group_by(
            source, event_year, site, series, series_id,
            lat, lon, depth_m, depth_class
        ) |>
        summarise(
            first_date = min(date), last_date = max(date),
            observed_days = n_distinct(date),
            mean_temperature_c = mean(temperature_c),
            max_temperature_c = max(temperature_c),
            .groups = 'drop'
        ),
    file.path(processed_dir, 'aims_temperature_depth_coverage.csv')
)

cat('Combined daily rows:', nrow(combined), '\n')
cat('Series:', n_distinct(combined$series), '\n')
cat('Sites:', n_distinct(combined$site), '\n')
cat('Depth classes:\n')
print(count(combined, source, depth_class))
