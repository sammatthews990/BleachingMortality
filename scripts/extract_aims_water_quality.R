# Reproducible extraction of the public AIMS water-quality portal. The portal
# exposes continuous logger files and a separate bounding-box sample download.

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

refresh_raw <- identical(Sys.getenv('REFRESH_AIMS_WQ'), '1')
sites_url <- paste0(
    'https://data.aims.gov.au/wq-downloads/api/data/',
    'timeseries/sites'
)

sites_file <- file.path(raw_dir, 'aims_logger_sites.json')
if (refresh_raw || !file.exists(sites_file)) {
    download.file(sites_url, sites_file, mode = 'wb', quiet = TRUE)
}
sites <- fromJSON(sites_file) |>
    as_tibble() |>
    rename(site_latitude = latitude, site_longitude = longitude)

for (site_code in sites$site) {
    destination <- file.path(
        raw_dir, paste0('AIMS_TimeSeries_Data_', site_code, '_Days.csv')
    )
    if (refresh_raw || !file.exists(destination)) {
        link_url <- paste0(
            'https://data.aims.gov.au/wq-downloads/api/data/timeseries/',
            site_code, '/Days'
        )
        response <- fromJSON(link_url)
        if (is.null(response$downloadFileUrl)) {
            stop('No AIMS logger download URL returned for ', site_code)
        }
        download.file(
            response$downloadFileUrl, destination, mode = 'wb', quiet = TRUE
        )
    }
}

logger_rows <- vector('list', nrow(sites))
for (i in seq_len(nrow(sites))) {
    site_code <- sites$site[[i]]
    source_file <- file.path(
        raw_dir, paste0('AIMS_TimeSeries_Data_', site_code, '_Days.csv')
    )
    source_lines <- readLines(source_file, warn = FALSE)
    header_line <- which(grepl(
        '^"SITE","SAMPLE_TIME","PARAMETER","VALUE"', source_lines
    ))
    if (length(header_line) != 1L) {
        stop('Could not locate the AIMS logger header in ', source_file)
    }
    logger_rows[[i]] <- read_csv(
        source_file, skip = header_line - 1L, show_col_types = FALSE
    ) |>
        mutate(
            sample_date = dmy_hms(SAMPLE_TIME, tz = 'UTC'),
            VALUE = as.numeric(VALUE)
        ) |>
        select(site = SITE, sample_date, parameter = PARAMETER, value = VALUE)
}

logger_daily <- bind_rows(logger_rows) |>
    left_join(sites, by = 'site', relationship = 'many-to-one') |>
    arrange(site, sample_date, parameter)

longest_run <- function(condition) {
    condition[is.na(condition)] <- FALSE
    runs <- rle(condition)
    if (!any(runs$values)) return(0L)
    max(runs$lengths[runs$values])
}
finite_min <- function(x) if (any(is.finite(x))) min(x, na.rm = TRUE) else NA_real_
finite_max <- function(x) if (any(is.finite(x))) max(x, na.rm = TRUE) else NA_real_
finite_quantile <- function(x, probability) {
    if (!any(is.finite(x))) return(NA_real_)
    quantile(x, probability, na.rm = TRUE, names = FALSE)
}

logger_wide <- logger_daily |>
    mutate(
        parameter = recode(
            parameter,
            'Salinity (PSU)' = 'salinity_psu',
            'Temperature (C)' = 'temperature_c',
            'Conductivity (S/m)' = 'conductivity_sm'
        )
    ) |>
    pivot_wider(names_from = parameter, values_from = value) |>
    mutate(
        event_year = if_else(
            month(sample_date) >= 11L,
            year(sample_date) + 1L, year(sample_date)
        ),
        in_bleaching_summer = month(sample_date) %in% c(11L, 12L, 1:4)
    )

logger_event_summary <- logger_wide |>
    filter(in_bleaching_summer) |>
    arrange(site, event_year, sample_date) |>
    group_by(site, site_latitude, site_longitude, event_year) |>
    summarise(
        first_date = min(sample_date),
        last_date = max(sample_date),
        observed_days = n_distinct(as.Date(sample_date)),
        salinity_observed_days = sum(is.finite(salinity_psu)),
        salinity_min_psu = finite_min(salinity_psu),
        salinity_p10_psu = finite_quantile(salinity_psu, 0.1),
        salinity_days_below30 = sum(salinity_psu < 30, na.rm = TRUE),
        salinity_days_below25 = sum(salinity_psu < 25, na.rm = TRUE),
        salinity_longest_spell_below30 = longest_run(salinity_psu < 30),
        freshwater_deficit30_psu_days = sum(
            pmax(30 - salinity_psu, 0), na.rm = TRUE
        ),
        temperature_min_c = finite_min(temperature_c),
        temperature_max_c = finite_max(temperature_c),
        .groups = 'drop'
    )

# The GUI has two controls that materially change the discrete download:
# MMP-only versus all programmes, and depth average (da) versus depth-weighted
# average (dwa). Use all programmes so reef-based LTMP-DIC, MBM and legacy
# observations are not silently excluded. A slightly wider western bound gives
# these requests a distinct server filename from the old MMP-only cache.
discrete_event_years <- 2016:2025

read_discrete_archive <- function(discrete_event_year, aggregation) {
    nutrient_from <- paste0(discrete_event_year - 1L, '-11-01')
    nutrient_thru <- paste0(discrete_event_year, '-05-31')
    criterion <- paste0(
        '141.9000/154.0000/-25.0000/-10.0000/from/', nutrient_from,
        '/thru/', nutrient_thru, '/', aggregation, '?mmpOnly=false'
    )
    nutrient_api <- paste0(
        'https://data.aims.gov.au/wq-downloads/api/data/nut/bbox/',
        criterion
    )
    nutrient_zip <- file.path(
        raw_dir, paste0(
            'AIMS_NutData_GBR_event', discrete_event_year,
            '_all_', aggregation, '.zip'
        )
    )
    archive_name <- paste0(
        'AIMS_NutData_14191540-250-100',
        gsub('-', '', nutrient_from), gsub('-', '', nutrient_thru),
        aggregation, '.zip'
    )
    archive_url <- paste0(
        'https://data.aims.gov.au/data-download/rwqpp/wqnut/', archive_name
    )

    if (refresh_raw || !file.exists(nutrient_zip)) {
        direct_status <- try(
            download.file(
                archive_url, nutrient_zip, mode = 'wb', quiet = TRUE
            ),
            silent = TRUE
        )
        if (inherits(direct_status, 'try-error')) {
            response <- NULL
            for (attempt in 1:3) {
                response <- tryCatch(
                    fromJSON(nutrient_api), error = function(error) NULL
                )
                if (!is.null(response$downloadFileUrl)) break
                Sys.sleep(attempt)
            }
            if (is.null(response$downloadFileUrl)) {
                stop(
                    'No AIMS ', aggregation, ' sample URL for event ',
                    discrete_event_year
                )
            }
            download.file(
                response$downloadFileUrl, nutrient_zip,
                mode = 'wb', quiet = TRUE
            )
        }
    }

    archive_files <- unzip(nutrient_zip, list = TRUE)$Name
    if (length(archive_files) != 1L ||
        !grepl('[.]csv$', archive_files, ignore.case = TRUE)) {
        stop('Unexpected AIMS archive contents for event ', discrete_event_year)
    }
    nutrient_lines <- readLines(
        unz(nutrient_zip, archive_files), warn = FALSE
    )
    nutrient_header <- which(grepl(
        '^"STATION_NAME","LOCATION_NAME","SHORT_NAME"', nutrient_lines
    ))
    if (length(nutrient_header) != 1L) {
        stop('No AIMS sample header for event ', discrete_event_year)
    }
    read_csv(
        I(paste(
            nutrient_lines[(nutrient_header):length(nutrient_lines)],
            collapse = '\n'
        )),
        show_col_types = FALSE
    ) |>
        mutate(download_aggregation = aggregation)
}

depth_average_parts <- lapply(
    discrete_event_years, read_discrete_archive, aggregation = 'da'
)
depth_weighted_parts <- lapply(
    discrete_event_years, read_discrete_archive, aggregation = 'dwa'
)

water_samples <- bind_rows(depth_average_parts) |>
    mutate(
        sample_datetime = dmy_hms(COLLECTION_START_DATE, tz = 'UTC'),
        event_year = if_else(
            month(sample_datetime) >= 11L,
            year(sample_datetime) + 1L, year(sample_datetime)
        ),
        is_mmp = startsWith(coalesce(PROJECT, ''), 'MMP-'),
        low_salinity_below30 = SAL < 30,
        low_salinity_below25 = SAL < 25
    ) |>
    arrange(sample_datetime, STATION_NAME, SAMPLE_DEPTH)

depth_weighted_samples <- bind_rows(depth_weighted_parts) |>
    mutate(
        sample_datetime = dmy_hms(COLLECTION_START_DATE, tz = 'UTC'),
        event_year = if_else(
            month(sample_datetime) >= 11L,
            year(sample_datetime) + 1L, year(sample_datetime)
        ),
        is_mmp = startsWith(coalesce(PROJECT, ''), 'MMP-'),
        low_salinity_below30 = SAL < 30,
        low_salinity_below25 = SAL < 25
    ) |>
    arrange(sample_datetime, STATION_NAME)

surface_event_summary <- water_samples |>
    filter(is.finite(SAMPLE_DEPTH), SAMPLE_DEPTH <= 1) |>
    group_by(
        event_year, PROJECT, SHORT_NAME, LOCATION_NAME,
        LATITUDE, LONGITUDE
    ) |>
    summarise(
        sample_dates = n_distinct(as.Date(sample_datetime)),
        sample_depth_min_m = finite_min(SAMPLE_DEPTH),
        sample_depth_max_m = finite_max(SAMPLE_DEPTH),
        salinity_min_psu = finite_min(SAL),
        salinity_below30_samples = sum(SAL < 30, na.rm = TRUE),
        salinity_below25_samples = sum(SAL < 25, na.rm = TRUE),
        suspended_solids_max = finite_max(SS),
        cdom_max = finite_max(CDOM),
        secchi_min_m = finite_min(SECCHI_DEPTH),
        chlorophyll_max = finite_max(CHL),
        .groups = 'drop'
    )

# Reef monitoring samples are commonly collected at 5 m. Retain a separate
# coral-depth table rather than discarding them with the <=1 m plume samples.
shallow_event_summary <- water_samples |>
    filter(is.finite(SAMPLE_DEPTH), SAMPLE_DEPTH <= 10) |>
    group_by(
        event_year, PROJECT, SHORT_NAME, LOCATION_NAME,
        LATITUDE, LONGITUDE
    ) |>
    summarise(
        sample_dates = n_distinct(as.Date(sample_datetime)),
        sample_depth_min_m = finite_min(SAMPLE_DEPTH),
        sample_depth_max_m = finite_max(SAMPLE_DEPTH),
        salinity_min_psu = finite_min(SAL),
        salinity_below30_samples = sum(SAL < 30, na.rm = TRUE),
        salinity_below25_samples = sum(SAL < 25, na.rm = TRUE),
        suspended_solids_max = finite_max(SS),
        cdom_max = finite_max(CDOM),
        secchi_min_m = finite_min(SECCHI_DEPTH),
        chlorophyll_max = finite_max(CHL),
        .groups = 'drop'
    )

depth_weighted_event_summary <- depth_weighted_samples |>
    group_by(
        event_year, PROJECT, SHORT_NAME, LOCATION_NAME,
        LATITUDE, LONGITUDE
    ) |>
    summarise(
        sample_dates = n_distinct(as.Date(sample_datetime)),
        salinity_min_psu = finite_min(SAL),
        salinity_below30_samples = sum(SAL < 30, na.rm = TRUE),
        salinity_below25_samples = sum(SAL < 25, na.rm = TRUE),
        suspended_solids_max = finite_max(SS),
        cdom_max = finite_max(CDOM),
        secchi_min_m = finite_min(SECCHI_DEPTH),
        chlorophyll_max = finite_max(CHL),
        .groups = 'drop'
    )

project_coverage <- water_samples |>
    mutate(PROJECT = coalesce(PROJECT, 'Unspecified')) |>
    group_by(event_year, PROJECT) |>
    summarise(
        rows = n(),
        stations = n_distinct(STATION_NAME),
        locations = n_distinct(LOCATION_NAME),
        salinity_rows = sum(is.finite(SAL)),
        temperature_rows = sum(is.finite(TEMP)),
        suspended_solids_rows = sum(is.finite(SS)),
        cdom_rows = sum(is.finite(CDOM)),
        secchi_rows = sum(is.finite(SECCHI_DEPTH)),
        .groups = 'drop'
    )

station_coverage <- logger_wide |>
    group_by(site, site_latitude, site_longitude) |>
    summarise(
        first_date = min(sample_date), last_date = max(sample_date),
        daily_rows = n(),
        salinity_days = sum(is.finite(salinity_psu)),
        temperature_days = sum(is.finite(temperature_c)),
        .groups = 'drop'
    )

write_csv(sites, file.path(processed_dir, 'aims_wq_logger_sites.csv'))
write_csv(logger_wide, file.path(processed_dir, 'aims_wq_logger_daily.csv'))
write_csv(
    logger_event_summary,
    file.path(processed_dir, 'aims_wq_logger_event_summary.csv')
)
write_csv(
    water_samples,
    file.path(processed_dir, 'aims_wq_discrete_samples_2015_2025.csv')
)
write_csv(
    depth_weighted_samples,
    file.path(
        processed_dir,
        'aims_wq_discrete_depth_weighted_2015_2025.csv'
    )
)
write_csv(
    surface_event_summary,
    file.path(processed_dir, 'aims_wq_surface_event_summary.csv')
)
write_csv(
    shallow_event_summary,
    file.path(processed_dir, 'aims_wq_shallow_event_summary.csv')
)
write_csv(
    depth_weighted_event_summary,
    file.path(processed_dir, 'aims_wq_depth_weighted_event_summary.csv')
)
write_csv(
    project_coverage,
    file.path(processed_dir, 'aims_wq_project_coverage.csv')
)
write_csv(
    station_coverage,
    file.path(processed_dir, 'aims_wq_station_coverage.csv')
)

cat('AIMS logger sites:', nrow(sites), '\n')
cat('Logger daily site-date rows:', nrow(logger_wide), '\n')
cat('All-program depth-average rows:', nrow(water_samples), '\n')
cat('All-program depth-weighted rows:', nrow(depth_weighted_samples), '\n')
cat('Shallow (<=10 m) station-event summaries:', nrow(shallow_event_summary), '\n')
cat('Discrete records below 30 PSU:', sum(water_samples$SAL < 30, na.rm = TRUE), '\n')
