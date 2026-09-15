# Link AIMS wet-season water-quality observations to the nearest GBR feature
# and assemble the same predictors for every reef. AIMS observations remain
# point measurements: distance to the nearest reef is retained so river-mouth
# samples cannot be mistaken for direct reef salinity observations.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(sf)
})

processed_dir <- 'data/processed'
reef_file <- paste0(
    'data/Great_Barrier_Reef_Features/',
    'Great_Barrier_Reef_Features.shp'
)

required_files <- c(
    aims = file.path(processed_dir, 'aims_wq_shallow_event_summary.csv'),
    rrn = file.path(processed_dir, 'rrn_pressure_reef_year.csv'),
    era5 = file.path(processed_dir, 'era5_weather_reef_year.csv'),
    ereefs = file.path(
        processed_dir, 'ereefs_salinity_imos_kd490_reef_year.csv'
    ),
    cyclone = file.path(processed_dir, 'bom_cyclone_reef_year.csv'),
    reefs = reef_file
)
if (!all(file.exists(required_files))) {
    stop('Run the AIMS, ERA5, eReefs, RRN and cyclone extractors first.')
}

reef_attributes <- st_read(reef_file, quiet = TRUE) |>
    st_drop_geometry() |>
    distinct(LABEL_ID, .keep_all = TRUE) |>
    transmute(
        ReefID = LABEL_ID, ReefName = LOC_NAME_S,
        reef_longitude = X_COORD, reef_latitude = Y_COORD
    )
# X_COORD/Y_COORD are the registry's representative reef points. They avoid
# treating multipart reef polygons as independent observations.
reef_points <- st_as_sf(
    reef_attributes,
    coords = c('reef_longitude', 'reef_latitude'), crs = 4283,
    remove = FALSE
) |>
    st_transform(3577)
reef_coordinates <- st_transform(reef_points, 4326) |>
    st_drop_geometry()

aims <- read_csv(required_files[['aims']], show_col_types = FALSE)
normalise_name <- function(x) {
    x |>
        toupper() |>
        gsub('[^A-Z0-9]+', ' ', x = _) |>
        trimws()
}

# Prefer unique survey-name lineage over spatial proximity. This prevents
# nearby multipart reefs (for example Linnet and Eyrie) being assigned to an
# adjacent registry feature solely because its representative point is closer.
transition_file <- file.path(processed_dir, 'annual_coral_transitions.csv')
survey_crosswalk <- read_csv(transition_file, show_col_types = FALSE) |>
    transmute(
        aims_name_key = normalise_name(raw_reef_name_key),
        name_match_ReefID = ReefID,
        name_match_ReefName = ReefName
    ) |>
    filter(nzchar(aims_name_key)) |>
    group_by(aims_name_key) |>
    filter(n_distinct(name_match_ReefID) == 1L) |>
    slice(1L) |>
    ungroup()

aims <- aims |>
    mutate(aims_name_key = normalise_name(LOCATION_NAME)) |>
    left_join(survey_crosswalk, by = 'aims_name_key', relationship = 'many-to-one')
aims_points <- st_as_sf(
    aims, coords = c('LONGITUDE', 'LATITUDE'), crs = 4326,
    remove = FALSE
) |>
    st_transform(3577)
nearest_index <- st_nearest_feature(aims_points, reef_points)
name_match_index <- match(aims$name_match_ReefID, reef_points$ReefID)
selected_index <- ifelse(
    is.finite(name_match_index), name_match_index, nearest_index
)
nearest_distance_km <- as.numeric(st_distance(
    aims_points, reef_points[selected_index, ], by_element = TRUE
)) / 1000

aims_linked <- aims |>
    mutate(
        ReefID = reef_points$ReefID[selected_index],
        ReefName = if_else(
            is.finite(name_match_index),
            name_match_ReefName, reef_points$ReefName[selected_index]
        ),
        reef_match_method = if_else(
            is.finite(name_match_index),
            'unique survey name', 'nearest registry point'
        ),
        nearest_reef_distance_km = nearest_distance_km,
        reef_proximal_5km = nearest_reef_distance_km <= 5,
        reef_proximal_20km = nearest_reef_distance_km <= 20,
        observation_depth_class = if_else(
            sample_depth_max_m <= 1,
            'surface/plume (<=1 m)', 'coral depth (1-10 m)'
        ),
        low_salinity_below30 = salinity_min_psu < 30,
        severe_low_salinity_below25 = salinity_min_psu < 25
    ) |>
    select(-name_match_ReefID, -name_match_ReefName)

rrn <- read_csv(required_files[['rrn']], show_col_types = FALSE) |>
    transmute(
        ReefID = LABEL_ID, event_year,
        rrn_coloured_water = wqc_freqcc12,
        rrn_coloured_water_decadal_percentile = wqc_prior10_percentile,
        rrn_cyclone_wave_hours4m = cyc_maxHrs4mw,
        rrn_cots_idw_per_tow = cot_idwmeanpertow
    )
era5 <- read_csv(required_files[['era5']], show_col_types = FALSE) |>
    transmute(
        ReefID = LABEL_ID, event_year = year,
        era5_coastal_rain_dec_mar_max30 =
            era5_coastal_rain_dec_mar_max_30day,
        era5_coastal_rain_q1_total = era5_coastal_rain_q1_total,
        era5_coastal_rain_q1_max7 = era5_coastal_rain_q1_max_7day,
        era5_reef_rain_q1_total = era5_reef_rain_q1_total,
        era5_reef_wind_q1_mean = era5_reef_wind_q1_mean
    ) |>
    group_by(ReefID, event_year) |>
    summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)),
              .groups = 'drop')
ereefs <- read_csv(required_files[['ereefs']], show_col_types = FALSE) |>
    transmute(
        ReefID = LABEL_ID, event_year = year,
        ereefs_salinity_min = salinity_min,
        ereefs_freshwater_exposure30 = freshwater_exposure_30,
        ereefs_days_below30 = days_below_30,
        imos_kd490_q90 = k490_q90,
        imos_secchi_p10 = secc3m_p10
    ) |>
    group_by(ReefID, event_year) |>
    summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)),
              .groups = 'drop')
cyclone <- read_csv(required_files[['cyclone']], show_col_types = FALSE) |>
    select(
        ReefID, event_year, tc_min_distance_km, tc_nearest_name,
        tc_nearest_max_wind_ms, tc_max_wind_within300km_ms,
        tc_min_pressure_within300km_hpa, tc_storms_within300km,
        tc_inside_reported_gale_radius
    )

event_years <- sort(unique(aims$event_year))
predictor_grid <- tidyr::crossing(
    reef_coordinates, event_year = event_years
) |>
    left_join(rrn, by = c('ReefID', 'event_year')) |>
    left_join(era5, by = c('ReefID', 'event_year')) |>
    left_join(ereefs, by = c('ReefID', 'event_year')) |>
    left_join(cyclone, by = c('ReefID', 'event_year'))

# The WMIP pilot is an optional acquisition sensitivity. Its March snapshot
# matches the existing Q1 calibration window. Absence of the file leaves the
# established calibration grid unchanged; it is not a production prerequisite.
wmip_file <- file.path(processed_dir, 'freshwater_reef_event.csv')
if (file.exists(wmip_file)) {
    wmip <- read_csv(wmip_file, show_col_types = FALSE) |>
        filter(
            cutoff_name == 'march',
            product_mode %in% c('initial_forecast', 'environmental_hindcast')
        ) |>
        select(
            ReefID, event_year,
            wmip_issue_time_utc = issue_time_utc,
            wmip_routed_discharge_total_ml,
            wmip_routed_discharge_max7_m3_s,
            wmip_connection_weight_sum,
            wmip_nearest_source_distance_km,
            wmip_dominant_source_id,
            wmip_routing_method,
            wmip_use_in_primary_mortality_model =
                use_in_primary_mortality_model
        )
    if (anyDuplicated(wmip[c('ReefID', 'event_year')])) {
        stop('WMIP March pilot has duplicate reef-event keys')
    }
    predictor_grid <- predictor_grid |>
        left_join(wmip, by = c('ReefID', 'event_year'),
                  relationship = 'one-to-one')
    message('Joined the optional WMIP March freshwater-source pilot')
}

training_points <- aims_linked |>
    left_join(
        predictor_grid,
        by = c('ReefID', 'ReefName', 'event_year'),
        relationship = 'many-to-one'
    )

write_csv(
    training_points,
    file.path(processed_dir, 'aims_freshwater_training_points.csv')
)
write_csv(
    predictor_grid,
    file.path(processed_dir, 'freshwater_predictor_grid_2016_2025.csv')
)

cat('AIMS point-event summaries:', nrow(training_points), '\n')
cat('Within 5 km of a reef:', sum(training_points$reef_proximal_5km), '\n')
cat('Within 20 km of a reef:', sum(training_points$reef_proximal_20km), '\n')
cat('GBR reef-event prediction rows:', nrow(predictor_grid), '\n')
cat(
    'Events with complete ERA5/eReefs predictor layers:',
    paste(sort(unique(predictor_grid$event_year[
        is.finite(predictor_grid$ereefs_salinity_min) &
            is.finite(predictor_grid$era5_coastal_rain_dec_mar_max30)
    ])), collapse = ', '), '\n'
)
