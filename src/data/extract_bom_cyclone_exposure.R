# Convert the Bureau of Meteorology Australian cyclone best-track database to
# reef-by-summer exposure fields. Track proximity complements, rather than
# replaces, RRN 4 m wave-hours and rainfall/freshwater predictors.

suppressPackageStartupMessages({
    library(dplyr)
    library(lubridate)
    library(matrixStats)
    library(readr)
    library(sf)
    library(units)
})

source_file <- 'data/raw/bom_cyclones/IDCKMSTM0S.csv'
reef_file <- 'data/Great_Barrier_Reef_Features/Great_Barrier_Reef_Features.shp'
output_file <- 'data/processed/bom_cyclone_reef_year.csv'
track_output <- 'data/processed/bom_cyclone_tracks_gbr.csv'
metadata_file <- 'data/processed/bom_cyclone_exposure_metadata.csv'

if (!file.exists(source_file)) stop('Missing BoM cyclone database: ', source_file)
if (!file.exists(reef_file)) stop('Missing GBR reef features: ', reef_file)

# The first four lines are the BoM copyright and generation preamble.
tracks <- read_csv(
    source_file, skip = 4, trim_ws = TRUE, na = c('', ' '),
    show_col_types = FALSE
) |>
    mutate(
        NAME = trimws(NAME),
        track_time_utc = as.POSIXct(TM, tz = 'UTC'),
        event_year = if_else(
            month(track_time_utc) >= 11L,
            year(track_time_utc) + 1L, year(track_time_utc)
        )
    ) |>
    filter(
        TYPE == 'T',
        month(track_time_utc) %in% c(11L, 12L, 1:4),
        is.finite(LAT), is.finite(LON),
        LAT >= -35, LAT <= -5, LON >= 135, LON <= 165
    ) |>
    arrange(event_year, DISTURBANCE_ID, track_time_utc)

reefs <- st_read(reef_file, quiet = TRUE) |>
    st_drop_geometry() |>
    transmute(
        ReefID = trimws(LABEL_ID), ReefName = LOC_NAME_S,
        reef_longitude = as.numeric(X_COORD),
        reef_latitude = as.numeric(Y_COORD)
    ) |>
    filter(
        nzchar(ReefID), is.finite(reef_longitude), is.finite(reef_latitude),
        reef_longitude > 130, reef_latitude < 0
    ) |>
    distinct(ReefID, .keep_all = TRUE) |>
    arrange(ReefID)

reef_points <- st_as_sf(
    reefs, coords = c('reef_longitude', 'reef_latitude'),
    crs = 4326, remove = FALSE
)
event_years <- 1985:max(tracks$event_year, na.rm = TRUE)
exposure_parts <- vector('list', length(event_years))

for (year_index in seq_along(event_years)) {
    current_year <- event_years[[year_index]]
    event_tracks <- tracks |>
        filter(event_year == current_year)

    if (nrow(event_tracks) == 0L) {
        exposure_parts[[year_index]] <- reefs |>
            transmute(
                ReefID, ReefName, reef_longitude, reef_latitude,
                event_year = current_year,
                tc_min_distance_km = NA_real_,
                tc_nearest_name = NA_character_,
                tc_nearest_id = NA_character_,
                tc_nearest_time_utc = as.POSIXct(NA),
                tc_nearest_max_wind_ms = NA_real_,
                tc_nearest_pressure_hpa = NA_real_,
                tc_max_wind_within300km_ms = NA_real_,
                tc_min_pressure_within300km_hpa = NA_real_,
                tc_storms_within300km = 0L,
                tc_inside_reported_gale_radius = FALSE
            )
        next
    }

    track_points <- st_as_sf(
        event_tracks, coords = c('LON', 'LAT'), crs = 4326,
        remove = FALSE
    )
    distance_km <- drop_units(st_distance(reef_points, track_points)) / 1000
    nearest_index <- max.col(-distance_km, ties.method = 'first')
    nearest_distance <- distance_km[
        cbind(seq_len(nrow(distance_km)), nearest_index)
    ]

    wind_matrix <- matrix(
        event_tracks$MAX_WIND_SPD,
        nrow = nrow(distance_km), ncol = ncol(distance_km), byrow = TRUE
    )
    wind_matrix[distance_km > 300] <- NA_real_
    max_wind_300 <- rowMaxs(wind_matrix, na.rm = TRUE)
    max_wind_300[!is.finite(max_wind_300)] <- NA_real_

    pressure_matrix <- matrix(
        event_tracks$CENTRAL_PRES,
        nrow = nrow(distance_km), ncol = ncol(distance_km), byrow = TRUE
    )
    pressure_matrix[distance_km > 300] <- NA_real_
    min_pressure_300 <- rowMins(pressure_matrix, na.rm = TRUE)
    min_pressure_300[!is.finite(min_pressure_300)] <- NA_real_

    storm_ids <- unique(event_tracks$DISTURBANCE_ID)
    minimum_storm_distance <- matrix(
        NA_real_, nrow = nrow(distance_km), ncol = length(storm_ids)
    )
    for (storm_index in seq_along(storm_ids)) {
        fixes <- event_tracks$DISTURBANCE_ID == storm_ids[[storm_index]]
        minimum_storm_distance[, storm_index] <- rowMins(
            distance_km[, fixes, drop = FALSE]
        )
    }

    # BoM commonly supplies quadrant-specific gale radii even when the mean
    # radius is blank. Select the radius in the bearing from each storm fix to
    # each reef, falling back to the mean radius where necessary.
    reef_lon <- reefs$reef_longitude * pi / 180
    reef_lat <- reefs$reef_latitude * pi / 180
    track_lon <- event_tracks$LON * pi / 180
    track_lat <- event_tracks$LAT * pi / 180
    delta_lon <- outer(reef_lon, track_lon, '-')
    reef_lat_matrix <- matrix(
        reef_lat, nrow = nrow(distance_km), ncol = ncol(distance_km)
    )
    track_lat_matrix <- matrix(
        track_lat, nrow = nrow(distance_km), ncol = ncol(distance_km),
        byrow = TRUE
    )
    bearing <- atan2(
        sin(delta_lon) * cos(reef_lat_matrix),
        cos(track_lat_matrix) * sin(reef_lat_matrix) -
            sin(track_lat_matrix) * cos(reef_lat_matrix) * cos(delta_lon)
    ) * 180 / pi
    bearing <- (bearing + 360) %% 360

    quadrant_radius <- matrix(
        NA_real_, nrow = nrow(distance_km), ncol = ncol(distance_km)
    )
    radius_ne <- matrix(event_tracks$MN_RADIUS_GF_SECNE,
        nrow = nrow(distance_km), ncol = ncol(distance_km), byrow = TRUE)
    radius_se <- matrix(event_tracks$MN_RADIUS_GF_SECSE,
        nrow = nrow(distance_km), ncol = ncol(distance_km), byrow = TRUE)
    radius_sw <- matrix(event_tracks$MN_RADIUS_GF_SECSW,
        nrow = nrow(distance_km), ncol = ncol(distance_km), byrow = TRUE)
    radius_nw <- matrix(event_tracks$MN_RADIUS_GF_SECNW,
        nrow = nrow(distance_km), ncol = ncol(distance_km), byrow = TRUE)
    quadrant_radius[bearing < 90] <- radius_ne[bearing < 90]
    quadrant_radius[bearing >= 90 & bearing < 180] <-
        radius_se[bearing >= 90 & bearing < 180]
    quadrant_radius[bearing >= 180 & bearing < 270] <-
        radius_sw[bearing >= 180 & bearing < 270]
    quadrant_radius[bearing >= 270] <- radius_nw[bearing >= 270]
    mean_radius <- matrix(
        event_tracks$MN_RADIUS_GF_WIND,
        nrow = nrow(distance_km), ncol = ncol(distance_km), byrow = TRUE
    )
    quadrant_radius[!is.finite(quadrant_radius)] <-
        mean_radius[!is.finite(quadrant_radius)]
    inside_gale <- rowAnys(
        is.finite(quadrant_radius) & distance_km <= quadrant_radius
    )

    nearest <- event_tracks[nearest_index, , drop = FALSE]
    exposure_parts[[year_index]] <- reefs |>
        transmute(
            ReefID, ReefName, reef_longitude, reef_latitude,
            event_year = current_year,
            tc_min_distance_km = nearest_distance,
            tc_nearest_name = nearest$NAME,
            tc_nearest_id = nearest$DISTURBANCE_ID,
            tc_nearest_time_utc = nearest$track_time_utc,
            tc_nearest_max_wind_ms = nearest$MAX_WIND_SPD,
            tc_nearest_pressure_hpa = nearest$CENTRAL_PRES,
            tc_max_wind_within300km_ms = max_wind_300,
            tc_min_pressure_within300km_hpa = min_pressure_300,
            tc_storms_within300km = rowSums(
                minimum_storm_distance <= 300, na.rm = TRUE
            ),
            tc_inside_reported_gale_radius = inside_gale
        )
    cat('Completed BoM cyclone exposure year:', current_year, '\n')
}

exposure <- bind_rows(exposure_parts) |>
    arrange(event_year, ReefID)

track_export <- tracks |>
    select(
        NAME, DISTURBANCE_ID, track_time_utc, event_year, TYPE, CYC_TYPE,
        LAT, LON, CENTRAL_PRES, MAX_WIND_SPD,
        MN_RADIUS_GF_WIND, MN_RADIUS_GF_SECNE, MN_RADIUS_GF_SECSE,
        MN_RADIUS_GF_SECSW, MN_RADIUS_GF_SECNW, COMMENT
    )

metadata <- tribble(
    ~variable, ~units, ~interpretation,
    'tc_min_distance_km', 'km',
    'Minimum reef-to-recorded-fix distance during November-April',
    'tc_max_wind_within300km_ms', 'm/s',
    'Strongest best-track maximum wind at a fix within 300 km',
    'tc_min_pressure_within300km_hpa', 'hPa',
    'Lowest central pressure at a fix within 300 km',
    'tc_storms_within300km', 'count',
    'Distinct BoM tropical cyclones with at least one fix within 300 km',
    'tc_inside_reported_gale_radius', 'logical',
    'Any fix whose reported mean gale radius reaches the reef'
) |>
    mutate(
        source = 'BoM Australian Tropical Cyclone Database IDCKMSTM0S',
        caveat = paste(
            'Track proximity is not cyclone rainfall, flooding or wave damage.',
            'Older tracks and wind radii are less complete; fixes are not',
            'interpolated between observation times.'
        )
    )

write_csv(exposure, output_file)
write_csv(track_export, track_output)
write_csv(metadata, metadata_file)

cat('BoM tropical-cyclone fixes retained:', nrow(tracks), '\n')
cat('GBR reef features:', nrow(reefs), '\n')
cat('Reef-year exposure rows:', nrow(exposure), '\n')
