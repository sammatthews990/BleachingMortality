# Build one-row-per-reef surface-salinity exposure metrics from the expanded
# reef-only eReefs GBR4 delivery. The source fields are already polygon-area
# means over valid 4 km cells; this script validates and harmonises them rather
# than re-aggregating unrelated features that happen to share a LABEL_ID.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(sf)
    library(stringr)
})

input_dir <- 'data/SalinityData/reefs'
output_file <- 'data/processed/ereefs_surface_salinity_reef_event.csv'
manifest_file <-
    'data/processed/ereefs_surface_salinity_reef_event_manifest.csv'

files <- list.files(
    input_dir,
    pattern = '^reefs_gbr4salinitylayers_[0-9]{4}-[0-9]{4}[.]shp$',
    full.names = TRUE
)
if (length(files) == 0L) {
    stop('No expanded reef salinity shapefiles found beneath ', input_dir)
}

parse_season <- function(path) {
    match <- str_match(basename(path), '([0-9]{4})-([0-9]{4})[.]shp$')
    if (any(is.na(match))) stop('Could not parse season from ', basename(path))
    list(
        label = paste0(match[, 2], '-', match[, 3]),
        start_year = as.integer(match[, 2]),
        end_year = as.integer(match[, 3])
    )
}

build_one_season <- function(path) {
    season_info <- parse_season(path)
    layer <- st_read(path, quiet = TRUE)
    required <- c(
        'LABEL_ID', 'FEAT_NAME', 'sss_min_me', 'exp30psu_m', 'exp26psu_m'
    )
    missing <- setdiff(required, names(layer))
    if (length(missing) > 0L) {
        stop(basename(path), ' lacks: ', paste(missing, collapse = ', '))
    }
    if (is.na(st_crs(layer))) stop(basename(path), ' has no CRS')

    layer <- layer |>
        mutate(
            ReefID = str_to_upper(str_trim(as.character(LABEL_ID))),
            source_feature_type = str_trim(as.character(FEAT_NAME)),
            reef_area_m2 = as.numeric(st_area(st_transform(geometry, 3577))),
            sss_min_area_mean_psu = as.numeric(sss_min_me),
            hours_below30_area_mean = as.numeric(exp30psu_m),
            hours_below26_area_mean = as.numeric(exp26psu_m)
        )

    if (anyDuplicated(layer[['ReefID']])) {
        stop(basename(path), ' has duplicate ReefID values')
    }
    if (any(layer[['source_feature_type']] != 'Reef', na.rm = TRUE)) {
        stop(basename(path), ' contains non-reef feature types')
    }
    if (any(is.na(layer[['ReefID']]) | layer[['ReefID']] == '')) {
        stop(basename(path), ' contains blank reef identifiers')
    }
    if (any(!is.finite(layer[['reef_area_m2']]) |
        layer[['reef_area_m2']] <= 0)) {
        stop(basename(path), ' contains invalid reef polygon areas')
    }
    if (any(
        is.finite(layer[['sss_min_area_mean_psu']]) &
            (layer[['sss_min_area_mean_psu']] < 0 |
                layer[['sss_min_area_mean_psu']] > 45)
    )) stop(basename(path), ' contains implausible minimum SSS')
    if (any(
        is.finite(layer[['hours_below30_area_mean']]) &
            layer[['hours_below30_area_mean']] < 0
    ) || any(
        is.finite(layer[['hours_below26_area_mean']]) &
            layer[['hours_below26_area_mean']] < 0
    )) stop(basename(path), ' contains negative exposure durations')
    if (any(
        layer[['hours_below26_area_mean']] >
            layer[['hours_below30_area_mean']] + 1e-6,
        na.rm = TRUE
    )) stop(basename(path), ' has hours below 26 exceeding hours below 30')

    layer |>
        st_drop_geometry() |>
        transmute(
            ReefID,
            event_year = season_info$end_year,
            season = season_info$label,
            season_start_year = season_info$start_year,
            season_end_year = season_info$end_year,
            source_reef_rows = 1L,
            reef_area_m2,
            # The delivery does not provide the valid-cell area denominator.
            # It is audited separately against the 4 km grid and must not be
            # inferred from the single reef-level record.
            valid_area_fraction = NA_real_,
            sss_min_area_mean_psu,
            hours_below30_area_mean,
            hours_below26_area_mean,
            salinity_deficit30_from_min_psu = pmax(
                30 - sss_min_area_mean_psu, 0
            ),
            salinity_deficit26_from_min_psu = pmax(
                26 - sss_min_area_mean_psu, 0
            ),
            log1p_hours_below30_area_mean = log1p(
                hours_below30_area_mean
            ),
            log1p_hours_below26_area_mean = log1p(
                hours_below26_area_mean
            ),
            source_feature_type,
            source_file = basename(path),
            source_metric_scope = 'two-dimensional surface salinity',
            source_aggregation = paste(
                'supplied reef-polygon area mean over valid GBR4 cells'
            ),
            event_year_rule = 'season-ending calendar year'
        )
}

metrics <- bind_rows(lapply(files, build_one_season)) |>
    arrange(event_year, ReefID)

if (anyDuplicated(metrics[c('ReefID', 'event_year')])) {
    stop('Surface salinity output has duplicate ReefID-event_year keys')
}

manifest <- metrics |>
    group_by(event_year, season, source_file) |>
    summarise(
        reefs = n(),
        reefs_with_metrics = sum(is.finite(sss_min_area_mean_psu)),
        reefs_missing_metrics = sum(!is.finite(sss_min_area_mean_psu)),
        reefs_with_hours_below30 = sum(
            hours_below30_area_mean > 0, na.rm = TRUE
        ),
        reefs_with_hours_below26 = sum(
            hours_below26_area_mean > 0, na.rm = TRUE
        ),
        minimum_sss_area_mean_psu = min(
            sss_min_area_mean_psu, na.rm = TRUE
        ),
        maximum_hours_below30_area_mean = max(
            hours_below30_area_mean, na.rm = TRUE
        ),
        maximum_hours_below26_area_mean = max(
            hours_below26_area_mean, na.rm = TRUE
        ),
        .groups = 'drop'
    )

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write_csv(metrics, output_file, na = '')
write_csv(manifest, manifest_file, na = '')

cat('Wrote', nrow(metrics), 'reef-season rows to', output_file, '\n')
print(manifest)
