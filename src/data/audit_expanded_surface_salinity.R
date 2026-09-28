# Audit the expanded reef-level salinity delivery against its supplied GBR4 grid.
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(sf)
    library(stringr)
})

reef_dir <- 'data/SalinityData/reefs'
grid_dir <- 'data/SalinityData/gbr4'
output_dir <- 'output/surface_salinity_mortality'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

season_from_path <- function(path) {
    str_match(basename(path), '([0-9]{4}-[0-9]{4})[.]shp$')[, 2]
}
reef_files <- list.files(
    reef_dir,
    pattern = '^reefs_gbr4salinitylayers_[0-9]{4}-[0-9]{4}[.]shp$',
    full.names = TRUE
)
grid_files <- list.files(
    grid_dir,
    pattern = '^gbr4_freshwater_exposure30psu_hours_[0-9]{4}-[0-9]{4}[.]shp$',
    full.names = TRUE
)
pairs <- inner_join(
    tibble(reef_path = reef_files, season = season_from_path(reef_files)),
    tibble(grid_path = grid_files, season = season_from_path(grid_files)),
    by = 'season'
) |>
    arrange(season)
if (nrow(pairs) != length(reef_files) || nrow(pairs) != length(grid_files)) {
    stop('Reef and grid seasons do not pair one-to-one')
}

source_audit <- tibble()
grid_audit <- tibble()
reconstruction <- tibble()
reference_reef_ids <- reference_grid_ids <- reference_missing_ids <- NULL
reference_missing_grid_ids <- NULL
cached_grid_geometry <- NULL

for (i in seq_len(nrow(pairs))) {
    season <- pairs$season[[i]]
    message('Auditing ', season)
    reef <- st_read(pairs$reef_path[[i]], quiet = TRUE) |>
        mutate(ReefID = str_to_upper(str_trim(as.character(LABEL_ID))))
    grid <- st_read(pairs$grid_path[[i]], quiet = TRUE)
    reef_values <- st_drop_geometry(reef)
    grid_values <- st_drop_geometry(grid)

    reef_required <- c(
        'ReefID', 'FEAT_NAME', 'sss_min_me', 'exp30psu_m', 'exp26psu_m'
    )
    grid_required <- c(
        'SP_ID', 'sss_min', 'exposure_h', 'exp30psu_h', 'exp26psu_h'
    )
    if (length(setdiff(reef_required, names(reef_values))) ||
        length(setdiff(grid_required, names(grid_values)))) {
        stop('Missing required fields in ', season)
    }

    reef_ids <- sort(reef_values$ReefID)
    grid_ids <- sort(as.character(grid_values$SP_ID))
    missing_ids <- sort(
        reef_values$ReefID[!is.finite(reef_values$sss_min_me)]
    )
    missing_grid_ids <- sort(as.character(
        grid_values$SP_ID[!is.finite(grid_values$sss_min)]
    ))
    if (is.null(reference_reef_ids)) {
        reference_reef_ids <- reef_ids
        reference_grid_ids <- grid_ids
        reference_missing_ids <- missing_ids
        reference_missing_grid_ids <- missing_grid_ids
    }
    set_difference_count <- function(x, reference) {
        length(setdiff(union(x, reference), intersect(x, reference)))
    }

    source_audit <- bind_rows(source_audit, tibble(
        season,
        event_year = as.integer(substr(season, 6, 9)),
        reef_rows = nrow(reef_values),
        unique_reef_ids = n_distinct(reef_values$ReefID),
        duplicate_reef_ids = sum(duplicated(reef_values$ReefID)),
        nonreef_feature_rows = sum(
            reef_values$FEAT_NAME != 'Reef', na.rm = TRUE
        ),
        reef_id_set_diff_from_first = set_difference_count(
            reef_ids, reference_reef_ids
        ),
        missing_metric_rows = sum(!is.finite(reef_values$sss_min_me)),
        missing_id_set_diff_from_first = set_difference_count(
            missing_ids, reference_missing_ids
        ),
        negative_duration_rows = sum(
            reef_values$exp30psu_m < 0 | reef_values$exp26psu_m < 0,
            na.rm = TRUE
        ),
        hours26_exceeds_hours30_rows = sum(
            reef_values$exp26psu_m > reef_values$exp30psu_m + 1e-6,
            na.rm = TRUE
        ),
        positive_hours30_reefs = sum(
            reef_values$exp30psu_m > 0, na.rm = TRUE
        ),
        positive_hours26_reefs = sum(
            reef_values$exp26psu_m > 0, na.rm = TRUE
        ),
        minimum_sss_psu = min(reef_values$sss_min_me, na.rm = TRUE),
        maximum_hours30 = max(reef_values$exp30psu_m, na.rm = TRUE),
        maximum_hours26 = max(reef_values$exp26psu_m, na.rm = TRUE)
    ))
    grid_audit <- bind_rows(grid_audit, tibble(
        season,
        grid_rows = nrow(grid_values),
        unique_grid_ids = n_distinct(grid_values$SP_ID),
        duplicate_grid_ids = sum(duplicated(grid_values$SP_ID)),
        grid_id_set_diff_from_first = set_difference_count(
            grid_ids, reference_grid_ids
        ),
        missing_grid_id_set_diff_from_first = set_difference_count(
            missing_grid_ids, reference_missing_grid_ids
        ),
        exposure_alias_mismatch_rows = sum(
            abs(grid_values$exposure_h - grid_values$exp30psu_h) > 1e-6,
            na.rm = TRUE
        ),
        hours26_exceeds_hours30_rows = sum(
            grid_values$exp26psu_h > grid_values$exp30psu_h + 1e-6,
            na.rm = TRUE
        ),
        missing_sss_rows = sum(!is.finite(grid_values$sss_min)),
        missing_hours30_rows = sum(!is.finite(grid_values$exp30psu_h)),
        missing_hours26_rows = sum(!is.finite(grid_values$exp26psu_h)),
        minimum_sss_psu = min(grid_values$sss_min, na.rm = TRUE),
        maximum_hours30 = max(grid_values$exp30psu_h, na.rm = TRUE),
        maximum_hours26 = max(grid_values$exp26psu_h, na.rm = TRUE),
        bathymetry_is_grid_context_not_salinity_depth = TRUE
    ))

    target_rows <- unique(c(
        which.min(reef_values$sss_min_me),
        which.max(reef_values$exp30psu_m),
        which.max(reef_values$exp26psu_m),
        match('18-014', reef_values$ReefID),
        match('17-053', reef_values$ReefID)
    ))
    target_rows <- target_rows[is.finite(target_rows)]
    targets <- reef[target_rows, ] |>
        select(
            ReefID,
            supplied_sss = sss_min_me,
            supplied_hours30 = exp30psu_m,
            supplied_hours26 = exp26psu_m
        ) |>
        st_transform(3577)

    if (is.null(cached_grid_geometry)) {
        cached_grid_geometry <- st_geometry(st_transform(grid, 3577))
    } else if (!identical(grid_ids, reference_grid_ids)) {
        stop('Grid IDs changed; cached geometry cannot be reused in ', season)
    }
    grid_metric <- st_sf(
        grid_values |>
            select(SP_ID, sss_min, exp30psu_h, exp26psu_h),
        geometry = cached_grid_geometry,
        crs = 3577
    )
    hits <- suppressWarnings(st_intersection(grid_metric, targets)) |>
        mutate(overlap_area_m2 = as.numeric(st_area(geometry))) |>
        st_drop_geometry()
    reconstructed <- hits |>
        group_by(
            ReefID, supplied_sss, supplied_hours30, supplied_hours26
        ) |>
        summarise(
            reef_grid_intersections = n(),
            valid_grid_intersections = sum(is.finite(sss_min)),
            total_overlap_area_m2 = sum(overlap_area_m2),
            valid_overlap_area_m2 = sum(
                overlap_area_m2[is.finite(sss_min)]
            ),
            reconstructed_sss = weighted.mean(
                sss_min[is.finite(sss_min)],
                overlap_area_m2[is.finite(sss_min)]
            ),
            reconstructed_hours30 = weighted.mean(
                exp30psu_h[is.finite(exp30psu_h)],
                overlap_area_m2[is.finite(exp30psu_h)]
            ),
            reconstructed_hours26 = weighted.mean(
                exp26psu_h[is.finite(exp26psu_h)],
                overlap_area_m2[is.finite(exp26psu_h)]
            ),
            .groups = 'drop'
        ) |>
        mutate(
            season,
            valid_area_fraction =
                valid_overlap_area_m2 / total_overlap_area_m2,
            sss_difference = reconstructed_sss - supplied_sss,
            hours30_difference = reconstructed_hours30 - supplied_hours30,
            hours26_difference = reconstructed_hours26 - supplied_hours26
        )
    reconstruction <- bind_rows(reconstruction, reconstructed)
}

# The grid and missing-value mask are stable, so one exact spatial overlay gives
# a reusable valid-cell coverage fraction for every reef in all 14 seasons.
coverage_reef <- st_read(pairs$reef_path[[1]], quiet = TRUE) |>
    mutate(
        ReefID = str_to_upper(str_trim(as.character(LABEL_ID))),
        reef_area_m2 = as.numeric(st_area(st_transform(geometry, 3577)))
    ) |>
    select(ReefID, reef_area_m2) |>
    st_transform(3577)
coverage_grid <- st_read(pairs$grid_path[[1]], quiet = TRUE) |>
    select(SP_ID, sss_min) |>
    st_transform(3577)
coverage_hits <- suppressWarnings(st_intersection(
    coverage_grid, coverage_reef |> select(ReefID)
)) |>
    mutate(overlap_area_m2 = as.numeric(st_area(geometry))) |>
    st_drop_geometry()
coverage <- coverage_hits |>
    group_by(ReefID) |>
    summarise(
        grid_overlap_area_m2 = sum(overlap_area_m2),
        valid_grid_area_m2 = sum(overlap_area_m2[is.finite(sss_min)]),
        intersecting_grid_cells = n(),
        valid_grid_cells = sum(is.finite(sss_min)),
        .groups = 'drop'
    ) |>
    right_join(
        coverage_reef |> st_drop_geometry(),
        by = 'ReefID', relationship = 'one-to-one'
    ) |>
    mutate(
        across(
            c(grid_overlap_area_m2, valid_grid_area_m2,
              intersecting_grid_cells, valid_grid_cells),
            ~ coalesce(.x, 0)
        ),
        grid_overlap_fraction = pmin(grid_overlap_area_m2 / reef_area_m2, 1),
        valid_area_fraction = pmin(valid_grid_area_m2 / reef_area_m2, 1),
        coverage_class = case_when(
            valid_area_fraction >= 0.8 ~ 'high_ge_0.8',
            valid_area_fraction >= 0.5 ~ 'moderate_0.5_to_0.8',
            valid_area_fraction > 0 ~ 'low_below_0.5',
            TRUE ~ 'none'
        ),
        coverage_reference_season = pairs$season[[1]],
        coverage_static_across_seasons = TRUE
    ) |>
    arrange(ReefID)
write_csv(
    coverage,
    'data/processed/ereefs_surface_salinity_reef_grid_coverage.csv'
)
legacy_files <- list.files(
    'data/SalinityData',
    pattern = '^reefs_sssmin_exposure_mean_[0-9]{4}-[0-9]{4}[.]shp$',
    full.names = TRUE
)
legacy_audit <- bind_rows(lapply(legacy_files, function(path) {
    x <- st_read(path, quiet = TRUE) |>
        st_drop_geometry() |>
        mutate(ReefID = str_to_upper(str_trim(as.character(LABEL_ID))))
    duplicate_ids <- x |>
        count(ReefID) |>
        filter(n > 1) |>
        pull(ReefID)
    tibble(
        season = season_from_path(path),
        source_file = basename(path),
        rows = nrow(x),
        unique_reef_ids = n_distinct(x$ReefID),
        duplicate_reef_id_count = length(duplicate_ids),
        nonreef_rows = sum(x$FEAT_NAME != 'Reef', na.rm = TRUE),
        duplicated_ids_with_nonreef = n_distinct(
            x$ReefID[x$ReefID %in% duplicate_ids & x$FEAT_NAME != 'Reef']
        )
    )
}))

legacy_example <- tibble()
legacy_2024 <- legacy_files[grepl('2023-2024', legacy_files)]
new_2024 <- reef_files[grepl('2023-2024', reef_files)]
if (length(legacy_2024) == 1L && length(new_2024) == 1L) {
    legacy_example <- st_read(legacy_2024, quiet = TRUE) |>
        filter(str_to_upper(str_trim(as.character(LABEL_ID))) == '18-014') |>
        mutate(area_m2 = as.numeric(st_area(st_transform(geometry, 3577)))) |>
        st_drop_geometry() |>
        transmute(
            ReefID = LABEL_ID, feature_type = FEAT_NAME,
            feature_name = GBR_NAME, sss_min_psu = sss_min_me,
            hours_below30 = exposure_m, area_m2,
            delivery = 'superseded mixed-feature'
        ) |>
        bind_rows(
            st_read(new_2024, quiet = TRUE) |>
                filter(
                    str_to_upper(str_trim(as.character(LABEL_ID))) == '18-014'
                ) |>
                mutate(
                    area_m2 = as.numeric(st_area(st_transform(geometry, 3577)))
                ) |>
                st_drop_geometry() |>
                transmute(
                    ReefID = LABEL_ID, feature_type = FEAT_NAME,
                    feature_name = GBR_NAME, sss_min_psu = sss_min_me,
                    hours_below30 = exp30psu_m, area_m2,
                    delivery = 'expanded reef-only'
                )
        )
}

contract <- tibble(
    check = c(
        'paired_seasons', 'reef_ids_stable', 'grid_ids_stable',
        'grid_missing_mask_stable',
        'reef_ids_unique', 'reef_features_only', 'durations_nonnegative',
        'hours26_nested_within_hours30', 'grid_exposure_alias_consistent',
        'sampled_reef_means_reconstruct_from_grid',
        'temporal_window_dates_documented_in_delivery'
    ),
    passed = c(
        nrow(pairs) == 14L,
        all(source_audit$reef_id_set_diff_from_first == 0),
        all(grid_audit$grid_id_set_diff_from_first == 0),
        all(grid_audit$missing_grid_id_set_diff_from_first == 0),
        all(source_audit$duplicate_reef_ids == 0),
        all(source_audit$nonreef_feature_rows == 0),
        all(source_audit$negative_duration_rows == 0),
        all(source_audit$hours26_exceeds_hours30_rows == 0) &&
            all(grid_audit$hours26_exceeds_hours30_rows == 0),
        all(grid_audit$exposure_alias_mismatch_rows == 0),
        all(abs(reconstruction$sss_difference) < 1e-4) &&
            all(abs(reconstruction$hours30_difference) < 1e-4) &&
            all(abs(reconstruction$hours26_difference) < 1e-4),
        FALSE
    ),
    implication = c(
        'All 2010-11 through 2023-24 summers are present in both forms.',
        'The reef key is stable across seasons.',
        'The GBR4 cell key is stable across seasons.',
        'The GBR4 missing-cell mask is stable across seasons.',
        'Each expanded file has one record per ReefID.',
        'The expanded files exclude islands, rocks, cays and mainland.',
        'No invalid negative exposure durations were found.',
        '<26 PSU exposure is a valid nested subset of <30 PSU exposure.',
        'The legacy exposure_h field equals exp30psu_h.',
        'Selected extremes and known influential reefs reproduce exactly.',
        paste(
            'No start/end dates accompany the summary layers; absolute hours',
            'cannot be converted to seasonal fractions without provenance.'
        )
    )
)

write_csv(source_audit, file.path(output_dir, 'salinity_source_audit.csv'))
write_csv(grid_audit, file.path(output_dir, 'salinity_grid_audit.csv'))
write_csv(
    reconstruction,
    file.path(output_dir, 'salinity_grid_reconstruction_audit.csv')
)
write_csv(
    legacy_audit,
    file.path(output_dir, 'salinity_legacy_feature_mixing_audit.csv')
)
write_csv(
    legacy_example,
    file.path(output_dir, 'salinity_legacy_18_014_example.csv')
)
write_csv(contract, file.path(output_dir, 'salinity_source_contract.csv'))

required_failures <- contract |>
    filter(!passed, check != 'temporal_window_dates_documented_in_delivery')
if (nrow(required_failures) > 0L) {
    print(required_failures)
    stop('One or more required salinity source-contract checks failed')
}
cat('Audited', nrow(pairs), 'paired seasons.\n')
print(contract)
