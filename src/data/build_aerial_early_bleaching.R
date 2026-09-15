# Harmonise historical and contemporary GBR aerial bleaching surveys for the
# explicitly timed within-event update. Raw source files are never modified.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(sf)
    library(stringr)
})

root <- normalizePath('.', winslash = '/', mustWork = TRUE)
source_dir <- file.path(root, 'data', 'AerialBleach2022_25')
hughes_dir <- file.path(root, 'data', 'Hughes2021', 'data')
processed_dir <- file.path(root, 'data', 'processed')
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)

crosswalk_file <- file.path(root, 'config', 'aerial_bleaching_score_crosswalk.csv')
score_crosswalk <- read_csv(crosswalk_file, show_col_types = FALSE)
required_crosswalk <- c(
    'source_scale', 'native_category', 'midpoint_bleached_cover',
    'high_bleaching', 'review_action'
)
if (!all(required_crosswalk %in% names(score_crosswalk))) {
    stop('Aerial score crosswalk is missing required columns')
}
if (anyDuplicated(score_crosswalk[c('source_scale', 'native_category')])) {
    stop('Aerial score crosswalk has duplicate scale/category keys')
}

score_midpoint <- function(score, source_scale) {
    map <- score_crosswalk |>
        filter(.data$source_scale == .env$source_scale) |>
        arrange(native_category)
    if (source_scale == 'contemporary_aims_sop11') {
        return(approx(
            x = map$native_category, y = map$midpoint_bleached_cover,
            xout = score, method = 'linear', rule = 1
        )$y)
    }
    map$midpoint_bleached_cover[match(score, map$native_category)]
}

normalise_reef_id <- function(x) str_to_upper(str_trim(as.character(x)))
mean_or_na <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
min_or_na <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
date_min_or_na <- function(x) if (all(is.na(x))) as.Date(NA) else min(x, na.rm = TRUE)
date_max_or_na <- function(x) if (all(is.na(x))) as.Date(NA) else max(x, na.rm = TRUE)

parse_dmy <- function(x) as.Date(x, format = '%d/%m/%Y')

require_fields <- function(x, fields, source_name) {
    missing <- setdiff(fields, names(x))
    if (length(missing)) {
        stop(source_name, ' is missing fields: ', paste(missing, collapse = ', '))
    }
}

modern_status <- function(score, context) {
    context <- str_to_lower(coalesce(as.character(context), ''))
    case_when(
        is.finite(score) & between(score, 0, 5) ~ 'scored',
        str_detect(context, 'sediment') ~ 'not_scored_no_view_sediment',
        str_detect(context, 'deep') ~ 'not_scored_no_view_depth',
        str_detect(context, 'low|no live') ~ 'not_scored_low_or_no_live_coral_cover',
        str_detect(context, 'no.?view') ~ 'not_scored_no_view',
        TRUE ~ 'not_scored_unknown_status'
    )
}

source_files <- c(
    `2022` = file.path(source_dir, '2022_AIMSGBRMPA_Cantin_GBRAerial_Final.csv'),
    `2024` = file.path(source_dir, '2024_AIMSGBRMPA_Cantin_GBRAerial_Final.csv'),
    `2025` = file.path(source_dir, '2025_AIMSGBRMPA_Cantin_GBRAerial_Final.csv')
)
missing_files <- source_files[!file.exists(source_files)]
if (length(missing_files)) {
    stop('Missing contemporary aerial files: ', paste(missing_files, collapse = ', '))
}

reef_coordinates <- st_read(
    file.path(root, 'data', 'Great_Barrier_Reef_Features',
              'Great_Barrier_Reef_Features.shp'), quiet = TRUE
) |>
    st_drop_geometry() |>
    transmute(
        ReefID = normalise_reef_id(LABEL_ID),
        lon = as.numeric(X_COORD), lat = as.numeric(Y_COORD)
    ) |>
    filter(between(lon, 142, 155), between(lat, -26, -9)) |>
    group_by(ReefID) |>
    summarise(lon = mean(lon), lat = mean(lat), .groups = 'drop')

# The Hughes source retains its native category and binary response. Survey
# dates are absent, and no continuous modern severity is manufactured from it.
hughes_observations <- bind_rows(lapply(c(2016L, 2017L, 2020L), function(year) {
    source_file <- file.path(hughes_dir, paste0('bleach_', year, '.csv'))
    raw <- read_csv(source_file, show_col_types = FALSE)
    require_fields(raw, c('ReefID', 'score', 'bin.score'), basename(source_file))
    raw |>
        transmute(
            native_unit_id = paste0(year, ':', normalise_reef_id(ReefID)),
            ReefID = normalise_reef_id(ReefID),
            ReefName = NA_character_, event_year = year,
            survey_date = as.Date(NA),
            lon = NA_real_, lat = NA_real_,
            native_score_mean = as.numeric(score),
            native_score_rounded_up = as.numeric(score),
            native_score_min = as.numeric(score),
            native_score_max = as.numeric(score),
            native_score_sd = NA_real_, native_observation_count = NA_real_,
            native_status = NA_character_,
            observation_status = 'scored_date_unavailable',
            observation_usable = TRUE,
            source_score_scale = 'historical_hughes',
            aerial_bleached_cover_midpoint = score_midpoint(
                as.numeric(score), 'historical_hughes'
            ),
            aerial_severity_01 = aerial_bleached_cover_midpoint,
            aerial_high_bleaching = as.numeric(bin.score),
            score_review_flag = as.numeric(score) == 5,
            source_family = 'Hughes historical aerial survey',
            source_file = file.path('data', 'Hughes2021', 'data',
                                    basename(source_file)),
            source_provider = 'ARC Centre of Excellence for Coral Reef Studies',
            source_method = 'Historical reef-level aerial bleaching category',
            score_definition = paste(
                'Native category retained; bin.score is the supplied indicator',
                'for category 3 or higher. Survey date is unavailable.'
            ),
            score_harmonisation = paste(
                'Documented cover-bin midpoint; score 5 top-coded to the',
                'historical >60% class'
            ),
            timing_provenance = 'Survey date absent in supplied Hughes source',
            dhw_native = as.numeric(if ('DHW_5km_max' %in% names(raw))
                raw$DHW_5km_max else NA_real_)
        ) |>
        left_join(reef_coordinates, by = 'ReefID', suffix = c('', '_reference')) |>
        mutate(
            lon = coalesce(lon, lon_reference),
            lat = coalesce(lat, lat_reference)
        ) |>
        select(-lon_reference, -lat_reference)
}))

raw_2022 <- read_csv(source_files[['2022']], show_col_types = FALSE)
require_fields(
    raw_2022,
    c('LOC_NAME_S', 'LABEL_ID', 'Y_COORD', 'X_COORD', 'AIR_BLEACH_CAT',
      'year', 'surveyDate', 'AerialBLC_ReefAvg',
      'AerialBLC_ReefAvgRoundedUP'),
    basename(source_files[['2022']])
)
observations_2022 <- raw_2022 |>
    transmute(
        native_unit_id = paste(
            '2022', LOC_NAME_S, LABEL_ID, X_COORD, Y_COORD, surveyDate,
            sep = ':'
        ),
        ReefID = normalise_reef_id(LABEL_ID), ReefName = LOC_NAME_S,
        event_year = as.integer(year), survey_date = parse_dmy(surveyDate),
        lon = as.numeric(X_COORD), lat = as.numeric(Y_COORD),
        native_score_mean = as.numeric(AerialBLC_ReefAvg),
        native_score_rounded_up = as.numeric(AerialBLC_ReefAvgRoundedUP),
        native_score_min = NA_real_, native_score_max = NA_real_,
        native_score_sd = NA_real_, native_observation_count = NA_real_,
        native_status = as.character(AIR_BLEACH_CAT),
        dhw_native = as.numeric(SAMPLE_DHW1)
    ) |>
    distinct(native_unit_id, .keep_all = TRUE) |>
    mutate(
        observation_status = modern_status(
            native_score_rounded_up, native_status
        ),
        observation_usable = observation_status == 'scored',
        source_score_scale = 'contemporary_aims_sop11',
        aerial_bleached_cover_midpoint = if_else(
            observation_usable,
            score_midpoint(native_score_mean, 'contemporary_aims_sop11'), NA_real_
        ),
        aerial_severity_01 = aerial_bleached_cover_midpoint,
        score_review_flag = FALSE,
        aerial_high_bleaching = if_else(
            observation_usable,
            as.numeric(native_score_rounded_up >= 3), NA_real_
        ),
        source_family = 'AIMS-GBRMPA contemporary aerial survey',
        source_file = file.path('data', 'AerialBleach2022_25',
                                basename(source_files[['2022']])),
        source_provider = paste(
            'Australian Institute of Marine Science and Great Barrier Reef',
            'Marine Park Authority'
        ),
        source_method = 'SOP 11 aerial survey of shallow visible coral cover',
        score_definition = paste(
            'Valid categories: 0 <1%, 1 1-10%, 2 11-30%, 3 31-60%,',
            '4 61-90%, 5 >90% bleached. Codes above 5 are status codes.'
        ),
        score_harmonisation = paste(
            'Piecewise interpolation between SOP 11 cover-bin midpoints'
        ),
        score_harmonisation = 'SOP 11 cover-bin midpoint for supplied rounded category',
        timing_provenance = 'Survey date supplied in source; interpreted as AEST'
    )

raw_2024 <- read_csv(source_files[['2024']], show_col_types = FALSE)
require_fields(
    raw_2024,
    c('Reef_Count', 'LOC_NAME_S', 'LABEL_ID', 'Lat', 'Long',
      'AerialCat_RoundUP', 'AerialBleaching_Category_mean',
      'Aerial_OBS_Count', 'AerialOBS_Min', 'AerialOBS_Max',
      'AerialOBS_StDev', 'Date_AEST_min', 'NoBleachingCATEGORY'),
    basename(source_files[['2024']])
)
observations_2024 <- raw_2024 |>
    transmute(
        native_unit_id = paste0('2024:', Reef_Count),
        ReefID = normalise_reef_id(LABEL_ID), ReefName = LOC_NAME_S,
        event_year = 2024L, survey_date = parse_dmy(Date_AEST_min),
        lon = as.numeric(Long), lat = as.numeric(Lat),
        native_score_mean = as.numeric(AerialBleaching_Category_mean),
        native_score_rounded_up = as.numeric(AerialCat_RoundUP),
        native_score_min = as.numeric(AerialOBS_Min),
        native_score_max = as.numeric(AerialOBS_Max),
        native_score_sd = as.numeric(AerialOBS_StDev),
        native_observation_count = as.numeric(Aerial_OBS_Count),
        native_status = as.character(NoBleachingCATEGORY),
        dhw_native = as.numeric(SAMPLE_DHW1)
    ) |>
    distinct(native_unit_id, .keep_all = TRUE) |>
    mutate(
        observation_status = modern_status(
            native_score_rounded_up, native_status
        ),
        observation_usable = observation_status == 'scored',
        source_score_scale = 'contemporary_aims_sop11',
        aerial_bleached_cover_midpoint = if_else(
            observation_usable,
            score_midpoint(native_score_mean, 'contemporary_aims_sop11'), NA_real_
        ),
        aerial_severity_01 = aerial_bleached_cover_midpoint,
        score_review_flag = FALSE,
        aerial_high_bleaching = if_else(
            observation_usable,
            as.numeric(native_score_rounded_up >= 3), NA_real_
        ),
        source_family = 'AIMS-GBRMPA contemporary aerial survey',
        source_file = file.path('data', 'AerialBleach2022_25',
                                basename(source_files[['2024']])),
        source_provider = paste(
            'Australian Institute of Marine Science and Great Barrier Reef',
            'Marine Park Authority'
        ),
        source_method = 'SOP 11 aerial survey of shallow visible coral cover',
        score_definition = paste(
            'Valid categories: 0 <1%, 1 1-10%, 2 11-30%, 3 31-60%,',
            '4 61-90%, 5 >90% bleached. Codes 6-8 are status codes.'
        ),
        score_harmonisation = paste(
            'Piecewise interpolation between SOP 11 cover-bin midpoints'
        ),
        timing_provenance = 'Minimum survey date and time supplied in AEST'
    )

raw_2025 <- read_csv(source_files[['2025']], show_col_types = FALSE)
require_fields(
    raw_2025,
    c('fid', 'LOC_NAME_S', 'LABEL_ID', 'X_COORD', 'Y_COORD',
      'survey_date', 'AerialBLC_CAT_RoundUP'),
    basename(source_files[['2025']])
)
observations_2025 <- raw_2025 |>
    mutate(source_row = row_number()) |>
    transmute(
        native_unit_id = paste0('2025:', fid, ':', source_row),
        ReefID = normalise_reef_id(LABEL_ID), ReefName = LOC_NAME_S,
        event_year = 2025L, survey_date = parse_dmy(survey_date),
        lon = as.numeric(X_COORD), lat = as.numeric(Y_COORD),
        native_score_mean = as.numeric(AerialBLC_CAT_RoundUP),
        native_score_rounded_up = as.numeric(AerialBLC_CAT_RoundUP),
        native_score_min = as.numeric(AerialBLC_CAT_RoundUP),
        native_score_max = as.numeric(AerialBLC_CAT_RoundUP),
        native_score_sd = NA_real_, native_observation_count = 1,
        native_status = NA_character_, dhw_native = as.numeric(DHW)
    ) |>
    mutate(
        observation_status = modern_status(
            native_score_rounded_up, native_status
        ),
        observation_usable = observation_status == 'scored',
        source_score_scale = 'contemporary_aims_sop11',
        aerial_bleached_cover_midpoint = if_else(
            observation_usable,
            score_midpoint(native_score_mean, 'contemporary_aims_sop11'), NA_real_
        ),
        aerial_severity_01 = aerial_bleached_cover_midpoint,
        score_review_flag = FALSE,
        aerial_high_bleaching = if_else(
            observation_usable,
            as.numeric(native_score_rounded_up >= 3), NA_real_
        ),
        source_family = 'AIMS-GBRMPA contemporary aerial survey',
        source_file = file.path('data', 'AerialBleach2022_25',
                                basename(source_files[['2025']])),
        source_provider = paste(
            'Australian Institute of Marine Science and Great Barrier Reef',
            'Marine Park Authority'
        ),
        source_method = 'SOP 11 aerial survey of shallow visible coral cover',
        score_definition = paste(
            'Rounded category supplied: 0 <1%, 1 1-10%, 2 11-30%,',
            '3 31-60%, 4 61-90%, 5 >90% bleached.'
        ),
        score_harmonisation = paste(
            'Piecewise interpolation between SOP 11 cover-bin midpoints'
        ),
        score_harmonisation = 'SOP 11 cover-bin midpoint for supplied rounded category',
        timing_provenance = 'Survey date supplied in source; interpreted as AEST'
    )

observations <- bind_rows(
    hughes_observations, observations_2022, observations_2024, observations_2025
) |>
    mutate(
        available_by_march_31 = case_when(
            is.na(survey_date) ~ NA,
            TRUE ~ survey_date <= as.Date(paste0(event_year, '-03-31'))
        ),
        available_by_april_30 = case_when(
            is.na(survey_date) ~ NA,
            TRUE ~ survey_date <= as.Date(paste0(event_year, '-04-30'))
        ),
        product_mode = 'within_event_update',
        use_in_initial_forecast = FALSE,
        use_in_selected_model = FALSE
    ) |>
    arrange(event_year, ReefID, native_unit_id)

if (anyDuplicated(observations[c('event_year', 'native_unit_id')])) {
    stop('Duplicate aerial native-unit keys remain after source-specific cleaning')
}
if (any(observations$observation_usable &
        !between(observations$native_score_rounded_up, 0, 5))) {
    stop('A status code above 5 was incorrectly retained as bleaching severity')
}
if (any(observations$observation_usable &
        !between(observations$aerial_high_bleaching, 0, 1))) {
    stop('Harmonised high-bleaching indicator is outside 0-1')
}
if (any(observations$observation_usable &
        !between(observations$aerial_bleached_cover_midpoint, 0, 1))) {
    stop('Harmonised aerial cover midpoint is outside 0-1')
}

summarise_cutoff <- function(rows, cutoff_field, suffix) {
    rows |>
        filter(observation_usable, .data[[cutoff_field]] %in% TRUE) |>
        group_by(ReefID, event_year) |>
        summarise(
            severity = mean_or_na(aerial_severity_01),
            high = max_or_na(aerial_high_bleaching),
            units = n(), .groups = 'drop'
        ) |>
        rename_with(~ paste0(.x, '_', suffix), c(severity, high, units))
}

reef_event <- observations |>
    group_by(ReefID, event_year) |>
    summarise(
        ReefName = first(na.omit(ReefName), default = NA_character_),
        lon = mean_or_na(lon), lat = mean_or_na(lat),
        survey_date = date_min_or_na(survey_date),
        first_survey_date = date_min_or_na(survey_date),
        last_survey_date = date_max_or_na(survey_date),
        native_score_mean = mean_or_na(native_score_mean[observation_usable]),
        native_score_rounded_up = if_else(
            any(observation_usable),
            ceiling(mean(native_score_mean[observation_usable], na.rm = TRUE)),
            NA_real_
        ),
        native_score_min = min_or_na(native_score_min[observation_usable]),
        native_score_max = max_or_na(native_score_max[observation_usable]),
        aerial_bleached_cover_midpoint = mean_or_na(
            aerial_bleached_cover_midpoint
        ),
        aerial_severity_01 = mean_or_na(aerial_severity_01),
        aerial_high_bleaching = max_or_na(aerial_high_bleaching),
        observation_usable = any(observation_usable),
        n_native_units = n(),
        n_scored_units = sum(observation_usable),
        n_native_observations = sum(native_observation_count, na.rm = TRUE),
        score_review_flag = any(score_review_flag),
        source_score_scale = paste(
            sort(unique(source_score_scale)), collapse = ' | '
        ),
        identifier_collision = n_distinct(ReefName) > 1 |
            n_distinct(round(lon, 4), round(lat, 4)) > 1,
        source_family = paste(sort(unique(source_family)), collapse = ' | '),
        source_file = paste(sort(unique(source_file)), collapse = ' | '),
        source_provider = paste(sort(unique(source_provider)), collapse = ' | '),
        source_method = paste(sort(unique(source_method)), collapse = ' | '),
        score_definition = paste(sort(unique(score_definition)), collapse = ' | '),
        score_harmonisation = paste(
            sort(unique(score_harmonisation)), collapse = ' | '
        ),
        timing_provenance = paste(
            sort(unique(timing_provenance)), collapse = ' | '
        ),
        product_mode = 'within_event_update',
        use_in_initial_forecast = FALSE,
        use_in_selected_model = FALSE,
        .groups = 'drop'
    ) |>
    left_join(
        summarise_cutoff(observations, 'available_by_march_31', 'march_31'),
        by = c('ReefID', 'event_year'), relationship = 'one-to-one'
    ) |>
    left_join(
        summarise_cutoff(observations, 'available_by_april_30', 'april_30'),
        by = c('ReefID', 'event_year'), relationship = 'one-to-one'
    ) |>
    mutate(
        available_by_march_31 = units_march_31 > 0,
        available_by_april_30 = units_april_30 > 0
    ) |>
    arrange(event_year, ReefID)

if (anyDuplicated(reef_event[c('ReefID', 'event_year')])) {
    stop('Harmonised aerial table has duplicate reef-event keys')
}

conflict_audit <- observations |>
    group_by(ReefID, event_year) |>
    summarise(
        n_native_units = n(), n_names = n_distinct(ReefName),
        coordinate_span_lon = max_or_na(lon) - min_or_na(lon),
        coordinate_span_lat = max_or_na(lat) - min_or_na(lat),
        score_range = max_or_na(native_score_mean[observation_usable]) -
            min_or_na(native_score_mean[observation_usable]),
        .groups = 'drop'
    ) |>
    filter(n_native_units > 1 | n_names > 1 | score_range > 0)

file_audit <- bind_rows(
    tibble(
        event_year = c(2016L, 2017L, 2020L),
        source_file = file.path(
            'data', 'Hughes2021', 'data',
            paste0('bleach_', c(2016L, 2017L, 2020L), '.csv')
        ),
        raw_rows = sapply(c(2016L, 2017L, 2020L), function(year) nrow(
            read_csv(file.path(hughes_dir, paste0('bleach_', year, '.csv')),
                     show_col_types = FALSE)
        ))
    ),
    tibble(
        event_year = c(2022L, 2024L, 2025L),
        source_file = file.path(
            'data', 'AerialBleach2022_25', basename(source_files)
        ),
        raw_rows = c(nrow(raw_2022), nrow(raw_2024), nrow(raw_2025))
    )
) |>
    left_join(
        observations |>
            group_by(event_year, source_file) |>
            summarise(
                native_units = n(), unique_reef_ids = n_distinct(ReefID),
                scored_units = sum(observation_usable),
                excluded_status_units = sum(!observation_usable),
                score_review_units = sum(score_review_flag),
                missing_survey_dates = sum(is.na(survey_date)),
                first_survey_date = date_min_or_na(survey_date),
                last_survey_date = date_max_or_na(survey_date),
                march_31_scored_units = sum(
                    observation_usable & available_by_march_31 %in% TRUE
                ),
                april_30_scored_units = sum(
                    observation_usable & available_by_april_30 %in% TRUE
                ),
                .groups = 'drop'
            ),
        by = c('event_year', 'source_file'), relationship = 'one-to-one'
    ) |>
    mutate(
        duplicate_gis_join_rows_removed = raw_rows - native_units,
        file_md5 = as.character(tools::md5sum(file.path(root, source_file))),
        provider = if_else(
            event_year >= 2022,
            'Australian Institute of Marine Science and GBRMPA',
            'ARC Centre of Excellence / Hughes historical source'
        ),
        method = if_else(
            event_year >= 2022,
            'AIMS SOP 11 aerial bleaching survey',
            'Historical aerial bleaching survey; date unavailable'
        ),
        operational_role = if_else(
            event_year == 2025L,
            'within-event update validation with post-May 2025 mortality outcomes',
            'within-event update validation only'
        )
    )

provenance <- tibble(
    contract_version = 2L,
    product_mode = 'within_event_update',
    initial_forecast_allowed = FALSE,
    selected_model_changed = FALSE,
    category_definition = paste(
        'Historical: 0 <1%, 1 1-10%, 2 10-30%, 3 30-60%, 4 >60%;',
        'contemporary: 0 <1%, 1 1-10%, 2 11-30%, 3 31-60%,',
        '4 61-90%, 5 >90% of visible shallow living coral cover bleached'
    ),
    common_severity_definition = paste(
        'Cover-bin midpoint on a 0-1 scale. Contemporary mean categories use',
        'piecewise interpolation; unexpected historical 2020 score 5 is',
        'flagged and top-coded to the documented historical >60% class.'
    ),
    high_bleaching_definition = paste(
        'Historical-compatible indicator: native or rounded category >=3;',
        'retained as a binary sensitivity beside the common midpoint severity'
    ),
    excluded_status_rule = paste(
        'Codes above 5 and corresponding no-view, sediment, depth or',
        'low/no-live-coral statuses are not bleaching severity'
    ),
    crosswalk_file = 'config/aerial_bleaching_score_crosswalk.csv',
    crosswalk_md5 = as.character(tools::md5sum(crosswalk_file)),
    historical_reference = paste(
        'Hughes et al. methods; dataset doi:10.17632/tncdys47mh.1'
    ),
    method_reference = 'AIMS SOP 11 version 3 (2022), doi:10.25845/n00q-z603',
    source_url = paste0(
        'https://www.aims.gov.au/sites/default/files/2022-06/',
        'AIMS_SOP11v3_Aerial-Surveys-Coral-Bleaching_202206.pdf'
    )
)

write_csv(
    observations,
    file.path(processed_dir, 'aerial_early_bleaching_observation.csv'), na = ''
)
write_csv(
    reef_event,
    file.path(processed_dir, 'aerial_early_bleaching_reef_event.csv'), na = ''
)
write_csv(
    file_audit,
    file.path(processed_dir, 'aerial_early_bleaching_source_audit.csv'), na = ''
)
write_csv(
    conflict_audit,
    file.path(processed_dir, 'aerial_early_bleaching_conflict_audit.csv'), na = ''
)
write_csv(
    provenance,
    file.path(processed_dir, 'aerial_early_bleaching_provenance.csv'), na = ''
)

message('Aerial observations: ', nrow(observations))
message('Harmonised aerial reef-events: ', nrow(reef_event))
message('Scored 2022 reef-events: ', sum(
    reef_event$event_year == 2022L & reef_event$observation_usable
))
message('Scored 2024 reef-events: ', sum(
    reef_event$event_year == 2024L & reef_event$observation_usable
))
message('Scored 2025 reef-events: ', sum(
    reef_event$event_year == 2025L & reef_event$observation_usable
))
