# Build approximately annual coral-cover transitions from the raw AIMS survey
# series. Gains and losses are retained. This is a net-change dataset and does
# not apply the bleaching-only disturbance restriction used by the mortality
# analysis.

suppressPackageStartupMessages({
    library(dplyr)
    library(lubridate)
    library(readr)
})

workspace <- 'data/processed/01_exploratory_workspace.RData'
output_dir <- 'data/processed'
if (!file.exists(workspace)) stop('Missing exploratory workspace: ', workspace)
load(workspace)

decimal_year_to_date <- function(x) {
    x <- as.numeric(x)
    year <- floor(x)
    fraction <- x - year
    start <- ymd(paste0(year, '-01-01'))
    next_start <- ymd(paste0(year + 1L, '-01-01'))
    start + round(fraction * as.numeric(next_start - start))
}

collapse_source_text <- function(x) {
    x <- trimws(as.character(x))
    x <- unique(x[!is.na(x) & nzchar(x) & toupper(x) != 'NA'])
    if (length(x) == 0L) NA_character_ else paste(x, collapse = ' | ')
}

reef_reference <- read_csv(
    'data/AIMS-Reef_Reference.csv', show_col_types = FALSE
)
unique_reference <- reef_reference |>
    transmute(
        raw_reef_name_key = toupper(trimws(AIMS_REEF_NAME)),
        ReefID, ReefName, SECTOR, SECT_NAME
    ) |>
    add_count(raw_reef_name_key, name = 'reference_matches') |>
    filter(reference_matches == 1L) |>
    select(-reference_matches)

attach_reef_lineage <- function(data, raw_name, programme) {
    data <- data |>
        mutate(
            raw_reef_name_key = toupper(trimws(.data[[raw_name]])),
            programme_key = programme
        ) |>
        left_join(
            unique_reference, by = 'raw_reef_name_key',
            relationship = 'many-to-one'
        )

    # These two AIMS names are genuinely non-unique in the reference table.
    # The explicit rules reproduce their monitored sampling units without
    # reintroducing the cleaned-name collision that affected Mackay.
    data <- data |>
        mutate(
            ReefID = case_when(
                raw_reef_name_key == 'FRANKLAND ISLANDS' &
                    coalesce(reef_zone, '_') == 'West' ~ '17-013',
                raw_reef_name_key == 'FRANKLAND ISLANDS' ~ '17-012',
                raw_reef_name_key == 'SNAKE REEF' &
                    programme_key == 'manta' ~ '22-088a',
                raw_reef_name_key == 'SNAKE REEF' ~ '14-087',
                TRUE ~ ReefID
            )
        ) |>
        select(-ReefName, -SECTOR, -SECT_NAME) |>
        left_join(
            select(reef_reference, ReefID, ReefName, SECTOR, SECT_NAME) |>
                distinct(ReefID, .keep_all = TRUE),
            by = 'ReefID', relationship = 'many-to-one'
        )
    if (any(is.na(data$ReefID))) {
        stop('Annual survey rows remain unmatched to the reef reference')
    }
    data
}

manta <- df.AIMS.Mant |>
    mutate(reef_zone = '_') |>
    attach_reef_lineage('Reef_Name', 'manta') |>
    mutate(sampling_unit = paste(
        programme_key, ReefID, depth, reefpage_category, sep = '__'
    ))

benthic <- bind_rows(
    df.AIMS.Bent |> filter(project_code == 'LTMP') |>
        attach_reef_lineage('Reef_Name', 'ltmp'),
    df.AIMS.Bent |> filter(project_code == 'MMP') |>
        attach_reef_lineage('Reef_Name', 'mmp')
) |>
    mutate(sampling_unit = paste(
        programme_key, ReefID, depth, reef_zone,
        reefpage_category, sep = '__'
    ))

surveys <- bind_rows(manta, benthic) |>
    transmute(
        programme_key, source_observation_id = as.character(id),
        ReefID, ReefName, SECTOR, SECT_NAME, raw_reef_name_key,
        report_year = as.integer(report_year),
        survey_date = decimal_year_to_date(date),
        depth = as.numeric(depth), reef_zone, reefpage_category,
        sampling_unit, coral_cover = as.numeric(mean)
    ) |>
    arrange(sampling_unit, survey_date) |>
    group_by(sampling_unit) |>
    mutate(
        baseline_source_observation_id = lag(source_observation_id),
        baseline_report_year = lag(report_year),
        baseline_survey_date = lag(survey_date),
        pre_cover = lag(coral_cover),
        prior_change_pp = lag(100 * (coral_cover - lag(coral_cover)))
    ) |>
    ungroup() |>
    mutate(
        interval_days = as.integer(survey_date - baseline_survey_date),
        interval_years = interval_days / 365.25,
        post_cover = coral_cover,
        cover_change_pp = 100 * (post_cover - pre_cover),
        annualised_change_pp = cover_change_pp / interval_years,
        event_year = if_else(
            month(survey_date) >= 5L,
            year(survey_date), year(survey_date) - 1L
        ),
        available_space = 1 - pre_cover,
        cover_declined = cover_change_pp < 0,
        cover_increased = cover_change_pp > 0,
        region_block = case_when(
            SECTOR %in% c('CG', 'CL', 'PB') ~ 'Northern GBR',
            SECTOR %in% c('CA', 'CU', 'IN', 'TO', 'WH') ~ 'Central GBR',
            SECTOR %in% c('CB', 'PO', 'SW') ~ 'Southern GBR',
            TRUE ~ NA_character_
        )
    ) |>
    filter(
        is.finite(pre_cover), is.finite(post_cover),
        interval_days >= 240L, interval_days <= 550L,
        event_year >= 1985L, event_year <= 2025L
    ) |>
    mutate(transition_id = row_number())

# Pre-transition Acropora composition. Exact depth/year and reef/year values
# are preferred; a same-reef value up to two years old is the only fallback.
composition <- df.AIMS.full |>
    filter(
        data_type == 'photo-transect', domain_category == 'reef',
        purpose == 'COMPOSITION', variable == 'HARD CORAL'
    ) |>
    mutate(Reef_Name = domain_name) |>
    attach_reef_lineage('Reef_Name', 'ltmp') |>
    group_by(ReefID, report_year, depth) |>
    summarise(
        total_hard_coral = sum(mean, na.rm = TRUE),
        acropora_cover = sum(
            mean[grepl('Acropora', reefpage_category, ignore.case = TRUE)],
            na.rm = TRUE
        ),
        .groups = 'drop'
    ) |>
    mutate(
        prop_acropora = if_else(
            total_hard_coral > 0,
            pmin(pmax(acropora_cover / total_hard_coral, 0), 1),
            NA_real_
        )
    )
composition_reef <- composition |>
    group_by(ReefID, report_year) |>
    summarise(
        prop_acropora_reef = median(prop_acropora, na.rm = TRUE),
        .groups = 'drop'
    ) |>
    mutate(
        prop_acropora_reef = if_else(
            is.nan(prop_acropora_reef), NA_real_, prop_acropora_reef
        )
    )

transitions <- surveys |>
    left_join(
        select(
            composition, ReefID, report_year, depth,
            prop_acropora_depth = prop_acropora
        ),
        by = c(
            'ReefID', 'baseline_report_year' = 'report_year', 'depth'
        ),
        relationship = 'many-to-one'
    ) |>
    left_join(
        rename(composition_reef, prop_acropora_exact_reef = prop_acropora_reef),
        by = c('ReefID', 'baseline_report_year' = 'report_year'),
        relationship = 'many-to-one'
    ) |>
    mutate(
        prop_acropora_pre = coalesce(
            prop_acropora_depth, prop_acropora_exact_reef
        ),
        acropora_source = case_when(
            is.finite(prop_acropora_depth) ~ 'observed_depth_year',
            is.finite(prop_acropora_exact_reef) ~ 'observed_reef_year',
            TRUE ~ NA_character_
        )
    )

missing_composition <- transitions |>
    filter(!is.finite(prop_acropora_pre)) |>
    select(transition_id, ReefID, baseline_report_year) |>
    inner_join(composition_reef, by = 'ReefID', relationship = 'many-to-many') |>
    filter(
        report_year < baseline_report_year,
        report_year >= baseline_report_year - 2L,
        is.finite(prop_acropora_reef)
    ) |>
    arrange(transition_id, desc(report_year)) |>
    group_by(transition_id) |>
    slice_head(n = 1L) |>
    ungroup() |>
    select(
        transition_id,
        prop_acropora_past = prop_acropora_reef,
        acropora_latest_source_year = report_year
    )

transitions <- transitions |>
    left_join(missing_composition, by = 'transition_id') |>
    mutate(
        prop_acropora_pre = coalesce(
            prop_acropora_pre, prop_acropora_past
        ),
        acropora_source = case_when(
            !is.na(acropora_source) ~ acropora_source,
            is.finite(prop_acropora_past) ~ 'past_same_reef',
            TRUE ~ 'fold_median_required'
        ),
        acropora_latest_source_year = coalesce(
            acropora_latest_source_year,
            if_else(
                acropora_source %in% c(
                    'observed_depth_year', 'observed_reef_year'
                ),
                baseline_report_year, NA_integer_
            )
        )
    )

disturbance <- read_csv(
    'data/aims_ltmp/reef_disturbance.csv', show_col_types = FALSE
) |>
    mutate(
        raw_reef_name_key = toupper(trimws(aims_reef_name)),
        report_year = as.integer(year),
        survey_sample_type = toupper(sample_type)
    ) |>
    group_by(raw_reef_name_key, report_year, survey_sample_type) |>
    summarise(
        DISTURBANCE_TYPE = collapse_source_text(disturbance),
        storm_name = collapse_source_text(storm_name),
        disturbance_description = collapse_source_text(description),
        tooltip = collapse_source_text(tooltip),
        .groups = 'drop'
    )

transitions <- transitions |>
    mutate(
        survey_sample_type = if_else(
            programme_key == 'manta', 'MANTA', 'PPOINT'
        )
    ) |>
    left_join(
        disturbance,
        by = c('raw_reef_name_key', 'report_year', 'survey_sample_type'),
        relationship = 'many-to-one'
    ) |>
    mutate(
        disturbance_text = paste(
            coalesce(disturbance_description, ''),
            coalesce(storm_name, ''), coalesce(tooltip, '')
        ),
        disturbance_has_bleaching = grepl(
            'bleach|blch', disturbance_text, ignore.case = TRUE
        ),
        disturbance_has_cyclone = grepl(
            'cyclone|storm|jasper', disturbance_text, ignore.case = TRUE
        ),
        disturbance_has_flood = grepl(
            'flood|freshwater|low salinity',
            disturbance_text, ignore.case = TRUE
        ),
        disturbance_has_cots = grepl(
            'cots|crown-of-thorns', disturbance_text, ignore.case = TRUE
        )
    )

rrn <- read_csv(
    'data/processed/rrn_pressure_reef_year.csv', show_col_types = FALSE
)

# Carry the same thermal-history, optical and weather covariates used by the
# relative-mortality models into the annual transition table. The independently
# reconstructed Cheung-style fields and ERA5 extraction currently cover the
# focal bleaching years (2016, 2017, 2020, 2022 and 2024). Other years are
# retained and marked unavailable; imputation is performed inside each model
# training fold rather than here.
mean_or_na <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

mortality_environment <- read_csv(
    'data/processed/cheung_recreated_gbr_full.csv',
    show_col_types = FALSE, guess_max = Inf
) |>
    transmute(
        ReefID = LABEL_ID, event_year = as.integer(year),
        histmDHW6, yrsince6, histmDHW4, yrsince4,
        dhw10_load4, dhw_novelty10, ann_maxsst,
        winyear_mean, winyear_sd, mcur_90, dist_to_er_km,
        secc3m, secc3m_p10, cloudp_90
    ) |>
    group_by(ReefID, event_year) |>
    summarise(across(everything(), mean_or_na), .groups = 'drop')

weather <- read_csv(
    'data/processed/era5_weather_reef_year.csv',
    show_col_types = FALSE, guess_max = Inf
) |>
    transmute(
        ReefID = LABEL_ID, event_year = as.integer(year_x),
        log_coastal_rain30 = log1p(era5_coastal_rain_dec_mar_max_30day),
        era5_wind_mean = era5_reef_wind_q1_mean,
        era5_wind_calm_fraction = era5_reef_wind_q1_fraction_below_3,
        era5_coastal_distance_km = coastal_grid_distance_km
    ) |>
    group_by(ReefID, event_year) |>
    summarise(across(everything(), mean_or_na), .groups = 'drop')

transitions <- transitions |>
    left_join(
        rrn,
        by = c('ReefID' = 'LABEL_ID', 'event_year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        mortality_environment,
        by = c('ReefID', 'event_year'), relationship = 'many-to-one'
    ) |>
    left_join(
        weather,
        by = c('ReefID', 'event_year'), relationship = 'many-to-one'
    ) |>
    mutate(
        wqc_available = is.finite(wqc_freqcc12),
        cots_available = is.finite(cot_idwmeanpertow),
        mortality_environment_available = is.finite(histmDHW6),
        weather_available = is.finite(log_coastal_rain30),
        previous_change_available = is.finite(prior_change_pp),
        programme_factor = factor(
            programme_key, levels = c('ltmp', 'manta', 'mmp')
        ),
        reef_effect = factor(ReefID)
    ) |>
    arrange(event_year, ReefID, programme_key, survey_date)

if (any(
    transitions$acropora_latest_source_year >
        transitions$baseline_report_year,
    na.rm = TRUE
)) {
    stop('Annual Acropora reconstruction used a future observation')
}
if (any(transitions$post_cover < 0 | transitions$post_cover > 1)) {
    stop('Annual post-cover values fall outside 0--1')
}
if (any(
    transitions$raw_reef_name_key == 'MACKAY REEF' &
        transitions$ReefID != '16-015'
)) {
    stop('Singular Mackay Reef was assigned to the wrong reef lineage')
}
if (any(
    transitions$raw_reef_name_key == 'MACKAY REEFS' &
        transitions$ReefID != '15-024'
)) {
    stop('Plural Mackay Reefs was assigned to the wrong reef lineage')
}

manifest <- transitions |>
    group_by(programme_key) |>
    summarise(
        transitions = n(), reefs = n_distinct(ReefID),
        first_event_year = min(event_year),
        last_event_year = max(event_year),
        positive_changes = sum(cover_change_pp > 0),
        negative_changes = sum(cover_change_pp < 0),
        zero_changes = sum(cover_change_pp == 0),
        median_interval_days = median(interval_days),
        acropora_observed = sum(
            acropora_source %in% c(
                'observed_depth_year', 'observed_reef_year'
            )
        ),
        .groups = 'drop'
    )

write_csv(
    transitions,
    file.path(output_dir, 'annual_coral_transitions.csv'), na = ''
)
saveRDS(transitions, file.path(output_dir, 'annual_coral_transitions.rds'))
write_csv(manifest, file.path(output_dir, 'annual_coral_transitions_manifest.csv'))
print(manifest)
