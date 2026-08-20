# Build the fixed outer resampling design used by the mortality models.
#
# This script deliberately keeps the workflow linear and inspectable. It joins
# the canonical outcomes to the reconstructed environmental table, assigns
# three complementary blocked validation schemes, checks them, and writes the
# row-level assignments. The selected production model is later refitted using
# every event, including 2024.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

set.seed(202408L)

outcome_files <- c(
    manta = "data/processed/mortality_outcomes_manta.rds",
    ltmp = "data/processed/mortality_outcomes_ltmp.rds",
    mmp = "data/processed/mortality_outcomes_mmp.rds"
)
environment_file <- "data/processed/cheung_recreated_gbr_full.csv"
output_dir <- "data/processed"

acropora_workspace_file <- 'data/processed/01_exploratory_workspace.RData'
missing_inputs <- c(outcome_files, environment_file, acropora_workspace_file)[
    !file.exists(c(outcome_files, environment_file, acropora_workspace_file))
]
if (length(missing_inputs) > 0L) {
    stop("Missing validation inputs: ", paste(missing_inputs, collapse = ", "))
}

# Primary modelling is restricted to years covered by the independently
# reconstructed environmental table. Historical outcomes remain available in
# the outcome files for simpler DHW-only sensitivity analyses.
model_years <- c(2016L, 2017L, 2020L, 2022L, 2024L)

predictor_columns <- c(
    "ann_maxdhw", "histmDHW6", "yrsince6", "histmDHW4", "yrsince4",
    "ann_maxsst", "winyear_mean", "winyear_sd", "mcur_90",
    "dist_to_er_km", "secc3m", "cloudp_90"
)
predictor_columns <- c(
    predictor_columns,
    'dhw10_load4', 'dhw_novelty10', 'secc3m_p10'
)

# Read and stack outcomes. Programme identity is retained and the tables are
# separated again before writing.
outcomes <- bind_rows(lapply(names(outcome_files), function(name) {
    readRDS(outcome_files[[name]]) |>
        mutate(programme_key = name)
})) |>
    mutate(
        reef_name_key = toupper(trimws(as.character(ReefName))),
        reef_id_key = toupper(trimws(as.character(ReefID))),
        event_year = as.integer(event_year),
        legacy_survey_max_dhw = MaxDHW.mean
    )

if (anyDuplicated(paste(outcomes$programme_key, outcomes$source_observation_id))) {
    stop("Outcome observation IDs are not unique within programme")
}

environment <- read_csv(
    environment_file,
    show_col_types = FALSE,
    guess_max = Inf
) |>
    mutate(
        event_year = as.integer(year),
        reef_name_key = toupper(trimws(as.character(LOC_NAME_S))),
        reef_id_key = toupper(trimws(as.character(LABEL_ID)))
    )

missing_predictors <- setdiff(predictor_columns, names(environment))
if (length(missing_predictors) > 0L) {
    stop("Environmental table is missing predictors: ",
         paste(missing_predictors, collapse = ", "))
}

mean_or_na <- function(x) {
    if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# Reef names distinguish duplicated GBR feature IDs, so they are the primary
# join key. Grouping prevents any duplicated registry name from multiplying
# survey rows. Reef ID is used only as an explicit fallback for naming variants.
environment_by_name <- environment |>
    group_by(reef_name_key, event_year) |>
    summarise(
        across(all_of(c("lon", "lat", predictor_columns)), mean_or_na),
        environmental_feature_count = n(),
        cloud_product_version = paste(sort(unique(cloud_product_version)), collapse = ";"),
        cloud_aggregation = paste(sort(unique(cloud_aggregation)), collapse = ";"),
        .groups = "drop"
    )

environment_by_id <- environment |>
    group_by(reef_id_key, event_year) |>
    summarise(
        across(all_of(c("lon", "lat", predictor_columns)), mean_or_na),
        environmental_feature_count = n(),
        cloud_product_version = paste(sort(unique(cloud_product_version)), collapse = ";"),
        cloud_aggregation = paste(sort(unique(cloud_aggregation)), collapse = ";"),
        .groups = "drop"
    )

# Rebuild pre-event Acropora composition from the photo-transect observations.
# Exact baseline-year measurements are used first. Fallbacks are restricted to
# the baseline year or earlier so temporal assessment never sees future data.
load(acropora_workspace_file)
set.seed(202408L)
required_acropora_objects <- c('df.AIMS.full', 'dat.Ref.Clean')
if (!all(required_acropora_objects %in% ls())) {
    stop('The exploratory workspace lacks the AIMS composition objects')
}

acropora_depth <- df.AIMS.full |>
    filter(
        data_type == 'photo-transect', domain_category == 'reef',
        purpose == 'COMPOSITION', variable == 'HARD CORAL'
    ) |>
    mutate(
        reef_name_clean = gsub(
            ' REEF(S)?| ISLAND| IS', '', toupper(domain_name)
        )
    ) |>
    left_join(
        dat.Ref.Clean |>
            select(AIMS_REEF_NAME_Clean, ReefName),
        by = c('reef_name_clean' = 'AIMS_REEF_NAME_Clean')
    ) |>
    filter(!is.na(ReefName)) |>
    group_by(ReefName, report_year, depth) |>
    summarise(
        total_hard_coral = sum(mean, na.rm = TRUE),
        acropora_cover = sum(
            mean[grepl('Acropora', reefpage_category, ignore.case = TRUE)],
            na.rm = TRUE
        ),
        .groups = 'drop'
    ) |>
    mutate(
        prop_acropora_depth = if_else(
            total_hard_coral > 0,
            pmin(pmax(acropora_cover / total_hard_coral, 0), 1),
            NA_real_
        ),
        report_year = as.integer(report_year),
        depth = as.numeric(depth)
    ) |>
    select(ReefName, report_year, depth, prop_acropora_depth)

acropora_reef <- acropora_depth |>
    group_by(ReefName, report_year) |>
    summarise(
        prop_acropora_reef = mean(prop_acropora_depth, na.rm = TRUE),
        .groups = 'drop'
    ) |>
    mutate(
        prop_acropora_reef = if_else(
            is.nan(prop_acropora_reef), NA_real_, prop_acropora_reef
        ),
        reef_name_key = toupper(trimws(as.character(ReefName)))
    )

reef_coordinates <- environment |>
    group_by(reef_name_key) |>
    summarise(
        lon = mean_or_na(lon),
        lat = mean_or_na(lat),
        .groups = 'drop'
    )

acropora_reef_geo <- acropora_reef |>
    left_join(reef_coordinates, by = 'reef_name_key')

distance_km <- function(lon1, lat1, lon2, lat2) {
    radius <- 6371.0088
    lon1 <- lon1 * pi / 180
    lat1 <- lat1 * pi / 180
    lon2 <- lon2 * pi / 180
    lat2 <- lat2 * pi / 180
    delta_lon <- lon2 - lon1
    delta_lat <- lat2 - lat1
    angle <- sin(delta_lat / 2)^2 +
        cos(lat1) * cos(lat2) * sin(delta_lon / 2)^2
    2 * radius * atan2(sqrt(angle), sqrt(pmax(1 - angle, 0)))
}

estimate_past_acropora <- function(reef_name, baseline_year, lon, lat) {
    same_reef <- acropora_reef_geo |>
        filter(
            ReefName == reef_name,
            report_year < baseline_year,
            report_year >= baseline_year - 2L,
            is.finite(prop_acropora_reef)
        ) |>
        arrange(desc(report_year))
    if (nrow(same_reef) > 0L) {
        return(tibble(
            prop_acropora_fallback = same_reef$prop_acropora_reef[[1]],
            acropora_fallback_source = 'past_same_reef',
            acropora_latest_source_year = same_reef$report_year[[1]],
            acropora_source_reefs = 1L
        ))
    }

    candidates <- acropora_reef_geo |>
        filter(
            report_year <= baseline_year,
            report_year >= baseline_year - 2L,
            is.finite(prop_acropora_reef),
            is.finite(lon), is.finite(lat)
        )
    if (nrow(candidates) == 0L || !is.finite(lon) || !is.finite(lat)) {
        return(tibble(
            prop_acropora_fallback = NA_real_,
            acropora_fallback_source = 'fold_median_required',
            acropora_latest_source_year = NA_integer_,
            acropora_source_reefs = 0L
        ))
    }

    candidates <- candidates |>
        mutate(
            separation_km = distance_km(
                .env$lon, .env$lat, .data$lon, .data$lat
            ),
            year_lag = baseline_year - report_year
        ) |>
        arrange(separation_km, year_lag) |>
        slice_head(n = 8L) |>
        mutate(weight = 1 / ((separation_km + 5)^2 * (year_lag + 1)))

    tibble(
        prop_acropora_fallback = weighted.mean(
            candidates$prop_acropora_reef, candidates$weight
        ),
        acropora_fallback_source = 'past_only_spatiotemporal_idw',
        acropora_latest_source_year = max(candidates$report_year),
        acropora_source_reefs = n_distinct(candidates$ReefName)
    )
}

model_rows <- outcomes |>
    filter(event_year %in% model_years) |>
    left_join(
        rename_with(environment_by_name, ~ paste0(.x, "_name"),
                    -c(reef_name_key, event_year)),
        by = c("reef_name_key", "event_year")
    ) |>
    left_join(
        rename_with(environment_by_id, ~ paste0(.x, "_id"),
                    -c(reef_id_key, event_year)),
        by = c("reef_id_key", "event_year")
    )

# Coalesce the exact-name join with the reef-ID fallback without hiding how a
# row was matched. Predictor missing values are retained for fold-local handling
# during model fitting; they are never globally imputed here.
for (column in c("lon", "lat", predictor_columns, "environmental_feature_count")) {
    model_rows[[column]] <- coalesce(
        model_rows[[paste0(column, "_name")]],
        model_rows[[paste0(column, "_id")]]
    )
}
for (column in c("cloud_product_version", "cloud_aggregation")) {
    model_rows[[column]] <- coalesce(
        model_rows[[paste0(column, "_name")]],
        model_rows[[paste0(column, "_id")]]
    )
}

model_rows <- model_rows |>
    mutate(
        environment_match_method = case_when(
            !is.na(lon_name) ~ "exact_reef_name",
            is.na(lon_name) & !is.na(lon_id) ~ "reef_id_fallback",
            TRUE ~ "unmatched"
        ),
        region_block = case_when(
            SECTOR %in% c("CG", "CL", "PB") ~ "Northern GBR",
            SECTOR %in% c("CA", "CU", "IN", "TO", "WH") ~ "Central GBR",
            SECTOR %in% c("CB", "PO", "SW") ~ "Southern GBR",
            TRUE ~ NA_character_
        ),
        analysis_set = "model_development",
        temporal_fold = event_year
    ) |>
    select(
        -any_of(c(
            paste0(c(
                "lon", "lat", predictor_columns, "environmental_feature_count",
                "cloud_product_version", "cloud_aggregation"
            ), "_name"),
            paste0(c(
                "lon", "lat", predictor_columns, "environmental_feature_count",
                "cloud_product_version", "cloud_aggregation"
            ), "_id")
        )),
        -reef_name_key, -reef_id_key
    )

if (any(model_rows$environment_match_method == "unmatched")) {
    stop("Modern outcome rows remain unmatched to the environmental table")
}
model_rows <- model_rows |>
    left_join(
        acropora_depth,
        by = c(
            'ReefName', 'baseline_report_year' = 'report_year', 'depth'
        )
    ) |>
    left_join(
        acropora_reef |>
            select(ReefName, report_year, prop_acropora_reef),
        by = c('ReefName', 'baseline_report_year' = 'report_year')
    ) |>
    mutate(
        prop_acropora_pre = coalesce(
            prop_acropora_depth, prop_acropora_reef
        ),
        acropora_source = case_when(
            is.finite(prop_acropora_depth) ~ 'observed_depth_year',
            is.finite(prop_acropora_reef) ~ 'observed_reef_year',
            TRUE ~ NA_character_
        ),
        acropora_latest_source_year = if_else(
            !is.na(acropora_source), baseline_report_year, NA_integer_
        ),
        acropora_source_reefs = if_else(
            !is.na(acropora_source), 1L, NA_integer_
        )
    )

missing_acropora <- which(!is.finite(model_rows$prop_acropora_pre))
if (length(missing_acropora) > 0L) {
    fallbacks <- bind_rows(lapply(missing_acropora, function(row_index) {
        estimate_past_acropora(
            model_rows$ReefName[[row_index]],
            model_rows$baseline_report_year[[row_index]],
            model_rows$lon[[row_index]],
            model_rows$lat[[row_index]]
        ) |>
            mutate(row_index = row_index)
    }))
    model_rows$prop_acropora_pre[fallbacks$row_index] <-
        fallbacks$prop_acropora_fallback
    model_rows$acropora_source[fallbacks$row_index] <-
        fallbacks$acropora_fallback_source
    model_rows$acropora_latest_source_year[fallbacks$row_index] <-
        fallbacks$acropora_latest_source_year
    model_rows$acropora_source_reefs[fallbacks$row_index] <-
        fallbacks$acropora_source_reefs
}

model_rows <- model_rows |>
    mutate(
        acropora_cover_pre = prop_acropora_pre * observed_pre_cover,
        acropora_interpolated = !acropora_source %in%
            c('observed_depth_year', 'observed_reef_year')
    ) |>
    select(-prop_acropora_depth, -prop_acropora_reef)

if (any(
    model_rows$acropora_latest_source_year > model_rows$baseline_report_year,
    na.rm = TRUE
)) {
    stop('Pre-event Acropora reconstruction used a future observation')
}
if (any(
    model_rows$prop_acropora_pre < 0 | model_rows$prop_acropora_pre > 1,
    na.rm = TRUE
)) {
    stop('Pre-event Acropora proportions fall outside 0--1')
}

expected_modern_rows <- sum(
    outcomes$event_year %in% model_years
)
if (nrow(model_rows) != expected_modern_rows ||
    anyDuplicated(paste(model_rows$programme_key, model_rows$source_observation_id))) {
    stop("The environmental join changed the number or uniqueness of outcome rows")
}
if (any(is.na(model_rows$region_block))) {
    stop("Modern outcome rows have no spatial region block")
}

# Greedy allocation keeps all observations from a reef together while making
# the five reef folds reasonably similar in size. The seed above only resolves
# ties; saved assignments make later fits fully reproducible.
make_balanced_reef_folds <- function(data, folds = 5L) {
    reef_sizes <- data |>
        count(ReefID, name = "observations") |>
        mutate(tie_break = runif(n())) |>
        arrange(desc(observations), tie_break)

    fold_load <- rep(0L, min(folds, nrow(reef_sizes)))
    reef_sizes$reef_fold <- NA_integer_
    for (i in seq_len(nrow(reef_sizes))) {
        available <- which(fold_load == min(fold_load))
        chosen <- available[1]
        reef_sizes$reef_fold[i] <- chosen
        fold_load[chosen] <- fold_load[chosen] + reef_sizes$observations[i]
    }
    select(reef_sizes, ReefID, reef_fold)
}

model_rows$reef_fold <- NA_integer_
for (programme_name in names(outcome_files)) {
    development <- model_rows |>
        filter(programme_key == .env$programme_name)
    reef_map <- make_balanced_reef_folds(development)
    row_index <- which(
        model_rows$programme_key == programme_name
    )
    model_rows$reef_fold[row_index] <- reef_map$reef_fold[
        match(model_rows$ReefID[row_index], reef_map$ReefID)
    ]
}

# Create a compact manifest. Each row describes one outer assessment fold;
# models use all other development rows from that programme as analysis data.
manifest <- tibble()
for (programme_name in names(outcome_files)) {
    programme_rows <- filter(model_rows, programme_key == .env$programme_name)
    development <- programme_rows

    for (fold in sort(unique(development$reef_fold))) {
        assessment <- filter(development, reef_fold == fold)
        analysis <- filter(development, reef_fold != fold)
        manifest <- bind_rows(manifest, tibble(
            programme_key = programme_name,
            scheme = "reef_blocked_5fold",
            fold = as.character(fold),
            analysis_rows = nrow(analysis),
            assessment_rows = nrow(assessment),
            assessment_reefs = n_distinct(assessment$ReefID),
            assessment_events = paste(sort(unique(assessment$event_year)), collapse = ";"),
            assessment_regions = paste(sort(unique(assessment$region_block)), collapse = ";"),
            leakage_check = length(intersect(analysis$ReefID, assessment$ReefID)) == 0L
        ))
    }

    for (region in sort(unique(development$region_block))) {
        assessment <- filter(development, region_block == region)
        analysis <- filter(development, region_block != region)
        manifest <- bind_rows(manifest, tibble(
            programme_key = programme_name,
            scheme = "region_blocked",
            fold = region,
            analysis_rows = nrow(analysis),
            assessment_rows = nrow(assessment),
            assessment_reefs = n_distinct(assessment$ReefID),
            assessment_events = paste(sort(unique(assessment$event_year)), collapse = ";"),
            assessment_regions = region,
            leakage_check = length(intersect(analysis$region_block, assessment$region_block)) == 0L
        ))
    }

    for (year in sort(unique(development$event_year))) {
        assessment <- filter(development, event_year == year)
        analysis <- filter(development, event_year != year)
        manifest <- bind_rows(manifest, tibble(
            programme_key = programme_name,
            scheme = "leave_one_event_out",
            fold = as.character(year),
            analysis_rows = nrow(analysis),
            assessment_rows = nrow(assessment),
            assessment_reefs = n_distinct(assessment$ReefID),
            assessment_events = as.character(year),
            assessment_regions = paste(sort(unique(assessment$region_block)), collapse = ";"),
            leakage_check = length(intersect(analysis$event_year, assessment$event_year)) == 0L
        ))
    }

}

if (nrow(manifest) == 0L || any(!manifest$leakage_check)) {
    stop("Validation manifest failed its leakage checks")
}
if (any(manifest$analysis_rows == 0L | manifest$assessment_rows == 0L)) {
    empty_folds <- manifest |>
        filter(analysis_rows == 0L | assessment_rows == 0L) |>
        select(programme_key, scheme, fold, analysis_rows, assessment_rows)
    stop(
        "At least one validation fold has an empty analysis or assessment set:\n",
        paste(capture.output(print(empty_folds)), collapse = "\n")
    )
}
if (any(is.na(model_rows$reef_fold)) || any(is.na(model_rows$temporal_fold))) {
    stop("At least one modelling row has no reef or temporal fold assignment")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
for (programme_name in names(outcome_files)) {
    output <- model_rows |>
        filter(programme_key == .env$programme_name) |>
        arrange(analysis_set, event_year, ReefID, survey_date)
    saveRDS(output, file.path(output_dir, paste0("validation_rows_", programme_name, ".rds")))
    write_csv(output, file.path(output_dir, paste0("validation_rows_", programme_name, ".csv")), na = "")
}
write_csv(manifest, file.path(output_dir, "validation_split_manifest.csv"))

acropora_coverage <- model_rows |>
    count(programme_key, acropora_source, name = 'rows') |>
    group_by(programme_key) |>
    mutate(proportion = rows / sum(rows)) |>
    ungroup()
write_csv(
    acropora_coverage,
    file.path(output_dir, 'pre_event_acropora_coverage.csv')
)

summary <- model_rows |>
    count(programme_key, analysis_set, name = "rows") |>
    left_join(
        model_rows |>
            group_by(programme_key, analysis_set) |>
            summarise(
                reefs = n_distinct(ReefID),
                events = paste(sort(unique(event_year)), collapse = ";"),
                .groups = "drop"
            ),
        by = c("programme_key", "analysis_set")
    )
print(summary)
