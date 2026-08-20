# Harmonise growth-adjusted coral-cover outcomes across AIMS survey programmes.
#
# Run directly to rebuild the processed outcome tables from the exploratory
# workspace, or source this file from 01_DHW_Mortality_Exploratory.qmd.

suppressPackageStartupMessages({
    library(dplyr)
    library(lubridate)
    library(readr)
})

bleaching_event_years <- c(1998L, 2002L, 2016L, 2017L, 2020L, 2022L, 2024L)

decimal_year_to_date <- function(x) {
    x <- as.numeric(x)
    year <- floor(x)
    fraction <- x - year
    start <- ymd(paste0(year, "-01-01"))
    next_start <- ymd(paste0(year + 1L, "-01-01"))
    start + round(fraction * as.numeric(next_start - start))
}

add_baseline_dates <- function(outcomes, raw_surveys, group_columns) {
    if (all(c("Baseline.Date", "Interval.Days") %in% names(outcomes))) {
        return(outcomes)
    }

    missing_group_columns <- setdiff(c("id", "date", "report_year", group_columns), names(raw_surveys))
    if (length(missing_group_columns) > 0L) {
        stop("Cannot reconstruct baseline dates; raw survey columns are missing: ",
             paste(missing_group_columns, collapse = ", "))
    }

    baseline_lookup <- raw_surveys |>
        mutate(Survey.Date = decimal_year_to_date(date)) |>
        arrange(across(all_of(group_columns)), report_year, Survey.Date) |>
        group_by(across(all_of(group_columns))) |>
        mutate(
            Baseline.Date = lag(Survey.Date),
            Interval.Days = as.integer(Survey.Date - Baseline.Date)
        ) |>
        ungroup() |>
        select(id, Baseline.Date, Interval.Days)

    left_join(outcomes, baseline_lookup, by = "id")
}

prepare_mortality_outcomes <- function(data, programme) {
    required <- c(
        "id", "project_code", "report_year", "DHWYear", "Date", "Lag",
        "mean", "Prev.HC", "Pred.Grth2", "Abs.Change", "ReefName"
    )
    missing_required <- setdiff(required, names(data))
    if (length(missing_required) > 0L) {
        stop(programme, " data are missing required columns: ",
             paste(missing_required, collapse = ", "))
    }
    if (!all(c("Baseline.Date", "Interval.Days") %in% names(data))) {
        stop(programme, " data require Baseline.Date and Interval.Days")
    }

    output <- data |>
        ungroup() |>
        mutate(
            programme = programme,
            source_observation_id = as.character(id),
            survey_date = as.Date(Date),
            baseline_survey_date = as.Date(Baseline.Date),
            event_year = as.integer(DHWYear),
            baseline_report_year = as.integer(report_year - Lag),
            event_window_start = ymd(paste0(event_year, "-05-01")),
            event_window_end = ymd(paste0(event_year + 1L, "-04-30")),
            event_window_ok = survey_date >= event_window_start & survey_date <= event_window_end,
            interval_days = as.integer(Interval.Days),
            baseline_within_24_months = !is.na(baseline_survey_date) &
                baseline_survey_date >= (survey_date %m-% years(2)) &
                baseline_survey_date < survey_date,
            observed_pre_cover = as.numeric(Prev.HC),
            expected_cover_no_event = as.numeric(Prev.HC + Pred.Grth2),
            post_cover = as.numeric(mean),
            expected_growth = as.numeric(Pred.Grth2),
            mortality_prop_raw = (observed_pre_cover - post_cover) /
                observed_pre_cover,
            mortality_prop = pmin(pmax(mortality_prop_raw, 0), 1),
            cover_loss_pp = 100 * pmax(observed_pre_cover - post_cover, 0),
            observed_cover_change_pp = 100 * (post_cover - observed_pre_cover),
            growth_adjusted_mortality_prop_raw =
                (expected_cover_no_event - post_cover) / expected_cover_no_event,
            growth_adjusted_mortality_prop = pmin(pmax(
                growth_adjusted_mortality_prop_raw, 0
            ), 1),
            growth_adjusted_cover_loss_pp =
                100 * pmax(expected_cover_no_event - post_cover, 0),
            growth_adjusted_change_pp = 100 * as.numeric(Abs.Change),
            cover_gain_clamped_to_zero = is.finite(mortality_prop_raw) & mortality_prop_raw < 0,
            complete_loss = is.finite(mortality_prop_raw) & mortality_prop_raw >= 1,
            boundary_zero = is.finite(mortality_prop) & mortality_prop == 0,
            boundary_one = is.finite(mortality_prop) & mortality_prop == 1,
            target_event_year = event_year %in% bleaching_event_years,
            valid_cover_denominator = is.finite(observed_pre_cover) & observed_pre_cover > 0,
            eligible_primary = target_event_year & event_window_ok &
                baseline_within_24_months & Lag %in% c(1, 2) &
                valid_cover_denominator & is.finite(post_cover) &
                is.finite(mortality_prop_raw)
        ) |>
        arrange(event_year, ReefName, survey_date, source_observation_id)

    output
}

select_outcome_columns <- function(data) {
    core <- c(
        'baseline_report_year', 'Lag',
        "programme", "source_observation_id", "ReefName", "ReefID", "Region",
        "SECTOR", "report_year", "event_year", "survey_date",
        "baseline_survey_date", "interval_days", "event_window_start",
        "event_window_end", "depth", "reef_zone", "reefpage_category",
        "observed_pre_cover", "expected_growth", "expected_cover_no_event",
        "post_cover", "mortality_prop_raw", "mortality_prop", "cover_loss_pp",
        "observed_cover_change_pp", "growth_adjusted_mortality_prop_raw",
        "growth_adjusted_mortality_prop", "growth_adjusted_cover_loss_pp",
        "growth_adjusted_change_pp",
        "cover_gain_clamped_to_zero", "complete_loss", "boundary_zero",
        "boundary_one", "MaxDHW.mean", "DISTURBANCE_TYPE", "storm_name",
        "description"
    )
    select(data, any_of(core))
}

validate_mortality_outcomes <- function(data, programme) {
    failures <- character()
    if (nrow(data) == 0L) failures <- c(failures, 'contains no eligible rows')
    if (anyDuplicated(data$source_observation_id)) {
        failures <- c(failures, 'has duplicated source observation IDs')
    }
    if (any(!is.finite(data$mortality_prop)) ||
        any(data$mortality_prop < 0 | data$mortality_prop > 1)) {
        failures <- c(failures, 'has mortality proportions outside 0--1')
    }
    if (any(data$survey_date < data$event_window_start |
        data$survey_date > data$event_window_end)) {
        failures <- c(failures, 'has surveys outside the May--April event window')
    }
    if (any(is.na(data$baseline_survey_date)) ||
        any(data$baseline_survey_date >= data$survey_date) ||
        any((data$baseline_survey_date %m+% years(2)) < data$survey_date)) {
        failures <- c(failures, 'has an invalid or older-than-24-month baseline')
    }
    expected <- pmin(pmax(
        (data$observed_pre_cover - data$post_cover) /
            data$observed_pre_cover,
        0
    ), 1)
    if (!isTRUE(all.equal(data$mortality_prop, expected, tolerance = 1e-12))) {
        failures <- c(failures, 'does not reproduce the documented estimand')
    }
    if (length(failures) > 0L) {
        stop(programme, ': ', paste(failures, collapse = '; '))
    }
    invisible(data)
}

summarise_outcome_build <- function(all_rows, eligible_rows) {
    tibble(
        programme = first(all_rows$programme),
        input_rows = nrow(all_rows),
        target_event_rows = sum(all_rows$target_event_year, na.rm = TRUE),
        output_rows = nrow(eligible_rows),
        excluded_non_target_event = sum(!all_rows$target_event_year, na.rm = TRUE),
        excluded_outside_event_window = sum(all_rows$target_event_year & !all_rows$event_window_ok, na.rm = TRUE),
        excluded_baseline_over_24_months = sum(
            all_rows$target_event_year & all_rows$event_window_ok &
                !all_rows$baseline_within_24_months,
            na.rm = TRUE
        ),
        excluded_invalid_outcome = sum(
            all_rows$target_event_year & all_rows$event_window_ok &
                all_rows$baseline_within_24_months &
                (!all_rows$valid_cover_denominator | !is.finite(all_rows$post_cover) |
                    !is.finite(all_rows$mortality_prop_raw)),
            na.rm = TRUE
        ),
        zeros_retained = sum(eligible_rows$boundary_zero, na.rm = TRUE),
        ones_retained = sum(eligible_rows$boundary_one, na.rm = TRUE),
        cover_gains_set_to_zero = sum(eligible_rows$cover_gain_clamped_to_zero, na.rm = TRUE)
    )
}

build_mortality_outcome_tables <- function(
    dat_manta,
    dat_benthic,
    raw_manta = NULL,
    raw_benthic = NULL,
    output_dir = "data/processed",
    write_outputs = TRUE
) {
    if (!all(c("Baseline.Date", "Interval.Days") %in% names(dat_manta))) {
        if (is.null(raw_manta)) stop("raw_manta is required to recover baseline dates")
        dat_manta <- add_baseline_dates(
            dat_manta, raw_manta,
            c("Reef_Name", "depth", "reefpage_category")
        )
    }
    if (!all(c("Baseline.Date", "Interval.Days") %in% names(dat_benthic))) {
        if (is.null(raw_benthic)) stop("raw_benthic is required to recover baseline dates")
        dat_benthic <- add_baseline_dates(
            dat_benthic, raw_benthic,
            c("Reef_Name", "depth", "reef_zone", "reefpage_category")
        )
    }

    all_tables <- list(
        manta = prepare_mortality_outcomes(dat_manta, "LTMP manta tow"),
        ltmp = prepare_mortality_outcomes(
            filter(dat_benthic, project_code == "LTMP"), "LTMP benthic"
        ),
        mmp = prepare_mortality_outcomes(
            filter(dat_benthic, project_code == "MMP"), "MMP inshore benthic"
        )
    )
    eligible_tables <- lapply(all_tables, function(x) {
        x |>
            filter(eligible_primary) |>
            select_outcome_columns()
    })
    invisible(Map(validate_mortality_outcomes, eligible_tables, c(
        'LTMP manta tow', 'LTMP benthic', 'MMP inshore benthic'
    )))
    manifest <- bind_rows(Map(summarise_outcome_build, all_tables, eligible_tables))

    if (write_outputs) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        for (name in names(eligible_tables)) {
            saveRDS(
                eligible_tables[[name]],
                file.path(output_dir, paste0("mortality_outcomes_", name, ".rds"))
            )
            write_csv(
                eligible_tables[[name]],
                file.path(output_dir, paste0("mortality_outcomes_", name, ".csv")),
                na = ""
            )
        }
        write_csv(manifest, file.path(output_dir, "mortality_outcomes_manifest.csv"))
    }

    list(
        manta = eligible_tables$manta,
        ltmp = eligible_tables$ltmp,
        mmp = eligible_tables$mmp,
        manifest = manifest
    )
}

if (sys.nframe() == 0L) {
    workspace <- "data/processed/01_exploratory_workspace.RData"
    if (!file.exists(workspace)) stop("Missing exploratory workspace: ", workspace)
    load(workspace)
    required_objects <- c(
        "df.MANT", "dat.mod.bent.hc", "df.AIMS.Mant", "df.AIMS.Bent",
        "Disturbances_mant"
    )
    missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
    if (length(missing_objects) > 0L) {
        stop("Exploratory workspace is missing: ", paste(missing_objects, collapse = ", "))
    }
    # Recreate the manta candidate set without the former abs(Change) > 0.01
    # threshold, which selectively removed near-zero outcomes.
    manta_source <- df.MANT |>
        ungroup() |>
        filter(Lag < 3, project_code != "MMP") |>
        left_join(
            select(
                Disturbances_mant, ReefName, report_year, DISTURBANCE_TYPE,
                storm_name, description
            ),
            by = c("ReefName", "report_year")
        ) |>
        filter(
            is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n"),
            report_year %in% c(1999, 2003, 2017, 2018, 2021, 2023, 2024, 2025)
        )

    result <- build_mortality_outcome_tables(
        manta_source,
        dat.mod.bent.hc,
        raw_manta = df.AIMS.Mant,
        raw_benthic = df.AIMS.Bent
    )
    print(result$manifest)
}
