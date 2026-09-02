# Cause-aware competing-hazard candidate for relative coral mortality.
#
# The selected INLA response is refitted without cyclone/COTS predictors and
# without explicitly labelled non-thermal rows in its thermal training set.
# Separate two-part INLA hazards learn COTS and cyclone loss from the longer
# annual cover-transition record. Predictions are combined on the mortality
# scale as 1 - (1 - thermal) * (1 - COTS) * (1 - cyclone).

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(INLA)
    library(patchwork)
    library(readr)
    library(stringr)
    library(tidyr)
})

Sys.setenv(INLA_ST_RUN = '0')
source('scripts/fit_inla_spatiotemporal_screen.R')
source('scripts/model_registry.R')
source('scripts/model_diagnostics.R')

root <- project_root()
out_dir <- file.path(root, 'output', 'cause_aware_competing_hazards')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# 1. Thermal/freshwater INLA component
# -------------------------------------------------------------------------

thermal_core_id <-
    'persistent_rw1_local_dhw_decomposed_hazards_freshwater_partial_pool'
selected_id <- thermal_core_id
selected_row <- candidates |> filter(candidate == selected_id)
stopifnot(nrow(selected_row) == 1L)

selected_terms <- candidate_shared_terms(selected_row)
cause_terms <- c(
    'log1p_cyc_interval_maxHrs4mw_z',
    'tc_interval_proximity100_z', 'tc_interval_wind_distance_index_z',
    'log1p_cot_interval_idw_max_z', 'cots_outbreak_probability_z'
)
thermal_terms <- setdiff(selected_terms, cause_terms)

thermal_candidate <- selected_row |>
    mutate(
        candidate = 'test_cause_aware_competing_hazards',
        feature_set = 'thermal_freshwater'
    )

base_candidate_shared_terms <- candidate_shared_terms
candidate_shared_terms <- function(candidate_row) {
    if (identical(
        as.character(candidate_row$feature_set), 'thermal_freshwater'
    )) thermal_terms else base_candidate_shared_terms(candidate_row)
}

is_explicit_nonthermal <- function(rows) {
    coalesce(rows$disturbance_has_cots, FALSE) |
        coalesce(rows$disturbance_has_cyclone, FALSE) |
        coalesce(rows$disturbance_has_flood, FALSE)
}

base_fit_candidate <- fit_candidate
fit_candidate <- function(analysis, assessment, candidate,
                          compute_criteria = FALSE, seed = 1L) {
    if (identical(
        as.character(candidate$feature_set), 'thermal_freshwater'
    )) {
        analysis <- analysis[!is_explicit_nonthermal(analysis), , drop = FALSE]
    }
    base_fit_candidate(
        analysis, assessment, candidate,
        compute_criteria = compute_criteria, seed = seed
    )
}

full_thermal_path <- file.path(out_dir, 'full_thermal_component.rds')
if (file.exists(full_thermal_path)) {
    full_thermal <- readRDS(full_thermal_path)
} else {
    full_thermal <- fit_candidate(
        data, data[1, , drop = FALSE], thermal_candidate,
        compute_criteria = TRUE, seed = 20260903L
    )
    saveRDS(full_thermal, full_thermal_path)
}
thermal_summary <- extract_fit_summary(full_thermal, thermal_candidate)
write_csv(
    thermal_summary$criteria,
    file.path(out_dir, 'thermal_component_criteria.csv')
)
write_csv(
    thermal_summary$fixed,
    file.path(out_dir, 'thermal_component_fixed_effects.csv')
)

# -------------------------------------------------------------------------
# 2. Annual cause-labelled training data
# -------------------------------------------------------------------------

cots_hindcast_annual <- read_csv(
    file.path(root, 'data', 'gbrPredsAdj_20262408.csv'),
    show_col_types = FALSE
) |>
    transmute(
        ReefName = reefName, event_year = as.integer(year),
        cots_outbreak_probability = outbrProb
    )

cyclone_annual <- read_csv(
    file.path(root, 'data', 'processed', 'bom_cyclone_reef_year.csv'),
    show_col_types = FALSE
) |>
    transmute(
        ReefID, event_year,
        lon = reef_longitude, lat = reef_latitude,
        tc_proximity100 = exp(-pmin(tc_min_distance_km, 1000) / 100),
        tc_wind_distance_index = pmax(tc_nearest_max_wind_ms - 17, 0) *
            tc_proximity100
    )

annual <- read_csv(
    file.path(root, 'data', 'processed', 'annual_coral_transitions.csv'),
    show_col_types = FALSE
) |>
    left_join(
        cots_hindcast_annual,
        by = c('ReefName', 'event_year'), relationship = 'many-to-one'
    ) |>
    left_join(
        cyclone_annual,
        by = c('ReefID', 'event_year'), relationship = 'many-to-one'
    ) |>
    mutate(
        relative_loss = pmin(
            pmax(-cover_change_pp / (100 * pmax(pre_cover, 0.02)), 0),
            0.999
        ),
        cots_label = coalesce(disturbance_has_cots, FALSE) |
            str_detect(
                str_to_lower(coalesce(disturbance_text, '')),
                'cots|crown-of-thorns'
            ),
        cyclone_label = coalesce(disturbance_has_cyclone, FALSE) |
            coalesce(disturbance_has_flood, FALSE) |
            str_detect(
                str_to_lower(coalesce(disturbance_text, '')),
                'cyclone|storm|flood'
            ),
        cots_positive = cots_label & relative_loss > 0.02,
        cyclone_positive = cyclone_label & relative_loss > 0.02,
        cots_probability = cots_outbreak_probability,
        cots_weighted_excess = cots_probability *
            log1p(pmax(cot_idwmeanpertow - 0.22, 0)),
        wave_log = log1p(pmax(cyc_maxHrs4mw, 0)),
        wind_distance = tc_wind_distance_index,
        rain_log = log_coastal_rain30,
        acropora = prop_acropora_pre,
        is_manta = as.numeric(programme_key == 'manta'),
        is_mmp = as.numeric(programme_key == 'mmp')
    )

make_assessment_features <- function(rows) {
    rows |>
        transmute(
            source_observation_id = as.character(source_observation_id),
            programme_key, ReefID, ReefName, event_year,
            observed_mortality = mortality_prop,
            cots_probability = cots_outbreak_probability,
            cots_weighted_excess = cots_outbreak_probability *
                log1p(pmax(cot_interval_idw_max - 0.22, 0)),
            wave_log = log1p(pmax(cyc_interval_maxHrs4mw, 0)),
            wind_distance = tc_interval_wind_distance_index,
            rain_log = log_coastal_rain30,
            pre_cover = observed_pre_cover,
            acropora = prop_acropora_pre,
            lon, lat,
            is_manta = as.numeric(programme_key == 'manta'),
            is_mmp = as.numeric(programme_key == 'mmp'),
            cots_interval_max = cot_interval_idw_max,
            cots_outbreak_probability,
            cyclone_wave_hours = cyc_interval_maxHrs4mw,
            cyclone_name = tc_interval_peak_name,
            ann_maxdhw
        )
}

cots_features <- c(
    'cots_probability', 'cots_weighted_excess', 'pre_cover',
    'acropora', 'lon', 'lat', 'is_manta', 'is_mmp'
)
cyclone_features <- c(
    'wave_log', 'wind_distance', 'rain_log', 'pre_cover',
    'acropora', 'lon', 'lat', 'is_manta', 'is_mmp'
)

prepare_cause_fold <- function(training, assessment, features) {
    for (feature in features) {
        training_value <- training[[feature]]
        replacement <- median(training_value[is.finite(training_value)],
                              na.rm = TRUE)
        if (!is.finite(replacement)) replacement <- 0
        training_value[!is.finite(training_value)] <- replacement
        assessment_value <- assessment[[feature]]
        assessment_value[!is.finite(assessment_value)] <- replacement
        centre <- mean(training_value)
        spread <- sd(training_value)
        if (!is.finite(spread) || spread < 1e-8) spread <- 1
        training[[paste0(feature, '_z')]] <-
            (training_value - centre) / spread
        assessment[[paste0(feature, '_z')]] <-
            (assessment_value - centre) / spread
    }
    list(training = training, assessment = assessment)
}

fit_cause_hurdle <- function(training, assessment, cause, features) {
    prepared <- prepare_cause_fold(training, assessment, features)
    training <- prepared$training
    assessment <- prepared$assessment
    feature_terms <- paste0(features, '_z')
    formula <- as.formula(paste(
        'response ~ 1 +', paste(feature_terms, collapse = ' + ')
    ))

    occurrence_response <- training[[paste0(cause, '_positive')]]
    occurrence_data <- bind_rows(
        training |>
            transmute(
                response = as.numeric(occurrence_response),
                across(all_of(feature_terms))
            ),
        assessment |>
            transmute(
                response = NA_real_, across(all_of(feature_terms))
            )
    )
    occurrence_fit <- inla(
        formula, family = 'binomial', data = occurrence_data,
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
        verbose = FALSE
    )
    new_index <- seq.int(nrow(training) + 1L, nrow(occurrence_data))
    occurrence_prediction <-
        occurrence_fit$summary.fitted.values$mean[new_index]

    positive_training <- training |>
        filter(.data[[paste0(cause, '_positive')]]) |>
        transmute(
            response = pmin(pmax(relative_loss, 0.001), 0.999),
            across(all_of(feature_terms))
        )
    magnitude_data <- bind_rows(
        positive_training,
        assessment |>
            transmute(
                response = NA_real_, across(all_of(feature_terms))
            )
    )
    magnitude_fit <- inla(
        formula, family = 'beta', data = magnitude_data,
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE),
        verbose = FALSE
    )
    magnitude_index <- seq.int(
        nrow(positive_training) + 1L, nrow(magnitude_data)
    )
    magnitude_prediction <-
        magnitude_fit$summary.fitted.values$mean[magnitude_index]

    fixed <- bind_rows(
        as_tibble(occurrence_fit$summary.fixed, rownames = 'term') |>
            mutate(component = paste0(cause, '_occurrence')),
        as_tibble(magnitude_fit$summary.fixed, rownames = 'term') |>
            mutate(component = paste0(cause, '_magnitude'))
    )
    list(
        prediction = occurrence_prediction * magnitude_prediction,
        occurrence = occurrence_prediction,
        magnitude = magnitude_prediction,
        fixed = fixed
    )
}

# -------------------------------------------------------------------------
# 3. Matched validation and competing-risk combination
# -------------------------------------------------------------------------

baseline_predictions <- read_csv(
    file.path(root, 'output', 'inla_hinge_sensitivity', 'cv_predictions.csv'),
    show_col_types = FALSE
) |>
    filter(candidate == 'selected_hinge4_plus8') |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )

candidate_predictions <- tibble()
thermal_only_predictions <- tibble()

for (scheme in c('leave_one_event_out', 'reef_blocked_5fold')) {
    folds <- if (scheme == 'leave_one_event_out') {
        event_years
    } else sort(unique(data$joint_reef_fold))

    for (fold in folds) {
        if (scheme == 'leave_one_event_out') {
            assessment <- data |> filter(event_year == fold)
            analysis <- data |> filter(event_year != fold)
            annual_training <- annual |>
                filter(!(
                    baseline_report_year <= fold & report_year >= fold
                ))
        } else {
            held_reefs <- unique(data$ReefID[data$joint_reef_fold == fold])
            assessment <- data |> filter(ReefID %in% held_reefs)
            analysis <- data |> filter(!ReefID %in% held_reefs)
            annual_training <- annual |> filter(!ReefID %in% held_reefs)
        }

        thermal_cache <- file.path(
            out_dir, paste0('thermal_', scheme, '_', fold, '.rds')
        )
        if (file.exists(thermal_cache)) {
            thermal_result <- readRDS(thermal_cache)
        } else {
            fit <- fit_candidate(
                analysis, assessment, thermal_candidate,
                compute_criteria = FALSE,
                seed = 20260903L + as.integer(fold)
            )
            thermal_result <- list(predictions = fit$predictions)
            saveRDS(thermal_result, thermal_cache)
        }

        assessment_features <- make_assessment_features(assessment)
        cots_result <- fit_cause_hurdle(
            annual_training, assessment_features, 'cots', cots_features
        )
        cyclone_result <- fit_cause_hurdle(
            annual_training, assessment_features, 'cyclone', cyclone_features
        )

        thermal_rows <- thermal_result$predictions |>
            mutate(source_observation_id = as.character(source_observation_id))
        combined <- thermal_rows |>
            select(
                programme_key, source_observation_id, ReefID, ReefName,
                event_year, observed_mortality, observed_occurrence,
                thermal_prediction = predicted_mortality,
                thermal_occurrence = predicted_occurrence
            ) |>
            left_join(
                assessment_features,
                by = c(
                    'programme_key', 'source_observation_id', 'ReefID',
                    'ReefName', 'event_year', 'observed_mortality'
                ), relationship = 'one-to-one'
            ) |>
            mutate(
                cots_prediction = cots_result$prediction,
                cyclone_prediction = cyclone_result$prediction,
                cots_occurrence = cots_result$occurrence,
                cyclone_occurrence = cyclone_result$occurrence,
                predicted_mortality = 1 -
                    (1 - thermal_prediction) *
                    (1 - cots_prediction) *
                    (1 - cyclone_prediction),
                predicted_occurrence = 1 -
                    (1 - thermal_occurrence) *
                    (1 - cots_occurrence) *
                    (1 - cyclone_occurrence),
                residual = observed_mortality - predicted_mortality,
                candidate = 'cause_aware_competing_hazards',
                scheme = scheme, fold = as.character(fold)
            )
        candidate_predictions <- bind_rows(candidate_predictions, combined)
        thermal_only_predictions <- bind_rows(
            thermal_only_predictions,
            thermal_rows |>
                transmute(
                    programme_key, source_observation_id, ReefID, ReefName,
                    event_year, observed_mortality, observed_occurrence,
                    predicted_mortality, predicted_occurrence,
                    residual = observed_mortality -
                        predicted_mortality,
                    candidate = 'cause_aware_thermal_only',
                    scheme = scheme, fold = as.character(fold)
                )
        )
    }
}

# The cause fits estimate loss conditional on the longer annual-transition
# record. For operational use, they are activated only when the corresponding
# independent exposure is present. The COTS threshold is the ecological
# outbreak definition; the cyclone gate is the damaging-wave evidence tested
# here and remains provisional until the replacement cyclone layer arrives.
gated_predictions <- candidate_predictions |>
    mutate(
        cots_gate = as.numeric(coalesce(cots_interval_max, 0) > 0.22),
        cyclone_gate = as.numeric(coalesce(cyclone_wave_hours, 0) > 20),
        predicted_mortality = 1 -
            (1 - thermal_prediction) *
            (1 - cots_prediction * cots_gate) *
            (1 - cyclone_prediction * cyclone_gate),
        predicted_occurrence = 1 -
            (1 - thermal_occurrence) *
            (1 - cots_occurrence * cots_gate) *
            (1 - cyclone_occurrence * cyclone_gate),
        residual = observed_mortality - predicted_mortality,
        candidate = 'cause_aware_exposure_gated_competing_hazards'
    )

baseline_metric_rows <- baseline_predictions |>
    transmute(
        programme_key, source_observation_id, ReefID, ReefName, event_year,
        observed_mortality, observed_occurrence,
        predicted_mortality, predicted_occurrence,
        residual = observed_mortality - predicted_mortality,
        candidate = 'selected_current', scheme, fold
    )
metric_rows <- bind_rows(
    baseline_metric_rows, thermal_only_predictions,
    candidate_predictions |>
        select(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, observed_mortality, observed_occurrence,
            predicted_mortality, predicted_occurrence,
            residual, candidate, scheme, fold
        ),
    gated_predictions |>
        select(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, observed_mortality, observed_occurrence,
            predicted_mortality, predicted_occurrence,
            residual, candidate, scheme, fold
        )
)
comparison <- metric_rows |>
    group_by(candidate, scheme) |>
    metric_summary() |>
    arrange(scheme, rmse)

write_csv(candidate_predictions, file.path(out_dir, 'cv_predictions.csv'))
write_csv(gated_predictions, file.path(out_dir, 'gated_cv_predictions.csv'))
write_csv(comparison, file.path(out_dir, 'model_comparison.csv'))

event_metrics <- metric_rows |>
    filter(scheme == 'leave_one_event_out') |>
    group_by(candidate, event_year) |>
    metric_summary()
write_csv(event_metrics, file.path(out_dir, 'event_metrics.csv'))

focal_pattern <- 'Gannett|Chinaman|Taylor|Rib Reef|Penrith|Daydream|Double Cone|Shute'
focal_predictions <- bind_rows(candidate_predictions, gated_predictions) |>
    filter(
        scheme == 'leave_one_event_out',
        str_detect(ReefName, regex(focal_pattern, ignore_case = TRUE))
    )
write_csv(focal_predictions, file.path(out_dir, 'focal_predictions.csv'))

# Full-data cause coefficient summaries use all annual transitions and the
# selected response rows only as prediction scaffolding.
full_assessment <- make_assessment_features(data)
full_cots <- fit_cause_hurdle(annual, full_assessment, 'cots', cots_features)
full_cyclone <- fit_cause_hurdle(
    annual, full_assessment, 'cyclone', cyclone_features
)
write_csv(
    bind_rows(full_cots$fixed, full_cyclone$fixed),
    file.path(out_dir, 'cause_fixed_effects.csv')
)

# -------------------------------------------------------------------------
# 4. Audit tables and paper/report figures
# -------------------------------------------------------------------------

cots_audit <- annual |>
    filter(
        cots_label,
        event_year %in% 2016:2021,
        str_detect(ReefName, regex(
            'Gannett|Chinaman|Taylor|Rib Reef', ignore_case = TRUE
        ))
    ) |>
    select(
        programme_key, ReefID, ReefName, event_year,
        baseline_report_year, report_year, pre_cover, post_cover,
        cover_change_pp, relative_loss, cot_idwmeanpertow,
        cots_outbreak_probability, disturbance_text
    ) |>
    arrange(ReefName, event_year)
write_csv(cots_audit, file.path(out_dir, 'focal_cots_annual_losses.csv'))

cots_high_pressure_low_loss <- annual |>
    filter(cot_idwmeanpertow > 0.22, relative_loss <= 0.02) |>
    arrange(desc(cot_idwmeanpertow)) |>
    select(
        ReefID, ReefName, event_year, programme_key, relative_loss,
        cot_idwmeanpertow, cots_outbreak_probability, cots_label,
        disturbance_text, pre_cover, post_cover
    )
write_csv(
    cots_high_pressure_low_loss,
    file.path(out_dir, 'cots_high_pressure_low_loss_audit.csv')
)

debbie_annual <- annual |>
    filter(event_year == 2017, cyc_maxHrs4mw > 20) |>
    select(
        ReefID, ReefName, programme_key, lon, lat, relative_loss,
        cover_change_pp, cyc_maxHrs4mw, storm_name, disturbance_text
    )
write_csv(debbie_annual, file.path(out_dir, 'debbie_wave_loss_audit.csv'))

comparison_plot_data <- comparison |>
    select(candidate, scheme, rmse, severe_rmse, predictive_r2,
           false_extreme_rate) |>
    pivot_longer(
        c(rmse, severe_rmse, predictive_r2, false_extreme_rate),
        names_to = 'metric', values_to = 'value'
    ) |>
    mutate(
        model = recode(
            candidate,
            selected_current = 'Current selected',
            cause_aware_thermal_only = 'Cause-filtered thermal',
            cause_aware_competing_hazards = 'Ungated competing hazards',
            cause_aware_exposure_gated_competing_hazards =
                'Exposure-gated competing hazards'
        ),
        metric = recode(
            metric, rmse = 'RMSE', severe_rmse = 'Severe-event RMSE',
            predictive_r2 = 'Predictive R-squared',
            false_extreme_rate = 'False-extreme rate'
        )
    )
comparison_plot <- ggplot(comparison_plot_data, aes(model, value, fill = model)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = sprintf('%.3f', value)), vjust = -0.25, size = 3) +
    facet_grid(metric ~ scheme, scales = 'free_y') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.17))) +
    labs(
        title = 'Cause-aware competing hazards: held-out validation',
        x = NULL, y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_figure_bundle(
    comparison_plot, 'Fig-INLA-11_cause_aware_validation',
    comparison_plot_data,
    'Held-out performance of the current selected model, a cause-filtered thermal/freshwater INLA component, and its combination with independently trained COTS and cyclone hazards.',
    'Tests whether cause-labelled annual losses improve transfer without flattening the thermal response.',
    'Cause labels are incomplete and annual cover changes may combine multiple processes; event-overlapping intervals are excluded from each leave-event-out cause fit.',
    'cause_aware_competing_hazards', 'INLA composite', 'model_comparison',
    'operational_candidate_test', root, TRUE,
    code_source = 'scripts/fit_cause_aware_competing_hazards.R',
    width = 10, height = 8
)

cots_plot_data <- annual |>
    filter(cots_label, event_year %in% 2016:2021) |>
    mutate(
        focal = str_detect(ReefName, regex(
            'Gannett|Chinaman|Taylor|Rib Reef', ignore_case = TRUE
        )),
        label = if_else(focal, paste0(ReefName, ' ', event_year), NA_character_)
    )
cots_plot <- ggplot(
    cots_plot_data,
    aes(cot_idwmeanpertow, relative_loss, colour = factor(event_year))
) +
    geom_vline(xintercept = 0.22, linetype = 2, colour = 'grey45') +
    geom_point(aes(size = focal), alpha = 0.65) +
    ggrepel::geom_text_repel(
        data = cots_plot_data |> filter(focal),
        aes(label = label), size = 3, show.legend = FALSE, max.overlaps = Inf
    ) +
    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
    scale_size_manual(values = c('FALSE' = 1.7, 'TRUE' = 3.8),
                      guide = 'none') +
    labs(
        title = 'Cause-labelled COTS losses in the annual transition record',
        subtitle = 'Dashed line: 0.22 COTS/tow outbreak threshold',
        x = 'RRN COTS density per tow (log scale)',
        y = 'Observed relative coral-cover loss', colour = 'Year'
    ) +
    theme_bw(base_size = 11)
save_figure_bundle(
    cots_plot, 'Fig-DATA-02_cots_cause_evidence', cots_plot_data,
    'Relative coral-cover loss against RRN COTS density for cause-labelled annual transitions from 2016 to 2021, highlighting Gannett Cay, Chinaman, Taylor and Rib reefs.',
    'Shows the outbreak observations used to train the independent COTS hazard and the high-pressure/low-loss cases that constrain it.',
    'Cover-transition labels are observational and can contain multiple causes; the 0.22 line is the outbreak threshold, not a deterministic mortality threshold.',
    'cause_aware_cots_evidence', 'data', 'cause_evidence',
    'operational_candidate_test', root, TRUE,
    code_source = 'scripts/fit_cause_aware_competing_hazards.R'
)

debbie_track <- read_csv(
    file.path(root, 'data', 'processed', 'bom_cyclone_tracks_gbr.csv'),
    show_col_types = FALSE
) |>
    filter(NAME == 'Debbie', event_year == 2017) |>
    arrange(track_time_utc)

current_2017 <- baseline_metric_rows |>
    filter(scheme == 'leave_one_event_out', event_year == 2017) |>
    left_join(
        data |> distinct(
            source_observation_id = as.character(source_observation_id),
            lon, lat, cyc_interval_maxHrs4mw
        ),
        by = 'source_observation_id'
    )
candidate_2017 <- gated_predictions |>
    filter(scheme == 'leave_one_event_out', event_year == 2017)

map_limits <- list(x = c(146, 153), y = c(-24, -16))
debbie_map <- ggplot() +
    geom_path(
        data = debbie_track, aes(LON, LAT), colour = 'black',
        linewidth = 0.9, arrow = arrow(length = unit(0.12, 'cm'))
    ) +
    geom_point(
        data = debbie_annual,
        aes(lon, lat, size = cyc_maxHrs4mw, colour = relative_loss),
        alpha = 0.8
    ) +
    scale_colour_gradient2(
        low = '#2166AC', mid = 'white', high = '#B2182B', midpoint = 0.5,
        limits = c(0, 1), name = 'Relative loss'
    ) +
    scale_size_continuous(name = 'Wave hours >4 m', range = c(2.5, 7)) +
    coord_cartesian(xlim = map_limits$x, ylim = map_limits$y) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs(
        title = '2017 annual losses with >20 damaging-wave hours',
        x = 'Longitude', y = 'Latitude'
    ) + theme_bw(base_size = 10)

residual_plot_data <- bind_rows(
    current_2017 |>
        transmute(model = 'Current selected', ReefName, lon, lat,
                  residual, wave_hours = cyc_interval_maxHrs4mw),
    candidate_2017 |>
        transmute(model = 'Exposure-gated hazards', ReefName, lon, lat,
                  residual, wave_hours = cyclone_wave_hours)
)
residual_map <- ggplot(residual_plot_data, aes(lon, lat)) +
    geom_path(
        data = debbie_track, aes(LON, LAT), inherit.aes = FALSE,
        colour = 'black', linewidth = 0.8
    ) +
    geom_point(
        aes(colour = residual, size = pmax(abs(residual), 0.03)),
        alpha = 0.8
    ) +
    facet_wrap(~ model) +
    scale_colour_gradient2(
        low = '#2166AC', mid = 'white', high = '#B2182B', midpoint = 0,
        name = 'Observed - predicted'
    ) +
    scale_size_continuous(guide = 'none', range = c(2, 7)) +
    coord_cartesian(xlim = map_limits$x, ylim = map_limits$y) +
    scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs(
        title = '2017 event-held-out residuals and Cyclone Debbie track',
        x = 'Longitude', y = 'Latitude'
    ) + theme_bw(base_size = 10)

debbie_combined <- debbie_map / residual_map +
    plot_annotation(
        title = 'Cyclone Debbie evidence and the 2017 mortality residual field'
    )
save_figure_bundle(
    debbie_combined, 'Fig-INLA-12_debbie_residual_audit',
    bind_rows(
        debbie_annual |> mutate(panel = 'annual_wave_loss'),
        residual_plot_data |> mutate(panel = 'heldout_residual')
    ),
    'Cyclone Debbie track, 2017 annual coral losses at reefs with more than 20 hours of waves above 4 m, and held-out residuals before and after the competing cyclone hazard.',
    'Tests whether the spatial underprediction around Debbie is supported by independent wave-exposure and cover-loss evidence, with Penrith as the bleaching-mortality positive control.',
    'The annual transition and bleaching-mortality panels use different response records; cyclone labels are observational and the track layer remains provisional.',
    'cause_aware_debbie_audit', 'INLA composite', 'spatial_residual_audit',
    'operational_candidate_test', root, TRUE,
    code_source = 'scripts/fit_cause_aware_competing_hazards.R',
    width = 11, height = 10
)

write_figure_readme(root)
message('Wrote cause-aware competing-hazard candidate outputs')
