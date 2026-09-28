# Expanded multi-event BRT screen for reef-level eReefs surface salinity.
# Direct salinity is never imputed. The bleaching-event analysis preserves the
# current mortality outcome; an explicitly separate annual-transition analysis
# brings the 2010-11 wet season into scope without changing the selected model.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(ggplot2)
    library(readr)
    library(tidyr)
})

set.seed(20260928L)
output_dir <- 'output/surface_salinity_mortality'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

salinity <- read_csv(
    'data/processed/ereefs_surface_salinity_reef_event.csv',
    show_col_types = FALSE
) |>
    transmute(
        ReefID = toupper(trimws(ReefID)), event_year,
        sss_min = sss_min_area_mean_psu,
        hours_below30 = hours_below30_area_mean,
        hours_below26 = hours_below26_area_mean,
        log_hours_below30 = log1p_hours_below30_area_mean,
        log_hours_below26 = log1p_hours_below26_area_mean
    )
coverage <- read_csv(
    'data/processed/ereefs_surface_salinity_reef_grid_coverage.csv',
    show_col_types = FALSE
) |>
    mutate(ReefID = toupper(trimws(ReefID)))
if (anyDuplicated(salinity[c('ReefID', 'event_year')])) {
    stop('Salinity input is not unique by ReefID-event_year')
}

collapse_text <- function(x) {
    x <- unique(na.omit(x[x != '']))
    paste(x, collapse = ' | ')
}

model_rows <- read_csv(
    'output/explanatory_event_dhw/event_dhw_brt_data.csv',
    show_col_types = FALSE
) |>
    mutate(
        ReefID = toupper(trimws(ReefID)),
        across(
            starts_with('disturbance_has_'),
            ~ as.numeric(coalesce(.x, FALSE))
        )
    ) |>
    left_join(salinity, by = c('ReefID', 'event_year')) |>
    left_join(
        coverage |> select(ReefID, valid_area_fraction, coverage_class),
        by = 'ReefID', relationship = 'many-to-one'
    ) |>
    filter(
        programme_key %in% c('ltmp', 'manta', 'mmp'),
        is.finite(mortality_prop),
        is.finite(sss_min), is.finite(hours_below30),
        is.finite(hours_below26), valid_area_fraction >= 0.5
    )

bleaching_data <- model_rows |>
    group_by(ReefID, event_year, programme_key, depth, joint_reef_fold) |>
    summarise(
        ReefName = first(ReefName),
        mortality = mean(mortality_prop),
        source_rows = n(),
        thermal_dhw = median(ann_maxdhw, na.rm = TRUE),
        cots_hazard = median(cots_hazard_weight, na.rm = TRUE),
        cyclone_wave_hours = median(cyc_interval_maxHrs4mw, na.rm = TRUE),
        sss_min = median(sss_min),
        hours_below30 = median(hours_below30),
        hours_below26 = median(hours_below26),
        log_hours_below30 = median(log_hours_below30),
        log_hours_below26 = median(log_hours_below26),
        valid_area_fraction = first(valid_area_fraction),
        bleaching_label = max(disturbance_has_bleaching),
        flood_label = max(disturbance_has_flood),
        cyclone_label = max(disturbance_has_cyclone),
        cots_label = max(disturbance_has_cots),
        disturbance_text = collapse_text(disturbance_text),
        .groups = 'drop'
    ) |>
    mutate(
        analysis_context = 'bleaching_mortality',
        reef_fold = as.integer(joint_reef_fold),
        pre_cover = 0,
        observation_id = paste(
            analysis_context, ReefID, event_year, programme_key, depth,
            sep = '__'
        )
    )

existing_fold_map <- model_rows |>
    distinct(ReefID, joint_reef_fold) |>
    filter(is.finite(joint_reef_fold)) |>
    group_by(ReefID) |>
    summarise(reef_fold = first(as.integer(joint_reef_fold)), .groups = 'drop')
stable_fold <- function(id) {
    values <- utf8ToInt(id)
    as.integer(sum(values * seq_along(values)) %% 5L + 1L)
}

annual_rows <- read_csv(
    'data/processed/annual_coral_transitions.csv', show_col_types = FALSE
) |>
    mutate(
        ReefID = toupper(trimws(ReefID)),
        across(
            starts_with('disturbance_has_'),
            ~ as.numeric(coalesce(.x, FALSE))
        )
    ) |>
    filter(event_year >= 2011L, event_year <= 2024L) |>
    left_join(salinity, by = c('ReefID', 'event_year')) |>
    left_join(
        coverage |> select(ReefID, valid_area_fraction, coverage_class),
        by = 'ReefID', relationship = 'many-to-one'
    ) |>
    filter(
        programme_key %in% c('ltmp', 'manta', 'mmp'),
        is.finite(pre_cover), pre_cover >= 0.02,
        is.finite(post_cover),
        is.finite(sss_min), is.finite(hours_below30),
        is.finite(hours_below26), valid_area_fraction >= 0.5
    ) |>
    mutate(
        relative_loss = pmin(pmax((pre_cover - post_cover) / pre_cover, 0), 1)
    )

annual_data <- annual_rows |>
    group_by(ReefID, event_year, programme_key, depth) |>
    summarise(
        ReefName = first(ReefName),
        mortality = mean(relative_loss),
        source_rows = n(),
        pre_cover = mean(pre_cover),
        thermal_dhw = median(sst_maxdhw, na.rm = TRUE),
        cots_hazard = median(log1p(pmax(cot_idwmeanpertow, 0)), na.rm = TRUE),
        cyclone_wave_hours = median(cyc_maxHrs4mw, na.rm = TRUE),
        sss_min = median(sss_min),
        hours_below30 = median(hours_below30),
        hours_below26 = median(hours_below26),
        log_hours_below30 = median(log_hours_below30),
        log_hours_below26 = median(log_hours_below26),
        valid_area_fraction = first(valid_area_fraction),
        bleaching_label = max(disturbance_has_bleaching),
        flood_label = max(disturbance_has_flood),
        cyclone_label = max(disturbance_has_cyclone),
        cots_label = max(disturbance_has_cots),
        disturbance_text = collapse_text(disturbance_text),
        .groups = 'drop'
    ) |>
    left_join(existing_fold_map, by = 'ReefID') |>
    mutate(
        reef_fold = coalesce(
            reef_fold,
            vapply(ReefID, stable_fold, integer(1))
        ),
        analysis_context = 'annual_transition_sensitivity',
        observation_id = paste(
            analysis_context, ReefID, event_year, programme_key, depth,
            sep = '__'
        )
    )

make_analysis_fields <- function(x) {
    x |>
        mutate(
            depth = as.numeric(depth),
            is_manta = as.numeric(programme_key == 'manta'),
            is_mmp = as.numeric(programme_key == 'mmp'),
            no_recorded_cause = as.numeric(
                bleaching_label == 0 & flood_label == 0 &
                    cyclone_label == 0 & cots_label == 0
            ),
            stratum = case_when(
                programme_key == 'mmp' & depth <= 3 ~ 'MMP 2 m',
                programme_key == 'mmp' ~ 'MMP 5 m',
                programme_key == 'manta' ~ 'Manta 9 m',
                TRUE ~ 'LTMP 9 m'
            )
        )
}
bleaching_data <- make_analysis_fields(bleaching_data)
annual_data <- make_analysis_fields(annual_data)
write_csv(
    bleaching_data,
    file.path(output_dir, 'expanded_bleaching_salinity_data.csv'), na = ''
)
write_csv(
    annual_data,
    file.path(output_dir, 'expanded_annual_transition_salinity_data.csv'),
    na = ''
)

support <- bind_rows(bleaching_data, annual_data) |>
    group_by(analysis_context, event_year, programme_key, depth) |>
    summarise(
        rows = n(), reefs = n_distinct(ReefID),
        mean_mortality = mean(mortality),
        severe_rows = sum(mortality >= 0.2),
        flood_rows = sum(flood_label > 0),
        no_recorded_cause_rows = sum(no_recorded_cause > 0),
        minimum_sss = min(sss_min),
        maximum_hours_below30 = max(hours_below30),
        maximum_hours_below26 = max(hours_below26),
        .groups = 'drop'
    )
write_csv(support, file.path(output_dir, 'expanded_salinity_support.csv'))

correlations <- bind_rows(bleaching_data, annual_data) |>
    group_by(analysis_context) |>
    summarise(
        rows = n(), reefs = n_distinct(ReefID), events = n_distinct(event_year),
        spearman_sss_hours30 = cor(
            sss_min, hours_below30, method = 'spearman'
        ),
        spearman_sss_hours26 = cor(
            sss_min, hours_below26, method = 'spearman'
        ),
        spearman_hours30_hours26 = cor(
            hours_below30, hours_below26, method = 'spearman'
        ),
        hours30_positive_fraction = mean(hours_below30 > 0),
        hours26_positive_fraction = mean(hours_below26 > 0),
        .groups = 'drop'
    )
write_csv(correlations, file.path(output_dir, 'expanded_salinity_correlations.csv'))

survey_terms <- c('depth', 'is_manta', 'is_mmp')
exposure_terms <- c(
    survey_terms, 'thermal_dhw', 'cots_hazard', 'cyclone_wave_hours',
    'pre_cover'
)
label_terms <- c(
    'bleaching_label', 'flood_label', 'cyclone_label', 'cots_label'
)
candidates <- list(
    exposures = exposure_terms,
    exposures_sss = c(exposure_terms, 'sss_min'),
    exposures_h30 = c(exposure_terms, 'log_hours_below30'),
    exposures_h26 = c(exposure_terms, 'log_hours_below26'),
    exposures_sss_h26 = c(
        exposure_terms, 'sss_min', 'log_hours_below26'
    ),
    exposures_all = c(
        exposure_terms, 'sss_min', 'log_hours_below30',
        'log_hours_below26'
    ),
    labels = c(exposure_terms, label_terms),
    labels_sss = c(exposure_terms, label_terms, 'sss_min'),
    labels_h30 = c(exposure_terms, label_terms, 'log_hours_below30'),
    labels_h26 = c(exposure_terms, label_terms, 'log_hours_below26'),
    labels_all = c(
        exposure_terms, label_terms, 'sss_min', 'log_hours_below30',
        'log_hours_below26'
    )
)

prepare_fold <- function(analysis, assessment, predictors) {
    active <- character()
    rules <- tibble(predictor = predictors, fill = NA_real_, varying = FALSE)
    for (i in seq_along(predictors)) {
        predictor <- predictors[[i]]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        fill <- if (length(observed)) median(observed) else 0
        analysis[[predictor]][!is.finite(analysis[[predictor]])] <- fill
        assessment[[predictor]][!is.finite(assessment[[predictor]])] <- fill
        varying <- n_distinct(analysis[[predictor]]) > 1L
        rules$fill[[i]] <- fill
        rules$varying[[i]] <- varying
        if (varying) active <- c(active, predictor)
    }
    list(
        analysis = analysis, assessment = assessment,
        active = active, rules = rules
    )
}

fit_brt <- function(data, predictors, seed, trees = 700L) {
    active <- predictors[vapply(
        predictors,
        function(x) n_distinct(data[[x]]) > 1L,
        logical(1)
    )]
    set.seed(seed)
    gbm(
        reformulate(active, response = 'mortality'),
        data = data, distribution = 'gaussian',
        n.trees = trees, interaction.depth = 2L,
        shrinkage = 0.02, n.minobsinnode = 4L,
        bag.fraction = 0.75, train.fraction = 1,
        keep.data = FALSE, verbose = FALSE
    )
}
predict_brt <- function(model, data) {
    pmin(pmax(as.numeric(predict(
        model, data, n.trees = model$n.trees
    )), 0), 1)
}

run_validation <- function(data, context_index) {
    predictions <- tibble()
    preprocessing <- tibble()
    schemes <- list(
        leave_one_event_out = 'event_year',
        reef_blocked_5fold = 'reef_fold'
    )
    for (scheme in names(schemes)) {
        fold_column <- schemes[[scheme]]
        folds <- sort(unique(data[[fold_column]]))
        for (candidate in names(candidates)) {
            for (fold_index in seq_along(folds)) {
                fold <- folds[[fold_index]]
                analysis <- data[data[[fold_column]] != fold, ]
                assessment <- data[data[[fold_column]] == fold, ]
                prepared <- prepare_fold(
                    analysis, assessment, candidates[[candidate]]
                )
                model <- fit_brt(
                    prepared$analysis, prepared$active,
                    20260928L + context_index * 10000L +
                        match(candidate, names(candidates)) * 100L +
                        fold_index
                )
                predictions <- bind_rows(
                    predictions,
                    prepared$assessment |>
                        transmute(
                            observation_id, ReefID, ReefName, event_year,
                            programme_key, depth, stratum, mortality,
                            flood_label, bleaching_label, cyclone_label,
                            cots_label, no_recorded_cause,
                            analysis_context, validation_scheme = scheme,
                            fold = as.character(fold), candidate,
                            predicted_mortality = predict_brt(
                                model, prepared$assessment
                            )
                        )
                )
                preprocessing <- bind_rows(
                    preprocessing,
                    prepared$rules |>
                        mutate(
                            analysis_context = first(data$analysis_context),
                            validation_scheme = scheme,
                            candidate, fold = as.character(fold)
                        )
                )
            }
        }
    }
    list(predictions = predictions, preprocessing = preprocessing)
}

validated <- list(
    run_validation(bleaching_data, 1L),
    run_validation(annual_data, 2L)
)
predictions <- bind_rows(lapply(validated, `[[`, 'predictions'))
preprocessing <- bind_rows(lapply(validated, `[[`, 'preprocessing'))
write_csv(
    predictions,
    file.path(output_dir, 'expanded_salinity_brt_predictions.csv'), na = ''
)
write_csv(
    preprocessing,
    file.path(output_dir, 'expanded_salinity_brt_preprocessing.csv'), na = ''
)

metric_summary <- function(x, groups) {
    x |>
        group_by(across(all_of(groups))) |>
        summarise(
            rows = n(), reefs = n_distinct(ReefID),
            events = n_distinct(event_year),
            rmse = sqrt(mean((mortality - predicted_mortality)^2)),
            mae = mean(abs(mortality - predicted_mortality)),
            predictive_r2 = 1 - sum(
                (mortality - predicted_mortality)^2
            ) / sum((mortality - mean(mortality))^2),
            bias = mean(predicted_mortality - mortality),
            severe_rows = sum(mortality >= 0.2),
            severe_rmse = if_else(
                severe_rows > 0,
                sqrt(mean(
                    (mortality[mortality >= 0.2] -
                        predicted_mortality[mortality >= 0.2])^2
                )),
                NA_real_
            ),
            false_extreme_rate = mean(
                predicted_mortality >= 0.5 & mortality < 0.2
            ),
            .groups = 'drop'
        )
}
metrics <- metric_summary(
    predictions,
    c('analysis_context', 'validation_scheme', 'candidate')
)
event_metrics <- metric_summary(
    predictions,
    c(
        'analysis_context', 'validation_scheme', 'candidate', 'event_year'
    )
)
write_csv(metrics, file.path(output_dir, 'expanded_salinity_brt_metrics.csv'))
write_csv(
    event_metrics,
    file.path(output_dir, 'expanded_salinity_brt_event_metrics.csv')
)

comparisons <- tribble(
    ~candidate, ~baseline,
    'exposures_sss', 'exposures',
    'exposures_h30', 'exposures',
    'exposures_h26', 'exposures',
    'exposures_sss_h26', 'exposures',
    'exposures_all', 'exposures',
    'labels_sss', 'labels',
    'labels_h30', 'labels',
    'labels_h26', 'labels',
    'labels_all', 'labels'
)

delta_rows <- list()
index <- 0L
for (context in unique(predictions$analysis_context)) {
    for (scheme in unique(predictions$validation_scheme)) {
        subset <- predictions |>
            filter(
                analysis_context == context,
                validation_scheme == scheme
            )
        for (i in seq_len(nrow(comparisons))) {
            index <- index + 1L
            candidate_name <- comparisons$candidate[[i]]
            baseline_name <- comparisons$baseline[[i]]
            paired <- subset |>
                filter(candidate == candidate_name) |>
                select(
                    observation_id, ReefID, event_year, mortality,
                    candidate_prediction = predicted_mortality
                ) |>
                inner_join(
                    subset |>
                        filter(candidate == baseline_name) |>
                        select(
                            observation_id,
                            baseline_prediction = predicted_mortality
                        ),
                    by = 'observation_id', relationship = 'one-to-one'
                )
            reefs <- unique(paired$ReefID)
            set.seed(20261028L + index)
            bootstrap_delta <- replicate(1000L, {
                sampled <- sample(reefs, length(reefs), replace = TRUE)
                row_index <- unlist(lapply(
                    sampled, function(id) which(paired$ReefID == id)
                ))
                sqrt(mean(
                    (paired$mortality[row_index] -
                        paired$candidate_prediction[row_index])^2
                )) - sqrt(mean(
                    (paired$mortality[row_index] -
                        paired$baseline_prediction[row_index])^2
                ))
            })
            delta_rows[[index]] <- tibble(
                analysis_context = context,
                validation_scheme = scheme,
                candidate = candidate_name, baseline = baseline_name,
                rows = nrow(paired), reefs = length(reefs),
                delta_rmse = sqrt(mean(
                    (paired$mortality - paired$candidate_prediction)^2
                )) - sqrt(mean(
                    (paired$mortality - paired$baseline_prediction)^2
                )),
                delta_rmse_q025 = quantile(bootstrap_delta, 0.025),
                delta_rmse_q975 = quantile(bootstrap_delta, 0.975),
                bootstrap_probability_improved = mean(bootstrap_delta < 0)
            )
        }
    }
}
deltas <- bind_rows(delta_rows)
write_csv(deltas, file.path(output_dir, 'expanded_salinity_brt_deltas.csv'))

full_data <- list(
    bleaching_mortality = bleaching_data,
    annual_transition_sensitivity = annual_data
)
metric_candidate <- c(
    sss_min = 'exposures_sss',
    log_hours_below30 = 'exposures_h30',
    log_hours_below26 = 'exposures_h26'
)
metric_label <- c(
    sss_min = 'Minimum surface salinity (PSU)',
    log_hours_below30 = 'log(1 + hours below 30 PSU)',
    log_hours_below26 = 'log(1 + hours below 26 PSU)'
)
strata <- tibble(
    stratum = c('LTMP 9 m', 'Manta 9 m', 'MMP 2 m', 'MMP 5 m'),
    depth = c(9, 9, 2, 5),
    is_manta = c(0, 1, 0, 0),
    is_mmp = c(0, 0, 1, 1)
)
importance <- tibble()
curves <- tibble()
curve_bootstrap <- tibble()

for (context_index in seq_along(full_data)) {
    context <- names(full_data)[[context_index]]
    data <- full_data[[context_index]]
    for (metric in names(metric_candidate)) {
        candidate <- metric_candidate[[metric]]
        predictors <- candidates[[candidate]]
        prepared <- prepare_fold(data, data, predictors)
        model <- fit_brt(
            prepared$analysis, prepared$active,
            20262028L + context_index * 100L + match(metric, names(metric_candidate))
        )
        importance <- bind_rows(
            importance,
            summary(model, plotit = FALSE) |>
                as_tibble() |>
                rename(predictor = var, relative_influence = rel.inf) |>
                mutate(
                    analysis_context = context,
                    metric, candidate,
                    rank = row_number()
                )
        )
        values <- data[[metric]]
        grid <- unique(as.numeric(quantile(
            values, probs = seq(0, 1, length.out = 31), na.rm = TRUE
        )))
        curve_for_model <- function(fitted_model, standardization_data) {
            bind_rows(lapply(seq_len(nrow(strata)), function(s) {
                bind_rows(lapply(grid, function(value) {
                    new_data <- standardization_data
                    new_data[[metric]] <- value
                    new_data$depth <- strata$depth[[s]]
                    new_data$is_manta <- strata$is_manta[[s]]
                    new_data$is_mmp <- strata$is_mmp[[s]]
                    tibble(
                        stratum = strata$stratum[[s]], metric_value = value,
                        predicted_mortality = mean(predict_brt(
                            fitted_model, new_data
                        ))
                    )
                }))
            }))
        }
        point <- curve_for_model(model, prepared$assessment) |>
            mutate(analysis_context = context, metric, candidate)
        curves <- bind_rows(curves, point)

        reefs <- unique(data$ReefID)
        for (b in seq_len(30L)) {
            set.seed(
                20263028L + context_index * 10000L +
                    match(metric, names(metric_candidate)) * 100L + b
            )
            sampled <- sample(reefs, length(reefs), replace = TRUE)
            rows <- unlist(lapply(
                sampled, function(id) which(data$ReefID == id)
            ))
            boot_data <- data[rows, ]
            boot_prepared <- prepare_fold(boot_data, data, predictors)
            boot_model <- tryCatch(
                fit_brt(
                    boot_prepared$analysis, boot_prepared$active,
                    20264028L + b, trees = 500L
                ),
                error = function(e) NULL
            )
            if (is.null(boot_model)) next
            curve_bootstrap <- bind_rows(
                curve_bootstrap,
                curve_for_model(boot_model, boot_prepared$assessment) |>
                    mutate(
                        analysis_context = context, metric, candidate,
                        bootstrap = b
                    )
            )
        }
    }
}
curve_intervals <- curve_bootstrap |>
    group_by(
        analysis_context, metric, candidate, stratum, metric_value
    ) |>
    summarise(
        lower = quantile(predicted_mortality, 0.025),
        upper = quantile(predicted_mortality, 0.975),
        successful_bootstraps = n_distinct(bootstrap),
        .groups = 'drop'
    )
curves <- curves |>
    left_join(
        curve_intervals,
        by = c(
            'analysis_context', 'metric', 'candidate', 'stratum',
            'metric_value'
        )
    ) |>
    mutate(metric_label = recode(metric, !!!metric_label))
write_csv(
    importance,
    file.path(output_dir, 'expanded_salinity_brt_importance.csv')
)
write_csv(
    curves,
    file.path(output_dir, 'expanded_salinity_brt_partial_dependence.csv')
)

curve_plot <- ggplot(
    curves,
    aes(metric_value, predicted_mortality, colour = stratum, fill = stratum)
) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.08, colour = NA) +
    geom_line(linewidth = 0.75) +
    facet_grid(analysis_context ~ metric_label, scales = 'free_x') +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
        x = NULL, y = 'Partial expected relative loss',
        colour = 'Programme/depth', fill = 'Programme/depth',
        title = 'Expanded multi-event salinity BRT relationships',
        subtitle = paste(
            'Bands are 30 reef-cluster screening bootstraps;',
            'each salinity metric is fit separately after heat, COTS and cyclone'
        )
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'bottom')
ggsave(
    file.path(output_dir, 'expanded_salinity_brt_partial_dependence.png'),
    curve_plot, width = 14, height = 8, dpi = 180
)

influential <- bind_rows(bleaching_data, annual_data) |>
    group_by(analysis_context) |>
    arrange(desc(hours_below26), desc(hours_below30), .by_group = TRUE) |>
    slice_head(n = 25L) |>
    ungroup() |>
    select(
        analysis_context, ReefID, ReefName, event_year, programme_key,
        depth, mortality, sss_min, hours_below30, hours_below26,
        thermal_dhw, cots_hazard, cyclone_wave_hours,
        bleaching_label, flood_label, cyclone_label, cots_label,
        disturbance_text
    )
write_csv(
    influential,
    file.path(output_dir, 'expanded_salinity_extreme_rows.csv'), na = ''
)

cat('Expanded salinity BRT screen complete.\n')
print(metrics |> filter(validation_scheme == 'leave_one_event_out') |>
    arrange(analysis_context, rmse))
print(deltas |> filter(validation_scheme == 'leave_one_event_out'))
