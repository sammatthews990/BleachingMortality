# Multi-programme BRT diagnostic for salinity after major mortality causes.
#
# The primary candidate uses only survey structure, local-first DHW, the
# selected-model COTS hazard index, cyclone wave exposure and the two supplied
# surface-salinity metrics. A separate retrospective sensitivity adds field
# bleaching/flood/cyclone/COTS labels. Those labels are explanatory only and
# must not be treated as operational forecast inputs.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(ggplot2)
    library(readr)
    library(tidyr)
})

set.seed(20260916L)

salinity_file <- 'data/processed/ereefs_surface_salinity_reef_event.csv'
model_rows_file <- 'output/explanatory_event_dhw/event_dhw_brt_data.csv'
output_dir <- 'output/surface_salinity_mortality'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!all(file.exists(c(salinity_file, model_rows_file)))) {
    stop('Missing salinity or explanatory-event-DHW input')
}

salinity <- read_csv(salinity_file, show_col_types = FALSE) |>
    select(
        ReefID, event_year, sss_min_area_mean_psu,
        hours_below30_area_mean
    )
if (anyDuplicated(salinity[c('ReefID', 'event_year')])) {
    stop('Salinity input is not unique by ReefID-event_year')
}

rows <- read_csv(model_rows_file, show_col_types = FALSE) |>
    mutate(
        ReefID = toupper(trimws(ReefID)),
        disturbance_has_bleaching = coalesce(
            disturbance_has_bleaching, FALSE
        ),
        disturbance_has_flood = coalesce(disturbance_has_flood, FALSE),
        disturbance_has_cyclone = coalesce(
            disturbance_has_cyclone, FALSE
        ),
        disturbance_has_cots = coalesce(disturbance_has_cots, FALSE)
    ) |>
    left_join(
        salinity,
        by = c('ReefID', 'event_year'),
        relationship = 'many-to-one'
    ) |>
    filter(
        event_year == 2024L,
        is.finite(mortality_prop),
        is.finite(sss_min_area_mean_psu),
        is.finite(hours_below30_area_mean),
        programme_key %in% c('mmp', 'manta', 'ltmp')
    )

# Equal-weight reef-programme-depth outcomes prevent duplicated transitions
# from giving a reef extra model weight. Exposure values should be constant
# within these groups; medians make that contract robust to tiny numeric noise.
data <- rows |>
    group_by(ReefID, programme_key, depth, joint_reef_fold) |>
    summarise(
        ReefName = first(ReefName),
        mortality = mean(mortality_prop),
        source_rows = n(),
        local_first_dhw = median(ann_maxdhw, na.rm = TRUE),
        cots_hazard = median(cots_hazard_weight, na.rm = TRUE),
        cyclone_wave_hours = median(
            cyc_interval_maxHrs4mw, na.rm = TRUE
        ),
        sss_min = median(sss_min_area_mean_psu, na.rm = TRUE),
        hours_below30 = median(hours_below30_area_mean, na.rm = TRUE),
        bleaching_label = any(disturbance_has_bleaching),
        flood_label = any(disturbance_has_flood),
        cyclone_label = any(disturbance_has_cyclone),
        cots_label = any(disturbance_has_cots),
        .groups = 'drop'
    ) |>
    mutate(
        observation_id = paste(ReefID, programme_key, depth, sep = '__'),
        is_manta = as.numeric(programme_key == 'manta'),
        is_mmp = as.numeric(programme_key == 'mmp'),
        depth = as.numeric(depth),
        across(
            c(
                bleaching_label, flood_label, cyclone_label, cots_label
            ),
            as.numeric
        ),
        no_recorded_cause = as.numeric(
            bleaching_label == 0 & flood_label == 0 &
                cyclone_label == 0 & cots_label == 0
        )
    ) |>
    arrange(programme_key, ReefID, depth)

if (anyDuplicated(data$observation_id)) {
    stop('Cause-adjusted BRT data are not unique by observation_id')
}
if (n_distinct(data$joint_reef_fold) < 5L) {
    stop('Expected five existing reef-blocked folds')
}

write_csv(
    data,
    file.path(output_dir, 'cause_adjusted_salinity_brt_data.csv'),
    na = ''
)

support <- data |>
    group_by(programme_key, depth) |>
    summarise(
        rows = n(), reefs = n_distinct(ReefID),
        mean_mortality = mean(mortality),
        severe_rows = sum(mortality >= 0.2),
        bleaching_rows = sum(bleaching_label),
        flood_rows = sum(flood_label),
        cyclone_rows = sum(cyclone_label),
        cots_rows = sum(cots_label),
        sss_minimum = min(sss_min), sss_maximum = max(sss_min),
        maximum_hours_below30 = max(hours_below30),
        maximum_local_first_dhw = max(local_first_dhw),
        maximum_cyclone_wave_hours = max(cyclone_wave_hours),
        .groups = 'drop'
    )
write_csv(
    support,
    file.path(output_dir, 'cause_adjusted_salinity_brt_support.csv'),
    na = ''
)

survey_predictors <- c('depth', 'is_manta', 'is_mmp')
cause_predictors <- c(
    survey_predictors, 'local_first_dhw', 'cots_hazard',
    'cyclone_wave_hours'
)
salinity_predictors <- c('sss_min', 'hours_below30')
label_predictors <- c(
    'bleaching_label', 'flood_label', 'cyclone_label', 'cots_label'
)
candidate_predictors <- list(
    survey_only = survey_predictors,
    cause_exposures = cause_predictors,
    cause_plus_sss = c(cause_predictors, 'sss_min'),
    cause_plus_hours = c(cause_predictors, 'hours_below30'),
    cause_plus_both = c(cause_predictors, salinity_predictors),
    cause_labels = c(cause_predictors, label_predictors),
    labels_plus_sss = c(
        cause_predictors, label_predictors, 'sss_min'
    ),
    labels_plus_hours = c(
        cause_predictors, label_predictors, 'hours_below30'
    ),
    labels_plus_both = c(
        cause_predictors, label_predictors, salinity_predictors
    )
)

prepare_fold <- function(analysis, assessment, predictors) {
    rules <- tibble(
        predictor = predictors, median = NA_real_, zero_variance = FALSE
    )
    active <- character()
    for (i in seq_along(predictors)) {
        predictor <- predictors[[i]]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        fill <- if (length(observed)) median(observed) else 0
        analysis[[predictor]][!is.finite(analysis[[predictor]])] <- fill
        assessment[[predictor]][!is.finite(assessment[[predictor]])] <- fill
        zero <- length(unique(analysis[[predictor]])) < 2L
        rules$median[[i]] <- fill
        rules$zero_variance[[i]] <- zero
        if (!zero) active <- c(active, predictor)
    }
    list(
        analysis = analysis, assessment = assessment,
        active = active, rules = rules
    )
}

fit_brt <- function(data, predictors, seed, n_trees = 800L) {
    active <- predictors[vapply(
        predictors,
        function(predictor) length(unique(data[[predictor]])) > 1L,
        logical(1)
    )]
    if (length(active) == 0L) stop('No varying predictors in BRT sample')
    set.seed(seed)
    gbm(
        reformulate(active, response = 'mortality'),
        data = data,
        distribution = 'gaussian',
        n.trees = n_trees,
        interaction.depth = 2L,
        shrinkage = 0.02,
        n.minobsinnode = 3L,
        bag.fraction = 0.75,
        train.fraction = 1,
        keep.data = FALSE,
        verbose = FALSE
    )
}

predict_brt <- function(model, new_data) {
    pmin(pmax(as.numeric(predict(
        model, new_data, n.trees = model$n.trees
    )), 0), 1)
}

predictions <- tibble()
preprocessing <- tibble()
for (candidate in names(candidate_predictors)) {
    predictors <- candidate_predictors[[candidate]]
    for (fold in sort(unique(data$joint_reef_fold))) {
        analysis <- data |> filter(joint_reef_fold != fold)
        assessment <- data |> filter(joint_reef_fold == fold)
        prepared <- prepare_fold(analysis, assessment, predictors)
        model <- fit_brt(
            prepared$analysis, prepared$active,
            20260916L + 100L * match(
                candidate, names(candidate_predictors)
            ) + fold
        )
        predictions <- bind_rows(
            predictions,
            prepared$assessment |>
                transmute(
                    observation_id, ReefID, ReefName, programme_key,
                    depth, joint_reef_fold,
                    observed_mortality = mortality,
                    bleaching_label, flood_label, cyclone_label,
                    cots_label, no_recorded_cause,
                    candidate,
                    predicted_mortality = predict_brt(
                        model, prepared$assessment
                    )
                )
        )
        preprocessing <- bind_rows(
            preprocessing,
            prepared$rules |> mutate(candidate, fold)
        )
    }
}

summarise_metrics <- function(x, grouping = 'candidate') {
    x |>
        group_by(across(all_of(grouping))) |>
        summarise(
            rows = n(), reefs = n_distinct(ReefID),
            rmse = sqrt(mean(
                (observed_mortality - predicted_mortality)^2
            )),
            mae = mean(abs(observed_mortality - predicted_mortality)),
            predictive_r2 = 1 - sum(
                (observed_mortality - predicted_mortality)^2
            ) / sum(
                (observed_mortality - mean(observed_mortality))^2
            ),
            bias = mean(predicted_mortality - observed_mortality),
            severe_rows = sum(observed_mortality >= 0.2),
            severe_rmse = if (severe_rows > 0L) sqrt(mean(
                (observed_mortality[observed_mortality >= 0.2] -
                    predicted_mortality[observed_mortality >= 0.2])^2
            )) else NA_real_,
            false_extreme_rate = mean(
                predicted_mortality >= 0.5 & observed_mortality < 0.2
            ),
            .groups = 'drop'
        )
}

metrics <- summarise_metrics(predictions)
metrics_by_programme_depth <- summarise_metrics(
    predictions, c('candidate', 'programme_key', 'depth')
)

label_predictions <- bind_rows(lapply(
    c(
        'bleaching_label', 'flood_label', 'cyclone_label',
        'cots_label', 'no_recorded_cause'
    ),
    function(label_name) {
        predictions |>
            filter(.data[[label_name]] == 1) |>
            mutate(label_group = label_name)
    }
))
metrics_by_label <- summarise_metrics(
    label_predictions, c('candidate', 'label_group')
)

comparison_specs <- tribble(
    ~candidate, ~baseline,
    'cause_plus_sss', 'cause_exposures',
    'cause_plus_hours', 'cause_exposures',
    'cause_plus_both', 'cause_exposures',
    'labels_plus_sss', 'cause_labels',
    'labels_plus_hours', 'cause_labels',
    'labels_plus_both', 'cause_labels'
)

bootstrap_deltas <- function(prediction_rows, replicates = 2000L) {
    bind_rows(lapply(seq_len(nrow(comparison_specs)), function(i) {
        spec <- comparison_specs[i, ]
        baseline <- prediction_rows |>
            filter(candidate == spec$baseline) |>
            select(
                observation_id, ReefID, observed_mortality,
                baseline_prediction = predicted_mortality
            )
        x <- prediction_rows |>
            filter(candidate == spec$candidate) |>
            inner_join(
                baseline,
                by = c('observation_id', 'ReefID', 'observed_mortality'),
                relationship = 'many-to-one'
            )
        reefs <- unique(x$ReefID)
        set.seed(20261016L + i)
        delta <- replicate(replicates, {
            sampled <- sample(reefs, length(reefs), replace = TRUE)
            index <- unlist(lapply(
                sampled, function(id) which(x$ReefID == id)
            ))
            sqrt(mean(
                (x$observed_mortality[index] -
                    x$predicted_mortality[index])^2
            )) - sqrt(mean(
                (x$observed_mortality[index] -
                    x$baseline_prediction[index])^2
            ))
        })
        tibble(
            candidate = spec$candidate, baseline = spec$baseline,
            rows = nrow(x), reefs = length(reefs),
            delta_rmse = sqrt(mean(
                (x$observed_mortality - x$predicted_mortality)^2
            )) - sqrt(mean(
                (x$observed_mortality - x$baseline_prediction)^2
            )),
            delta_rmse_q025 = quantile(delta, 0.025),
            delta_rmse_q975 = quantile(delta, 0.975),
            bootstrap_probability_improved = mean(delta < 0)
        )
    }))
}

metric_deltas <- bootstrap_deltas(predictions)
write_csv(
    predictions,
    file.path(output_dir, 'cause_adjusted_salinity_brt_predictions.csv'),
    na = ''
)
write_csv(
    preprocessing,
    file.path(output_dir, 'cause_adjusted_salinity_brt_preprocessing.csv'),
    na = ''
)
write_csv(
    metrics,
    file.path(output_dir, 'cause_adjusted_salinity_brt_metrics.csv'),
    na = ''
)
write_csv(
    metrics_by_programme_depth,
    file.path(
        output_dir,
        'cause_adjusted_salinity_brt_metrics_by_programme_depth.csv'
    ),
    na = ''
)
write_csv(
    metrics_by_label,
    file.path(
        output_dir, 'cause_adjusted_salinity_brt_metrics_by_label.csv'
    ),
    na = ''
)
write_csv(
    metric_deltas,
    file.path(output_dir, 'cause_adjusted_salinity_brt_metric_deltas.csv'),
    na = ''
)

model_specs <- list(
    cause_exposures = candidate_predictors$cause_plus_both,
    cause_exposures_and_labels = candidate_predictors$labels_plus_both
)
full_models <- list()
importance <- tibble()
for (model_context in names(model_specs)) {
    model <- fit_brt(
        data, model_specs[[model_context]],
        20261116L + match(model_context, names(model_specs))
    )
    full_models[[model_context]] <- model
    importance <- bind_rows(
        importance,
        summary(model, plotit = FALSE) |>
            as_tibble() |>
            transmute(
                model_context, predictor = var,
                relative_influence = rel.inf,
                rank = rank(-rel.inf, ties.method = 'min')
            )
    )
}
write_csv(
    importance,
    file.path(output_dir, 'cause_adjusted_salinity_brt_importance.csv'),
    na = ''
)
saveRDS(
    full_models,
    file.path(output_dir, 'cause_adjusted_salinity_brt_models.rds')
)

effect_metrics <- c(
    'local_first_dhw', 'cots_hazard', 'cyclone_wave_hours',
    'sss_min', 'hours_below30'
)
metric_labels <- c(
    local_first_dhw = 'Local-first DHW',
    cots_hazard = 'COTS hazard index',
    cyclone_wave_hours = 'Cyclone-wave exposure (hours >4 m)',
    sss_min = 'Minimum surface salinity (PSU)',
    hours_below30 = 'Hours below 30 PSU'
)
strata <- tribble(
    ~programme_key, ~depth, ~is_manta, ~is_mmp, ~programme_depth,
    'mmp', 2, 0, 1, 'MMP 2 m',
    'mmp', 5, 0, 1, 'MMP 5 m',
    'manta', 9, 1, 0, 'Manta 9 m',
    'ltmp', 9, 0, 0, 'LTMP 9 m'
)
metric_grids <- lapply(effect_metrics, function(metric) {
    unique(as.numeric(quantile(
        data[[metric]], probs = seq(0, 1, length.out = 35L),
        na.rm = TRUE, names = FALSE
    )))
})
names(metric_grids) <- effect_metrics

partial_curves <- function(model, integration_data, model_context) {
    bind_rows(lapply(effect_metrics, function(metric) {
        grid <- metric_grids[[metric]]
        bind_rows(lapply(seq_len(nrow(strata)), function(i) {
            stratum <- strata[i, ]
            bind_rows(lapply(grid, function(metric_value) {
                new_data <- integration_data
                new_data$depth <- stratum$depth
                new_data$is_manta <- stratum$is_manta
                new_data$is_mmp <- stratum$is_mmp
                new_data[[metric]] <- metric_value
                tibble(
                    model_context, metric, metric_value,
                    programme_depth = stratum$programme_depth,
                    predicted_mortality = mean(
                        predict_brt(model, new_data)
                    )
                )
            }))
        }))
    }))
}

label_effects <- function(model, integration_data) {
    bind_rows(lapply(label_predictors, function(label) {
        bind_rows(lapply(seq_len(nrow(strata)), function(i) {
            stratum <- strata[i, ]
            absent <- integration_data
            absent$depth <- stratum$depth
            absent$is_manta <- stratum$is_manta
            absent$is_mmp <- stratum$is_mmp
            present <- absent
            absent[[label]] <- 0
            present[[label]] <- 1
            prediction_absent <- mean(predict_brt(model, absent))
            prediction_present <- mean(predict_brt(model, present))
            tibble(
                label, programme_depth = stratum$programme_depth,
                prediction_absent, prediction_present,
                marginal_difference =
                    prediction_present - prediction_absent
            )
        }))
    }))
}

point_curves <- bind_rows(lapply(
    names(full_models),
    function(model_context) partial_curves(
        full_models[[model_context]], data, model_context
    )
))
point_label_effects <- label_effects(
    full_models$cause_exposures_and_labels, data
)

bootstrap_curves <- tibble()
bootstrap_labels <- tibble()
reef_ids <- unique(data$ReefID)
set.seed(20261216L)
for (model_context in names(model_specs)) {
    successes <- 0L
    for (b in seq_len(60L)) {
        sampled <- sample(reef_ids, length(reef_ids), replace = TRUE)
        index <- unlist(lapply(
            sampled, function(id) which(data$ReefID == id)
        ))
        boot_data <- data[index, , drop = FALSE]
        boot_model <- tryCatch(
            fit_brt(
                boot_data, model_specs[[model_context]],
                20262016L +
                    1000L * match(model_context, names(model_specs)) + b,
                n_trees = 500L
            ),
            error = function(e) NULL
        )
        if (is.null(boot_model)) next
        successes <- successes + 1L
        bootstrap_curves <- bind_rows(
            bootstrap_curves,
            partial_curves(
                boot_model, boot_data, model_context
            ) |>
                mutate(bootstrap = b)
        )
        if (model_context == 'cause_exposures_and_labels') {
            bootstrap_labels <- bind_rows(
                bootstrap_labels,
                label_effects(boot_model, boot_data) |>
                    mutate(bootstrap = b)
            )
        }
    }
    if (successes < 50L) {
        stop('Fewer than 50 successful effect bootstraps for ', model_context)
    }
}

curve_intervals <- bootstrap_curves |>
    group_by(model_context, metric, metric_value, programme_depth) |>
    summarise(
        lower = quantile(predicted_mortality, 0.025),
        upper = quantile(predicted_mortality, 0.975),
        .groups = 'drop'
    )
effect_curves <- point_curves |>
    left_join(
        curve_intervals,
        by = c(
            'model_context', 'metric', 'metric_value',
            'programme_depth'
        )
    ) |>
    mutate(
        metric_label = recode(metric, !!!metric_labels),
        model_context_label = recode(
            model_context,
            cause_exposures = 'Exposure covariates only',
            cause_exposures_and_labels =
                'Exposure covariates + retrospective labels'
        )
    )
write_csv(
    effect_curves,
    file.path(
        output_dir, 'cause_adjusted_salinity_brt_partial_dependence.csv'
    ),
    na = ''
)

label_intervals <- bootstrap_labels |>
    group_by(label, programme_depth) |>
    summarise(
        lower = quantile(marginal_difference, 0.025),
        upper = quantile(marginal_difference, 0.975),
        .groups = 'drop'
    )
label_marginal_effects <- point_label_effects |>
    left_join(
        label_intervals, by = c('label', 'programme_depth')
    ) |>
    mutate(
        label_name = recode(
            label,
            bleaching_label = 'Bleaching recorded',
            flood_label = 'Flood recorded',
            cyclone_label = 'Cyclone/storm recorded',
            cots_label = 'COTS recorded'
        )
    )
write_csv(
    label_marginal_effects,
    file.path(
        output_dir, 'cause_adjusted_salinity_brt_label_effects.csv'
    ),
    na = ''
)

effect_plot <- ggplot(
    effect_curves,
    aes(metric_value, predicted_mortality, colour = programme_depth,
        fill = programme_depth)
) +
    geom_ribbon(
        aes(ymin = lower, ymax = upper), alpha = 0.08, colour = NA
    ) +
    geom_line(linewidth = 0.75) +
    facet_grid(
        model_context_label ~ metric_label,
        scales = 'free_x'
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
        x = NULL, y = 'Partial expected relative mortality',
        colour = 'Programme/depth', fill = 'Programme/depth',
        title = paste(
            'Cause-adjusted BRT marginal relationships across programmes'
        ),
        subtitle = paste(
            'Bands are 60 reef-cluster bootstraps;',
            'field labels are retrospective explanatory sensitivities'
        )
    ) +
    theme_bw(base_size = 9) +
    theme(
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 35, hjust = 1)
    )
ggsave(
    file.path(
        output_dir, 'cause_adjusted_salinity_brt_partial_dependence.png'
    ),
    effect_plot, width = 18, height = 8, dpi = 180
)

label_plot <- ggplot(
    label_marginal_effects,
    aes(marginal_difference, programme_depth, colour = programme_depth)
) +
    geom_vline(xintercept = 0, colour = 'grey55', linetype = 2) +
    geom_errorbar(
        aes(xmin = lower, xmax = upper), orientation = 'y', width = 0.18
    ) +
    geom_point(size = 2) +
    facet_wrap(~ label_name, scales = 'free_x') +
    scale_x_continuous(labels = scales::percent) +
    labs(
        x = 'Marginal change in expected mortality when label is present',
        y = NULL, colour = 'Programme/depth',
        title = 'Retrospective field-label effects after continuous exposures'
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'none')
ggsave(
    file.path(output_dir, 'cause_adjusted_salinity_brt_label_effects.png'),
    label_plot, width = 12, height = 6, dpi = 180
)

cat(
    'Wrote cause-adjusted salinity BRT outputs using', nrow(data),
    'programme-depth rows from', n_distinct(data$ReefID), 'reefs.\n'
)
print(metrics |> arrange(rmse))
print(metric_deltas)
print(importance |> filter(predictor %in% c(
    'local_first_dhw', 'cots_hazard', 'cyclone_wave_hours',
    'sss_min', 'hours_below30'
)) |> arrange(model_context, rank))
