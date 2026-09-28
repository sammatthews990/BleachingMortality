# Minimal MMP-only direct BRT for the two supplied surface-salinity metrics.
#
# This is a deliberately small diagnostic, not a production-model candidate.
# It averages multiple mortality transitions to one ReefID-depth outcome, uses
# the existing reef-blocked folds, and includes no predictors other than depth,
# area-mean minimum SSS and area-mean hours below 30 PSU.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(ggplot2)
    library(readr)
    library(tidyr)
})

set.seed(20260916L)

input_file <- file.path(
    'output', 'surface_salinity_mortality',
    'cause_focused_reef_dataset.csv'
)
output_dir <- file.path('output', 'surface_salinity_mortality')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
    stop(
        'Missing cause-focused dataset. Run ',
        'src/evaluation/screen_surface_salinity_mortality.R first.'
    )
}

rows <- read_csv(input_file, show_col_types = FALSE) |>
    filter(
        programme_key == 'mmp', depth %in% c(2, 5),
        is.finite(mortality_prop),
        is.finite(sss_min_area_mean_psu),
        is.finite(hours_below30_area_mean)
    )

# Give every reef-depth combination equal weight. Some reefs have two source
# mortality transitions while others have one, but salinity is reef-level.
mmp <- rows |>
    group_by(
        ReefID, ReefName, depth, joint_reef_fold,
        sss_min_area_mean_psu, hours_below30_area_mean
    ) |>
    summarise(
        mortality = mean(mortality_prop),
        source_rows = n(),
        flood_label = any(reef_has_flood_description),
        description_classes = paste(
            sort(unique(target_description_class)), collapse = '; '
        ),
        .groups = 'drop'
    ) |>
    arrange(ReefID, depth)

if (n_distinct(mmp$ReefID) < 10L || !all(c(2, 5) %in% mmp$depth)) {
    stop('Insufficient paired MMP reef-depth support for the simple BRT')
}
if (anyDuplicated(mmp[c('ReefID', 'depth')])) {
    stop('MMP diagnostic data are not unique by ReefID-depth')
}

write_csv(
    mmp,
    file.path(output_dir, 'simple_mmp_salinity_brt_data.csv'), na = ''
)

candidate_predictors <- list(
    depth_only = 'depth',
    depth_plus_sss = c('depth', 'sss_min_area_mean_psu'),
    depth_plus_hours = c('depth', 'hours_below30_area_mean'),
    depth_plus_both = c(
        'depth', 'sss_min_area_mean_psu', 'hours_below30_area_mean'
    )
)

fit_brt <- function(data, predictors, seed, n_trees = 500L) {
    active_predictors <- predictors[vapply(
        predictors,
        function(predictor) length(unique(data[[predictor]])) > 1L,
        logical(1)
    )]
    if (length(active_predictors) == 0L) {
        stop('No varying predictors in BRT analysis sample')
    }
    set.seed(seed)
    gbm(
        reformulate(active_predictors, response = 'mortality'),
        data = data,
        distribution = 'gaussian',
        n.trees = n_trees,
        interaction.depth = 2L,
        shrinkage = 0.02,
        n.minobsinnode = 2L,
        bag.fraction = 0.75,
        train.fraction = 1,
        keep.data = FALSE,
        verbose = FALSE
    )
}

predict_brt <- function(model, new_data, n_trees = 500L) {
    pmin(pmax(
        as.numeric(predict(model, new_data, n.trees = n_trees)), 0
    ), 1)
}

# The original joint reef folds are preserved. No row from a held-out reef is
# used to fit its prediction at either depth.
predictions <- tibble()
for (candidate in names(candidate_predictors)) {
    predictors <- candidate_predictors[[candidate]]
    for (fold in sort(unique(mmp$joint_reef_fold))) {
        analysis <- mmp |> filter(joint_reef_fold != fold)
        assessment <- mmp |> filter(joint_reef_fold == fold)
        model <- fit_brt(
            analysis, predictors,
            20260916L + 100L * match(candidate, names(candidate_predictors)) +
                fold
        )
        predictions <- bind_rows(
            predictions,
            assessment |>
                transmute(
                    ReefID, ReefName, depth, joint_reef_fold,
                    observed_mortality = mortality,
                    flood_label,
                    candidate,
                    predicted_mortality = predict_brt(model, assessment)
                )
        )
    }
}

metrics <- predictions |>
    group_by(candidate) |>
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
        severe_rmse = sqrt(mean(
            (observed_mortality[observed_mortality >= 0.2] -
                predicted_mortality[observed_mortality >= 0.2])^2
        )),
        false_extreme_rate = mean(
            predicted_mortality >= 0.5 & observed_mortality < 0.2
        ),
        .groups = 'drop'
    )

metrics_by_flood_label <- predictions |>
    group_by(candidate, flood_label) |>
    summarise(
        rows = n(), reefs = n_distinct(ReefID),
        rmse = sqrt(mean(
            (observed_mortality - predicted_mortality)^2
        )),
        mae = mean(abs(observed_mortality - predicted_mortality)),
        bias = mean(predicted_mortality - observed_mortality),
        .groups = 'drop'
    )

bootstrap_metric_deltas <- function(prediction_rows, replicates = 2000L) {
    baseline <- prediction_rows |>
        filter(candidate == 'depth_only') |>
        select(
            ReefID, depth, observed_mortality,
            baseline_prediction = predicted_mortality
        )
    comparison <- prediction_rows |>
        filter(candidate != 'depth_only') |>
        inner_join(
            baseline,
            by = c('ReefID', 'depth', 'observed_mortality'),
            relationship = 'many-to-one'
        )
    bind_rows(lapply(
        setdiff(names(candidate_predictors), 'depth_only'),
        function(candidate_name) {
            x <- comparison |> filter(candidate == candidate_name)
            reefs <- unique(x$ReefID)
            set.seed(20261016L + match(
                candidate_name, names(candidate_predictors)
            ))
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
                candidate = candidate_name,
                delta_rmse = sqrt(mean(
                    (x$observed_mortality - x$predicted_mortality)^2
                )) - sqrt(mean(
                    (x$observed_mortality - x$baseline_prediction)^2
                )),
                delta_rmse_q025 = quantile(delta, 0.025),
                delta_rmse_q975 = quantile(delta, 0.975),
                bootstrap_probability_improved = mean(delta < 0)
            )
        }
    ))
}

metric_deltas <- bootstrap_metric_deltas(predictions)
write_csv(
    predictions,
    file.path(output_dir, 'simple_mmp_salinity_brt_predictions.csv'),
    na = ''
)
write_csv(
    metrics,
    file.path(output_dir, 'simple_mmp_salinity_brt_metrics.csv'),
    na = ''
)
write_csv(
    metrics_by_flood_label,
    file.path(
        output_dir, 'simple_mmp_salinity_brt_metrics_by_flood_label.csv'
    ),
    na = ''
)
write_csv(
    metric_deltas,
    file.path(output_dir, 'simple_mmp_salinity_brt_metric_deltas.csv'),
    na = ''
)

full_model <- fit_brt(
    mmp, candidate_predictors$depth_plus_both, 20261116L
)
importance <- summary(full_model, plotit = FALSE) |>
    as_tibble() |>
    transmute(
        predictor = var,
        relative_influence = rel.inf,
        rank = rank(-rel.inf, ties.method = 'min')
    )
write_csv(
    importance,
    file.path(output_dir, 'simple_mmp_salinity_brt_importance.csv'),
    na = ''
)
saveRDS(
    full_model,
    file.path(output_dir, 'simple_mmp_salinity_brt_model.rds')
)

partial_curve <- function(model, integration_data, metric) {
    metric_grid <- sort(unique(integration_data[[metric]]))
    bind_rows(lapply(c(2, 5), function(target_depth) {
        bind_rows(lapply(metric_grid, function(metric_value) {
            new_data <- integration_data
            new_data$depth <- target_depth
            new_data[[metric]] <- metric_value
            tibble(
                depth = target_depth,
                metric = metric,
                metric_value = metric_value,
                predicted_mortality = mean(predict_brt(model, new_data))
            )
        }))
    }))
}

metrics_for_curves <- c(
    'sss_min_area_mean_psu', 'hours_below30_area_mean'
)
point_curves <- bind_rows(lapply(
    metrics_for_curves,
    function(metric) partial_curve(full_model, mmp, metric)
))

set.seed(20261216L)
reef_ids <- unique(mmp$ReefID)
bootstrap_curves <- tibble()
successful_bootstraps <- 0L
for (b in seq_len(250L)) {
    sampled <- sample(reef_ids, length(reef_ids), replace = TRUE)
    index <- unlist(lapply(sampled, function(id) which(mmp$ReefID == id)))
    boot_data <- mmp[index, , drop = FALSE]
    boot_model <- tryCatch(
        fit_brt(
            boot_data, candidate_predictors$depth_plus_both,
            20262016L + b
        ),
        error = function(e) NULL
    )
    if (is.null(boot_model)) next
    successful_bootstraps <- successful_bootstraps + 1L
    bootstrap_curves <- bind_rows(
        bootstrap_curves,
        bind_rows(lapply(
            metrics_for_curves,
            function(metric) partial_curve(boot_model, boot_data, metric)
        )) |>
            mutate(bootstrap = b)
    )
}
if (successful_bootstraps < 200L) {
    stop('Fewer than 200 successful reef-bootstrap BRT fits')
}

curve_intervals <- bootstrap_curves |>
    group_by(depth, metric, metric_value) |>
    summarise(
        lower = quantile(predicted_mortality, 0.025),
        upper = quantile(predicted_mortality, 0.975),
        .groups = 'drop'
    )
partial_dependence <- point_curves |>
    left_join(
        curve_intervals,
        by = c('depth', 'metric', 'metric_value')
    ) |>
    mutate(
        metric_label = recode(
            metric,
            sss_min_area_mean_psu =
                'Area-mean minimum surface salinity (PSU)',
            hours_below30_area_mean =
                'Area-mean hours below 30 PSU'
        ),
        depth_label = paste0(depth, ' m')
    )
write_csv(
    partial_dependence,
    file.path(output_dir, 'simple_mmp_salinity_brt_partial_dependence.csv'),
    na = ''
)

raw_long <- mmp |>
    pivot_longer(
        all_of(metrics_for_curves),
        names_to = 'metric', values_to = 'metric_value'
    ) |>
    mutate(
        metric_label = recode(
            metric,
            sss_min_area_mean_psu =
                'Area-mean minimum surface salinity (PSU)',
            hours_below30_area_mean =
                'Area-mean hours below 30 PSU'
        ),
        depth_label = paste0(depth, ' m')
    )

partial_plot <- ggplot(
    partial_dependence,
    aes(metric_value, predicted_mortality, colour = depth_label,
        fill = depth_label)
) +
    geom_ribbon(
        aes(ymin = lower, ymax = upper), alpha = 0.12, colour = NA
    ) +
    geom_step(linewidth = 0.9) +
    geom_point(
        data = raw_long,
        aes(metric_value, mortality, colour = depth_label,
            shape = flood_label),
        inherit.aes = FALSE, alpha = 0.65, size = 2
    ) +
    facet_wrap(~ metric_label, scales = 'free_x', nrow = 1) +
    scale_colour_manual(values = c('2 m' = '#009EAA', '5 m' = '#8E5BD9')) +
    scale_fill_manual(values = c('2 m' = '#009EAA', '5 m' = '#8E5BD9')) +
    scale_y_continuous(
        limits = c(0, 1), labels = scales::percent
    ) +
    labs(
        x = NULL, y = 'Relative mortality',
        colour = 'Survey depth', fill = 'Survey depth',
        shape = 'Flood-labelled reef',
        title = 'Minimal MMP BRT: depth plus two surface-salinity metrics',
        subtitle = paste(
            'Lines are partial dependence; bands are 250 reef bootstraps;',
            'points are 14 reef-depth means'
        )
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = 'bottom')
ggsave(
    file.path(output_dir, 'simple_mmp_salinity_brt_partial_dependence.png'),
    partial_plot, width = 11, height = 5.5, dpi = 180
)

association <- tibble(
    variable_1 = c(
        'sss_min_area_mean_psu', 'sss_min_area_mean_psu',
        'hours_below30_area_mean'
    ),
    variable_2 = c(
        'mortality', 'hours_below30_area_mean', 'mortality'
    ),
    spearman_rho = c(
        cor(
            mmp$sss_min_area_mean_psu, mmp$mortality,
            method = 'spearman'
        ),
        cor(
            mmp$sss_min_area_mean_psu, mmp$hours_below30_area_mean,
            method = 'spearman'
        ),
        cor(
            mmp$hours_below30_area_mean, mmp$mortality,
            method = 'spearman'
        )
    )
)
write_csv(
    association,
    file.path(output_dir, 'simple_mmp_salinity_brt_associations.csv'),
    na = ''
)

cat(
    'Wrote simple MMP salinity BRT outputs using', nrow(mmp),
    'reef-depth rows from', n_distinct(mmp$ReefID), 'reefs.\n'
)
print(metrics |> arrange(rmse))
print(metric_deltas)
print(importance)
