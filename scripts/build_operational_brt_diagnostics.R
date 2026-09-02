# Operational BRT interpretation suite. Event identity is deliberately
# excluded. Bootstrap intervals refit the complete BRT using reef-event cluster
# resampling and summarise one-dimensional partial dependence uncertainty.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(ggplot2)
    library(purrr)
    library(readr)
    library(tidyr)
})
source('scripts/model_registry.R')
source('scripts/model_diagnostics.R')

root <- project_root()
output_dir <- file.path(root, 'output', 'operational_brt_diagnostics')
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(20260901L)

rows <- read_csv(
    file.path(root, 'output', 'explanatory_event_dhw', 'event_dhw_brt_data.csv'),
    show_col_types = FALSE
)

# Collinearity-screened operational inputs. Trees learn thresholds and
# interactions directly, so explicit DHW hinge/interactions are not duplicated.
predictors <- c(
    'ann_maxdhw', 'prop_acropora_pre', 'observed_pre_cover',
    'dhw_novelty10', 'dhw_events_since2016_n6',
    'dhw_years_since_last_n6_capped8', 'secc3m_p10', 'cloudp_90',
    'mcur_90', 'sst_summer_skewness', 'sst_summer_excess_kurtosis',
    'chla_wetseason_median', 'log_coastal_rain30',
    'wqc_freqcc12', 'wqc_prior10_percentile',
    'log1p_cyc_interval_maxHrs4mw', 'tc_interval_proximity100',
    'tc_interval_wind_distance_index', 'log1p_cot_interval_idw_max',
    'cots_outbreak_probability', 'programme_factor'
)
stopifnot(all(predictors %in% names(rows)))

model_data <- rows |>
    select(mortality_prop, reef_event_key, all_of(predictors))
for (variable in setdiff(predictors, 'programme_factor')) {
    replacement <- median(model_data[[variable]], na.rm = TRUE)
    model_data[[variable]][!is.finite(model_data[[variable]])] <- replacement
}
model_data$programme_factor <- droplevels(factor(model_data$programme_factor))

reef_event_weights <- function(data) {
    event_n <- table(data$reef_event_key)
    programme_n <- table(data$programme_factor)
    weights <- 1 / as.numeric(event_n[data$reef_event_key])
    weights <- weights / as.numeric(programme_n[data$programme_factor])
    weights / mean(weights)
}
fit_brt <- function(data, seed) {
    set.seed(seed)
    gbm(
        as.formula(paste('mortality_prop ~', paste(predictors, collapse = ' + '))),
        data = data, distribution = 'gaussian',
        weights = reef_event_weights(data),
        n.trees = 1200, interaction.depth = 3, shrinkage = 0.01,
        n.minobsinnode = 10, bag.fraction = 0.70,
        train.fraction = 1, keep.data = TRUE, verbose = FALSE
    )
}

production_fit <- fit_brt(model_data, 20260901L)
saveRDS(production_fit, file.path(output_dir, 'operational_brt.rds'))
influence <- summary(production_fit, n.trees = 1200, plotit = FALSE) |>
    as_tibble() |>
    transmute(variable = var, relative_influence = rel.inf) |>
    arrange(desc(relative_influence))
write_csv(influence, file.path(output_dir, 'variable_importance.csv'))

labels <- c(
    ann_maxdhw = 'Local-first maximum DHW',
    prop_acropora_pre = 'Pre-event Acropora proportion',
    observed_pre_cover = 'Pre-event coral cover',
    dhw_novelty10 = 'DHW novelty',
    dhw_events_since2016_n6 = 'Events above 6 DHW since 2016',
    dhw_years_since_last_n6_capped8 = 'Years since last event above 6 DHW',
    secc3m_p10 = 'Low-tail Secchi depth (p10)',
    cloudp_90 = 'High cloud cover (p90)',
    mcur_90 = 'High current speed (p90)',
    sst_summer_skewness = 'Summer SST skewness',
    sst_summer_excess_kurtosis = 'Summer SST excess kurtosis',
    chla_wetseason_median = 'Wet-season median chlorophyll-a',
    log_coastal_rain30 = 'Coastal 30-day rainfall (log)',
    wqc_freqcc12 = 'Current coloured-water frequency',
    wqc_prior10_percentile = 'Reef-relative coloured-water percentile',
    log1p_cyc_interval_maxHrs4mw = 'Interval damaging-wave exposure (log)',
    tc_interval_proximity100 = 'Interval cyclone proximity',
    tc_interval_wind_distance_index = 'Interval cyclone wind-distance index',
    log1p_cot_interval_idw_max = 'Interval modelled COTS pressure (log)',
    cots_outbreak_probability = 'COTS outbreak probability',
    programme_factor = 'Survey programme'
)
importance_plot_data <- influence |>
    mutate(
        label = unname(labels[variable]),
        label = coalesce(label, variable),
        label = factor(label, levels = rev(label))
    )
importance_plot <- ggplot(
    importance_plot_data,
    aes(relative_influence, label)
) +
    geom_col(fill = '#2A9D8F', alpha = 0.9) +
    geom_text(
        aes(label = sprintf('%.1f%%', relative_influence)),
        hjust = -0.1, size = 3
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
        title = 'Operational BRT variable importance',
        subtitle = 'No event identity; collinearity-screened, operationally obtainable predictors',
        x = 'Relative influence (%)', y = NULL
    ) +
    theme_bw(base_size = 11)
save_figure_bundle(
    importance_plot, 'Fig-BRT-03_operational_variable_importance',
    importance_plot_data,
    'Relative influence of predictors in the full-data operational BRT fitted without event identity.',
    'Shows which variables most often improve tree splits; correlated predictors can share influence.',
    'Descriptive full-data importance is not held-out skill and does not establish causal importance.',
    'operational_brt_collinearity_screened', 'BRT', 'variable_importance',
    'operational_diagnostic', root, TRUE,
    code_source = 'scripts/build_operational_brt_diagnostics.R',
    width = 10, height = 8
)

numeric_importance <- influence |>
    filter(variable != 'programme_factor')
pdp_variables <- head(numeric_importance$variable, 12)
grid_size <- 45L
bootstrap_n <- 50L
grids <- setNames(lapply(pdp_variables, function(variable) {
    limits <- quantile(model_data[[variable]], c(0.02, 0.98), na.rm = TRUE)
    seq(limits[[1]], limits[[2]], length.out = grid_size)
}), pdp_variables)
bootstrap_values <- setNames(lapply(pdp_variables, function(variable) {
    matrix(NA_real_, nrow = bootstrap_n, ncol = grid_size)
}), pdp_variables)

cluster_ids <- unique(model_data$reef_event_key)
for (bootstrap_index in seq_len(bootstrap_n)) {
    sampled_clusters <- sample(cluster_ids, length(cluster_ids), replace = TRUE)
    bootstrap_rows <- bind_rows(lapply(
        seq_along(sampled_clusters),
        function(copy_index) {
            model_data |>
                filter(reef_event_key == sampled_clusters[[copy_index]]) |>
                mutate(reef_event_key = paste0(reef_event_key, '__', copy_index))
        }
    ))
    fitted <- fit_brt(bootstrap_rows, 20260901L + bootstrap_index)
    for (variable in pdp_variables) {
        variable_index <- match(variable, fitted$var.names)
        partial <- gbm::plot.gbm(
            fitted, i.var = variable_index, n.trees = 1200,
            return.grid = TRUE
        )
        bootstrap_values[[variable]][bootstrap_index, ] <- approx(
            x = partial[[1]], y = partial[[2]], xout = grids[[variable]],
            rule = 2, ties = mean
        )$y
    }
}

pdp <- map_dfr(pdp_variables, function(variable) {
    values <- bootstrap_values[[variable]]
    tibble(
        variable = variable,
        label = unname(labels[variable]),
        predictor_value = grids[[variable]],
        median = apply(values, 2, median, na.rm = TRUE),
        lower95 = apply(values, 2, quantile, 0.025, na.rm = TRUE),
        upper95 = apply(values, 2, quantile, 0.975, na.rm = TRUE),
        bootstrap_fits = bootstrap_n
    )
})
write_csv(pdp, file.path(output_dir, 'bootstrapped_partial_dependence.csv'))
pdp_plot <- ggplot(pdp, aes(predictor_value, median)) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = '#56B4E9', alpha = 0.28) +
    geom_line(colour = '#0072B2', linewidth = 0.9) +
    facet_wrap(~ label, scales = 'free_x', ncol = 3) +
    labs(
        title = 'Operational BRT partial dependence with 95% bootstrap intervals',
        subtitle = 'Reef-event cluster bootstrap; B = 50 complete BRT refits',
        x = 'Predictor value', y = 'Partial dependence (mortality proportion)'
    ) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(size = 9))
save_figure_bundle(
    pdp_plot, 'Fig-BRT-04_operational_bootstrap_pdp', pdp,
    'One-dimensional partial-dependence curves with 95% reef-event cluster-bootstrap intervals for the twelve most influential numeric operational BRT predictors.',
    'Shows the fitted marginal shape and sampling stability of influential BRT relationships.',
    'Intervals quantify cluster-resampling uncertainty in the full-data fitted relationship; they are not prediction intervals or proof of causality.',
    'operational_brt_collinearity_screened', 'BRT', 'partial_dependence',
    'operational_diagnostic', root, TRUE,
    code_source = 'scripts/build_operational_brt_diagnostics.R',
    width = 12, height = 10
)
write_figure_readme(root)
message('Wrote operational BRT diagnostics')
