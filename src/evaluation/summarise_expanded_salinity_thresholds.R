# Posterior marginal contrasts and screening thresholds for the expanded INLA
# salinity curves. Thresholds are grid crossings relative to no exposure (or
# the highest SSS reference), not biological lethal limits.

suppressPackageStartupMessages({
    library(dplyr)
    library(INLA)
    library(readr)
    library(splines)
})

set.seed(20260928L)
output_dir <- 'output/surface_salinity_mortality'
models <- readRDS(file.path(output_dir, 'expanded_salinity_inla_models.rds'))
files <- c(
    bleaching_mortality = file.path(
        output_dir, 'expanded_bleaching_salinity_data.csv'
    ),
    annual_transition_sensitivity = file.path(
        output_dir, 'expanded_annual_transition_salinity_data.csv'
    )
)
base_names <- c(
    'thermal_dhw', 'cots_hazard', 'cyclone_wave_hours', 'pre_cover'
)
metric_names <- c('sss_min', 'log_hours_below30', 'log_hours_below26')
strata <- tibble(
    stratum = c('LTMP 9 m', 'Manta 9 m', 'MMP 2 m', 'MMP 5 m'),
    depth = c(9, 9, 2, 5),
    is_manta = c(0, 1, 0, 0),
    is_mmp = c(0, 0, 1, 1)
)

prepare_context <- function(data) {
    data <- data |>
        mutate(
            ReefID = toupper(trimws(ReefID)),
            depth = as.numeric(depth),
            is_manta = as.numeric(programme_key == 'manta'),
            is_mmp = as.numeric(programme_key == 'mmp'),
            reef_index = match(ReefID, unique(ReefID))
        )
    for (variable in base_names) {
        observed <- data[[variable]][is.finite(data[[variable]])]
        fill <- if (length(observed)) median(observed) else 0
        data[[variable]][!is.finite(data[[variable]])] <- fill
        center <- mean(data[[variable]])
        spread <- sd(data[[variable]])
        if (!is.finite(spread) || spread == 0) spread <- 1
        data[[paste0(variable, '_z')]] <-
            (data[[variable]] - center) / spread
    }
    data
}
add_metric_basis <- function(x, metric, spec) {
    basis <- ns(
        x[[metric]], knots = spec$knot,
        Boundary.knots = spec$bounds
    )
    basis <- sweep(basis, 2, spec$center, '-')
    basis <- sweep(basis, 2, spec$scale, '/')
    x[[paste0(metric, '_ns1')]] <- basis[, 1]
    x[[paste0(metric, '_ns2')]] <- basis[, 2]
    x
}
fixed_draws <- function(fit, samples = 1200L) {
    draws <- inla.posterior.sample(samples, fit)
    fixed_names <- rownames(fit$summary.fixed)
    latent_names <- paste0(fixed_names, ':1')
    beta <- vapply(
        draws,
        function(draw) as.numeric(draw$latent[latent_names, 1]),
        numeric(length(latent_names))
    )
    rownames(beta) <- fixed_names
    beta
}

contrasts <- tibble()
for (context in names(files)) {
    data <- prepare_context(read_csv(files[[context]], show_col_types = FALSE))
    for (metric in metric_names) {
        spec <- models[[context]]$spline_specs[[metric]]
        fixed_terms <- models[[context]]$model_specs[[metric]]
        fixed_formula <- reformulate(fixed_terms)
        fit <- models[[context]]$fits[[metric]]
        beta_occurrence <- fixed_draws(fit$occurrence)
        beta_magnitude <- fixed_draws(fit$magnitude)
        grid <- unique(as.numeric(quantile(
            data[[metric]], probs = seq(0, 1, length.out = 81),
            na.rm = TRUE
        )))
        reference_value <- if (metric == 'sss_min') max(grid) else 0
        for (s in seq_len(nrow(strata))) {
            reference_data <- data
            reference_data[[metric]] <- reference_value
            reference_data$depth <- strata$depth[[s]]
            reference_data$is_manta <- strata$is_manta[[s]]
            reference_data$is_mmp <- strata$is_mmp[[s]]
            reference_data <- add_metric_basis(
                reference_data, metric, spec
            )
            reference_design <- model.matrix(
                fixed_formula, reference_data
            )
            reference_draw <- colMeans(
                plogis(reference_design %*% beta_occurrence) *
                    plogis(reference_design %*% beta_magnitude)
            )
            for (value in grid) {
                new_data <- data
                new_data[[metric]] <- value
                new_data$depth <- strata$depth[[s]]
                new_data$is_manta <- strata$is_manta[[s]]
                new_data$is_mmp <- strata$is_mmp[[s]]
                new_data <- add_metric_basis(new_data, metric, spec)
                design <- model.matrix(fixed_formula, new_data)
                mortality_draw <- colMeans(
                    plogis(design %*% beta_occurrence) *
                        plogis(design %*% beta_magnitude)
                )
                difference <- mortality_draw - reference_draw
                contrasts <- bind_rows(
                    contrasts,
                    tibble(
                        analysis_context = context,
                        metric, metric_value = value,
                        original_value = if_else(
                            metric %in% c(
                                'log_hours_below30', 'log_hours_below26'
                            ),
                            expm1(value), value
                        ),
                        reference_value = if_else(
                            metric %in% c(
                                'log_hours_below30', 'log_hours_below26'
                            ),
                            expm1(reference_value), reference_value
                        ),
                        stratum = strata$stratum[[s]],
                        predicted_mortality = mean(mortality_draw),
                        predicted_lower = quantile(mortality_draw, 0.025),
                        predicted_upper = quantile(mortality_draw, 0.975),
                        difference_from_reference = mean(difference),
                        difference_lower = quantile(difference, 0.025),
                        difference_upper = quantile(difference, 0.975)
                    )
                )
            }
        }
    }
}

thresholds <- contrasts |>
    group_by(analysis_context, metric, stratum) |>
    group_modify(function(x, keys) {
        supported <- x |> filter(difference_lower > 0)
        if (nrow(supported) == 0L) {
            return(tibble(
                threshold_metric_value = NA_real_,
                threshold_original_value = NA_real_,
                difference_at_threshold = NA_real_,
                difference_lower = NA_real_, difference_upper = NA_real_
            ))
        }
        selected <- if (keys$metric == 'sss_min') {
            supported |> slice_max(metric_value, n = 1, with_ties = FALSE)
        } else {
            supported |> slice_min(metric_value, n = 1, with_ties = FALSE)
        }
        selected |>
            transmute(
                threshold_metric_value = metric_value,
                threshold_original_value = original_value,
                difference_at_threshold = difference_from_reference,
                difference_lower, difference_upper
            )
    }) |>
    ungroup()

support_rows <- bind_rows(lapply(names(files), function(context) {
    data <- read_csv(files[[context]], show_col_types = FALSE)
    thresholds |>
        filter(analysis_context == context) |>
        rowwise() |>
        mutate(
            empirical_rows_beyond_threshold = if (is.na(
                threshold_original_value
            )) NA_integer_ else if (metric == 'sss_min') {
                sum(data$sss_min <= threshold_original_value)
            } else if (metric == 'log_hours_below30') {
                sum(data$hours_below30 >= threshold_original_value)
            } else {
                sum(data$hours_below26 >= threshold_original_value)
            },
            empirical_reefs_beyond_threshold = if (is.na(
                threshold_original_value
            )) NA_integer_ else if (metric == 'sss_min') {
                n_distinct(data$ReefID[
                    data$sss_min <= threshold_original_value
                ])
            } else if (metric == 'log_hours_below30') {
                n_distinct(data$ReefID[
                    data$hours_below30 >= threshold_original_value
                ])
            } else {
                n_distinct(data$ReefID[
                    data$hours_below26 >= threshold_original_value
                ])
            }
        ) |>
        ungroup()
}))

write_csv(
    contrasts,
    file.path(output_dir, 'expanded_salinity_inla_marginal_contrasts.csv')
)
write_csv(
    support_rows,
    file.path(output_dir, 'expanded_salinity_inla_threshold_summary.csv')
)
cat('Wrote posterior marginal contrasts and screening thresholds.\n')
print(support_rows)
