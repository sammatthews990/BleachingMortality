# Additive two-part INLA diagnostic for the expanded salinity screen.
# This full-data explanatory fit smooths the BRT relationships; WAIC/DIC are
# supplementary and do not replace leave-one-event-out BRT evidence.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(INLA)
    library(readr)
    library(splines)
})

set.seed(20260928L)
output_dir <- 'output/surface_salinity_mortality'
files <- c(
    bleaching_mortality = file.path(
        output_dir, 'expanded_bleaching_salinity_data.csv'
    ),
    annual_transition_sensitivity = file.path(
        output_dir, 'expanded_annual_transition_salinity_data.csv'
    )
)
if (!all(file.exists(files))) stop('Run expanded salinity BRT screen first')

metric_names <- c('sss_min', 'log_hours_below30', 'log_hours_below26')
metric_labels <- c(
    sss_min = 'Minimum surface salinity (PSU)',
    log_hours_below30 = 'log(1 + hours below 30 PSU)',
    log_hours_below26 = 'log(1 + hours below 26 PSU)'
)
base_names <- c(
    'thermal_dhw', 'cots_hazard', 'cyclone_wave_hours', 'pre_cover'
)
strata <- tibble(
    stratum = c('LTMP 9 m', 'Manta 9 m', 'MMP 2 m', 'MMP 5 m'),
    depth = c(9, 9, 2, 5),
    is_manta = c(0, 1, 0, 0),
    is_mmp = c(0, 0, 1, 1)
)

make_spline_spec <- function(x) {
    bounds <- range(x, na.rm = TRUE)
    interior <- sort(unique(x[x > bounds[[1]] & x < bounds[[2]]]))
    if (!length(interior)) stop('Metric has no interior spline support')
    knot <- median(interior)
    basis <- ns(x, knots = knot, Boundary.knots = bounds)
    center <- colMeans(basis)
    scale <- apply(basis, 2, sd)
    scale[!is.finite(scale) | scale == 0] <- 1
    list(knot = knot, bounds = bounds, center = center, scale = scale)
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

prepare_context <- function(data) {
    data <- data |>
        mutate(
            ReefID = toupper(trimws(ReefID)),
            depth = as.numeric(depth),
            is_manta = as.numeric(programme_key == 'manta'),
            is_mmp = as.numeric(programme_key == 'mmp'),
            reef_index = match(ReefID, unique(ReefID))
        )
    scaling <- tibble(variable = character(), median = numeric(),
                      mean = numeric(), sd = numeric())
    active_base <- character()
    for (variable in base_names) {
        observed <- data[[variable]][is.finite(data[[variable]])]
        fill <- if (length(observed)) median(observed) else 0
        data[[variable]][!is.finite(data[[variable]])] <- fill
        center <- mean(data[[variable]])
        spread <- sd(data[[variable]])
        if (!is.finite(spread) || spread == 0) spread <- 1
        z_name <- paste0(variable, '_z')
        data[[z_name]] <- (data[[variable]] - center) / spread
        varying <- n_distinct(data[[variable]]) > 1L
        if (varying) active_base <- c(active_base, z_name)
        scaling <- bind_rows(
            scaling,
            tibble(variable, median = fill, mean = center, sd = spread)
        )
    }
    list(data = data, scaling = scaling, active_base = active_base)
}

random_term <- paste0(
    "f(reef_index, model = 'iid', constr = TRUE, ",
    "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
)
fit_hurdle <- function(data, fixed_terms) {
    rhs <- paste(c(fixed_terms, random_term), collapse = ' + ')
    occurrence_data <- data |>
        mutate(response = as.numeric(mortality > 0))
    occurrence <- inla(
        as.formula(paste('response ~', rhs)),
        family = 'binomial', Ntrials = 1, data = occurrence_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0,
            prec.intercept = 1
        ),
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(config = TRUE, waic = TRUE, dic = TRUE),
        verbose = FALSE
    )
    positive <- data |> filter(mortality > 0)
    n_positive <- nrow(positive)
    positive$response <- (
        positive$mortality * (n_positive - 1) + 0.5
    ) / n_positive
    magnitude <- inla(
        as.formula(paste('response ~', rhs)),
        family = 'beta', data = positive,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0,
            prec.intercept = 1
        ),
        control.family = list(
            hyper = list(
                theta = list(prior = 'loggamma', param = c(2, 0.1))
            )
        ),
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(config = TRUE, waic = TRUE, dic = TRUE),
        verbose = FALSE
    )
    list(occurrence = occurrence, magnitude = magnitude)
}
fixed_draws <- function(fit, samples = 800L) {
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

all_fits <- list()
fit_summary <- tibble()
fixed_effects <- tibble()
curves <- tibble()
scaling_output <- tibble()

for (context_index in seq_along(files)) {
    context <- names(files)[[context_index]]
    prepared <- prepare_context(read_csv(files[[context]], show_col_types = FALSE))
    data <- prepared$data
    scaling_output <- bind_rows(
        scaling_output,
        prepared$scaling |> mutate(analysis_context = context)
    )
    spline_specs <- lapply(data[metric_names], make_spline_spec)
    names(spline_specs) <- metric_names
    for (metric in metric_names) {
        data <- add_metric_basis(data, metric, spline_specs[[metric]])
    }
    base_terms <- c('depth', 'is_manta', 'is_mmp', prepared$active_base)
    model_specs <- list(baseline = base_terms)
    for (metric in metric_names) {
        model_specs[[metric]] <- c(
            base_terms, paste0(metric, c('_ns1', '_ns2'))
        )
    }

    context_fits <- list()
    posterior <- list()
    for (candidate_index in seq_along(model_specs)) {
        candidate <- names(model_specs)[[candidate_index]]
        message('Fitting INLA ', context, ': ', candidate)
        fitted <- fit_hurdle(data, model_specs[[candidate]])
        context_fits[[candidate]] <- fitted
        posterior[[candidate]] <- list(
            occurrence = fixed_draws(fitted$occurrence),
            magnitude = fixed_draws(fitted$magnitude)
        )
        fixed_effects <- bind_rows(
            fixed_effects,
            as_tibble(
                fitted$occurrence$summary.fixed, rownames = 'term'
            ) |>
                mutate(
                    analysis_context = context, candidate,
                    component = 'occurrence'
                ),
            as_tibble(
                fitted$magnitude$summary.fixed, rownames = 'term'
            ) |>
                mutate(
                    analysis_context = context, candidate,
                    component = 'positive_magnitude'
                )
        )
        fit_summary <- bind_rows(
            fit_summary,
            tibble(
                analysis_context = context, candidate,
                component = c('occurrence', 'positive_magnitude'),
                rows = c(nrow(data), sum(data$mortality > 0)),
                reefs = n_distinct(data$ReefID),
                events = n_distinct(data$event_year),
                waic = c(
                    fitted$occurrence$waic$waic,
                    fitted$magnitude$waic$waic
                ),
                dic = c(
                    fitted$occurrence$dic$dic,
                    fitted$magnitude$dic$dic
                )
            )
        )
    }

    for (metric in metric_names) {
        candidate <- metric
        fixed_terms <- model_specs[[candidate]]
        formula <- reformulate(fixed_terms)
        beta_occurrence <- posterior[[candidate]]$occurrence
        beta_magnitude <- posterior[[candidate]]$magnitude
        grid <- unique(as.numeric(quantile(
            data[[metric]], probs = seq(0, 1, length.out = 41),
            na.rm = TRUE
        )))
        for (s in seq_len(nrow(strata))) {
            for (value in grid) {
                new_data <- data
                new_data[[metric]] <- value
                new_data$depth <- strata$depth[[s]]
                new_data$is_manta <- strata$is_manta[[s]]
                new_data$is_mmp <- strata$is_mmp[[s]]
                new_data <- add_metric_basis(
                    new_data, metric, spline_specs[[metric]]
                )
                design <- model.matrix(formula, new_data)
                occurrence_probability <- plogis(
                    design %*% beta_occurrence
                )
                positive_magnitude <- plogis(
                    design %*% beta_magnitude
                )
                standardized_draw <- colMeans(
                    occurrence_probability * positive_magnitude
                )
                curves <- bind_rows(
                    curves,
                    tibble(
                        analysis_context = context,
                        metric, metric_value = value,
                        stratum = strata$stratum[[s]],
                        predicted_mortality = mean(standardized_draw),
                        lower = quantile(standardized_draw, 0.025),
                        upper = quantile(standardized_draw, 0.975),
                        framework = 'INLA two-part additive spline'
                    )
                )
            }
        }
    }
    all_fits[[context]] <- list(
        fits = context_fits, spline_specs = spline_specs,
        model_specs = model_specs
    )
}

fit_comparison <- fit_summary |>
    group_by(analysis_context, candidate) |>
    summarise(
        combined_waic = sum(waic), combined_dic = sum(dic),
        .groups = 'drop'
    ) |>
    group_by(analysis_context) |>
    mutate(
        delta_waic_from_baseline =
            combined_waic - combined_waic[candidate == 'baseline'],
        delta_dic_from_baseline =
            combined_dic - combined_dic[candidate == 'baseline']
    ) |>
    ungroup()
curves <- curves |>
    mutate(metric_label = recode(metric, !!!metric_labels))

saveRDS(
    all_fits,
    file.path(output_dir, 'expanded_salinity_inla_models.rds')
)
write_csv(
    scaling_output,
    file.path(output_dir, 'expanded_salinity_inla_preprocessing.csv')
)
write_csv(
    fixed_effects,
    file.path(output_dir, 'expanded_salinity_inla_fixed_effects.csv')
)
write_csv(
    fit_summary,
    file.path(output_dir, 'expanded_salinity_inla_fit_summary.csv')
)
write_csv(
    fit_comparison,
    file.path(output_dir, 'expanded_salinity_inla_fit_comparison.csv')
)
write_csv(
    curves,
    file.path(output_dir, 'expanded_salinity_inla_partial_dependence.csv')
)

brt_curves <- read_csv(
    file.path(output_dir, 'expanded_salinity_brt_partial_dependence.csv'),
    show_col_types = FALSE
) |>
    transmute(
        analysis_context, metric, metric_value, stratum,
        predicted_mortality, lower, upper,
        framework = 'BRT direct', metric_label
    )
comparison_curves <- bind_rows(brt_curves, curves) |>
    mutate(
        framework = factor(
            framework,
            levels = c('BRT direct', 'INLA two-part additive spline')
        )
    )
comparison_plot <- ggplot(
    comparison_curves,
    aes(
        metric_value, predicted_mortality,
        colour = stratum, fill = stratum, linetype = framework
    )
) +
    geom_ribbon(
        aes(ymin = lower, ymax = upper, linetype = NULL),
        alpha = 0.04, colour = NA
    ) +
    geom_line(linewidth = 0.75) +
    facet_grid(analysis_context ~ metric_label, scales = 'free_x') +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(
        x = NULL, y = 'Marginal expected relative loss',
        colour = 'Programme/depth', fill = 'Programme/depth',
        linetype = 'Framework',
        title = 'Expanded salinity relationships: BRT and additive INLA',
        subtitle = paste(
            'Each metric is fit separately after heat, COTS and cyclone;',
            'INLA bands are fixed-effect posterior intervals'
        )
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'bottom')
ggsave(
    file.path(
        output_dir, 'expanded_salinity_brt_inla_partial_dependence.png'
    ),
    comparison_plot, width = 15, height = 8, dpi = 180
)

cat('Expanded salinity INLA diagnostic complete.\n')
print(fit_comparison)
