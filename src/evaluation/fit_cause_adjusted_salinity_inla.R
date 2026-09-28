# Smooth cause-adjusted INLA diagnostic for the 2024 surface-salinity screen.
#
# This is an explanatory full-data fit, not a model-selection candidate. It uses
# the same equal-weight programme-depth outcomes as the cause-adjusted BRT and
# estimates population-level marginal curves from a Bernoulli-beta hurdle
# model. Reef random intercepts account for repeated programme/depth outcomes.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(INLA)
    library(readr)
    library(splines)
    library(tidyr)
})

set.seed(20260916L)

data_file <- paste0(
    'output/surface_salinity_mortality/',
    'cause_adjusted_salinity_brt_data.csv'
)
source_file <- 'output/explanatory_event_dhw/event_dhw_brt_data.csv'
brt_curve_file <- paste0(
    'output/surface_salinity_mortality/',
    'cause_adjusted_salinity_brt_partial_dependence.csv'
)
output_dir <- 'output/surface_salinity_mortality'

if (!all(file.exists(c(data_file, source_file, brt_curve_file)))) {
    stop('Missing cause-adjusted BRT data, source rows or BRT curves')
}

data <- read_csv(data_file, show_col_types = FALSE) |>
    mutate(
        stratum = case_when(
            programme_key == 'mmp' & depth == 2 ~ 'MMP 2 m',
            programme_key == 'mmp' & depth == 5 ~ 'MMP 5 m',
            programme_key == 'manta' ~ 'Manta 9 m',
            TRUE ~ 'LTMP 9 m'
        ),
        stratum = factor(
            stratum,
            levels = c('LTMP 9 m', 'Manta 9 m', 'MMP 2 m', 'MMP 5 m')
        ),
        stratum_manta9 = as.numeric(stratum == 'Manta 9 m'),
        stratum_mmp2 = as.numeric(stratum == 'MMP 2 m'),
        stratum_mmp5 = as.numeric(stratum == 'MMP 5 m'),
        across(
            c(bleaching_label, flood_label, cyclone_label, cots_label),
            as.numeric
        )
    )

if (nrow(data) != 178L || n_distinct(data$ReefID) != 102L) {
    stop('Expected 178 programme-depth outcomes from 102 reefs')
}

# Audit the BRT upturn at the high-salinity end. The 35.15-PSU boundary is the
# point at which the displayed BRT curve begins its final rise, not a proposed
# ecological threshold.
high_sss_boundary <- 35.15
source_rows <- read_csv(source_file, show_col_types = FALSE) |>
    mutate(ReefID = toupper(trimws(ReefID))) |>
    filter(event_year == 2024L)

reef_geography <- source_rows |>
    group_by(ReefID) |>
    summarise(
        Region = first(na.omit(Region), default = NA_character_),
        SECTOR = first(na.omit(SECTOR), default = NA_character_),
        .groups = 'drop'
    ) |>
    mutate(
        Region = case_when(
            !is.na(Region) ~ Region,
            SECTOR %in% c('CB', 'PO') ~ 'Southern GBR',
            TRUE ~ Region
        )
    )

source_audit <- source_rows |>
    group_by(ReefID, programme_key, depth) |>
    summarise(
        disturbance_text = paste(
            unique(na.omit(disturbance_text[disturbance_text != ''])),
            collapse = ' | '
        ),
        .groups = 'drop'
    )

high_sss_audit <- data |>
    filter(sss_min >= high_sss_boundary) |>
    left_join(reef_geography, by = 'ReefID', relationship = 'many-to-one') |>
    left_join(
        source_audit,
        by = c('ReefID', 'programme_key', 'depth'),
        relationship = 'many-to-one'
    ) |>
    arrange(desc(sss_min), ReefID, programme_key, depth) |>
    select(
        ReefID, ReefName, Region, SECTOR, programme_key, depth, mortality,
        sss_min, hours_below30, local_first_dhw, cots_hazard,
        cyclone_wave_hours, bleaching_label, flood_label, cyclone_label,
        cots_label, disturbance_text
    )

write_csv(
    high_sss_audit,
    file.path(output_dir, 'high_sss_uptick_reef_audit.csv'),
    na = ''
)

salinity_band_summary <- data |>
    mutate(
        salinity_band = case_when(
            sss_min >= high_sss_boundary ~ 'high_branch_ge_35.15',
            sss_min >= 34.5 ~ 'shoulder_34.5_to_35.15',
            TRUE ~ 'lower_than_34.5'
        )
    ) |>
    group_by(salinity_band) |>
    summarise(
        rows = n(),
        reefs = n_distinct(ReefID),
        mean_mortality = mean(mortality),
        median_mortality = median(mortality),
        mean_local_first_dhw = mean(local_first_dhw),
        median_local_first_dhw = median(local_first_dhw),
        bleaching_fraction = mean(bleaching_label),
        flood_fraction = mean(flood_label),
        cyclone_fraction = mean(cyclone_label),
        cots_fraction = mean(cots_label),
        mean_cots_hazard = mean(cots_hazard),
        mean_cyclone_wave_hours = mean(cyclone_wave_hours),
        .groups = 'drop'
    )

write_csv(
    salinity_band_summary,
    file.path(output_dir, 'high_sss_uptick_band_summary.csv'),
    na = ''
)

metric_names <- c(
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
stratum_terms <- c('stratum_manta9', 'stratum_mmp2', 'stratum_mmp5')
label_terms <- c(
    'bleaching_label', 'flood_label', 'cyclone_label', 'cots_label'
)

make_spline_spec <- function(x) {
    bounds <- range(x, na.rm = TRUE)
    interior <- sort(unique(x[x > bounds[[1]] & x < bounds[[2]]]))
    if (!length(interior)) stop('Metric has no interior support')
    knot <- median(interior)
    basis <- ns(x, knots = knot, Boundary.knots = bounds)
    center <- colMeans(basis)
    scale <- apply(basis, 2, sd)
    scale[!is.finite(scale) | scale == 0] <- 1
    list(
        knot = knot, bounds = bounds, center = center, scale = scale
    )
}

spline_specs <- lapply(data[metric_names], make_spline_spec)

add_spline_terms <- function(x) {
    for (metric in metric_names) {
        spec <- spline_specs[[metric]]
        basis <- ns(
            x[[metric]], knots = spec$knot,
            Boundary.knots = spec$bounds
        )
        basis <- sweep(basis, 2, spec$center, '-')
        basis <- sweep(basis, 2, spec$scale, '/')
        x[[paste0(metric, '_ns1')]] <- basis[, 1]
        x[[paste0(metric, '_ns2')]] <- basis[, 2]
    }
    x
}

data <- add_spline_terms(data)
spline_terms <- unlist(lapply(
    metric_names, function(x) paste0(x, c('_ns1', '_ns2'))
))
data$reef_index <- match(data$ReefID, unique(data$ReefID))

random_term <- paste0(
    "f(reef_index, model = 'iid', constr = TRUE, ",
    "hyper = list(prec = list(prior = 'pc.prec', param = c(1, 0.01))))"
)

model_specs <- list(
    cause_exposures = c(stratum_terms, spline_terms),
    cause_exposures_and_labels = c(
        stratum_terms, spline_terms, label_terms
    )
)

fit_hurdle <- function(x, fixed_terms) {
    rhs <- paste(c(fixed_terms, random_term), collapse = ' + ')
    occurrence_data <- x |>
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
    positive <- x |> filter(mortality > 0)
    positive_n <- nrow(positive)
    positive$response <- (
        positive$mortality * (positive_n - 1) + 0.5
    ) / positive_n
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

fits <- list()
posterior_fixed <- list()
fixed_effects <- tibble()
fit_summary <- tibble()
for (context in names(model_specs)) {
    message('Fitting cause-adjusted INLA: ', context)
    fitted <- fit_hurdle(data, model_specs[[context]])
    fits[[context]] <- fitted
    posterior_fixed[[context]] <- list(
        occurrence = fixed_draws(fitted$occurrence),
        magnitude = fixed_draws(fitted$magnitude)
    )
    fixed_effects <- bind_rows(
        fixed_effects,
        as_tibble(fitted$occurrence$summary.fixed, rownames = 'term') |>
            mutate(model_context = context, component = 'occurrence'),
        as_tibble(fitted$magnitude$summary.fixed, rownames = 'term') |>
            mutate(model_context = context, component = 'magnitude')
    )
    fit_summary <- bind_rows(
        fit_summary,
        tibble(
            model_context = context,
            component = c('occurrence', 'magnitude'),
            rows = c(nrow(data), sum(data$mortality > 0)),
            reefs = n_distinct(data$ReefID),
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

saveRDS(
    list(fits = fits, spline_specs = spline_specs, model_specs = model_specs),
    file.path(output_dir, 'cause_adjusted_salinity_inla_models.rds')
)
write_csv(
    fixed_effects,
    file.path(output_dir, 'cause_adjusted_salinity_inla_fixed_effects.csv'),
    na = ''
)
write_csv(
    fit_summary,
    file.path(output_dir, 'cause_adjusted_salinity_inla_fit_summary.csv'),
    na = ''
)

brt_curves <- read_csv(brt_curve_file, show_col_types = FALSE)
metric_grids <- brt_curves |>
    distinct(metric, metric_value) |>
    group_by(metric) |>
    summarise(grid = list(sort(unique(metric_value))), .groups = 'drop')

strata <- tibble(
    programme_depth = levels(data$stratum),
    stratum_manta9 = c(0, 1, 0, 0),
    stratum_mmp2 = c(0, 0, 1, 0),
    stratum_mmp5 = c(0, 0, 0, 1)
)

marginal_curves <- tibble()
for (context in names(model_specs)) {
    fixed_formula <- reformulate(model_specs[[context]])
    beta_occurrence <- posterior_fixed[[context]]$occurrence
    beta_magnitude <- posterior_fixed[[context]]$magnitude
    for (metric in metric_names) {
        grid <- metric_grids$grid[[match(metric, metric_grids$metric)]]
        for (stratum_index in seq_len(nrow(strata))) {
            for (metric_value in grid) {
                new_data <- data
                new_data[[metric]] <- metric_value
                for (term in stratum_terms) {
                    new_data[[term]] <- strata[[term]][[stratum_index]]
                }
                new_data <- add_spline_terms(new_data)
                design <- model.matrix(fixed_formula, new_data)
                occurrence_probability <- plogis(design %*% beta_occurrence)
                positive_magnitude <- plogis(design %*% beta_magnitude)
                standardized_draw <- colMeans(
                    occurrence_probability * positive_magnitude
                )
                marginal_curves <- bind_rows(
                    marginal_curves,
                    tibble(
                        model_context = context,
                        metric = metric,
                        metric_value = metric_value,
                        programme_depth =
                            strata$programme_depth[[stratum_index]],
                        predicted_mortality = mean(standardized_draw),
                        lower = quantile(standardized_draw, 0.025),
                        upper = quantile(standardized_draw, 0.975),
                        framework = 'INLA two-part additive spline'
                    )
                )
            }
        }
    }
}

marginal_curves <- marginal_curves |>
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
    marginal_curves,
    file.path(
        output_dir,
        'cause_adjusted_salinity_inla_partial_dependence.csv'
    ),
    na = ''
)

comparison_curves <- bind_rows(
    brt_curves |>
        mutate(framework = 'BRT direct'),
    marginal_curves
) |>
    filter(metric %in% c('sss_min', 'hours_below30')) |>
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
        colour = programme_depth, fill = programme_depth,
        linetype = framework
    )
) +
    geom_ribbon(
        aes(ymin = lower, ymax = upper, linetype = NULL),
        alpha = 0.045, colour = NA
    ) +
    geom_line(linewidth = 0.8) +
    facet_grid(
        model_context_label ~ metric_label,
        scales = 'free_x'
    ) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
        x = NULL,
        y = 'Marginal expected relative mortality',
        colour = 'Programme/depth', fill = 'Programme/depth',
        linetype = 'Framework',
        title = 'BRT and additive INLA salinity relationships',
        subtitle = paste(
            'INLA uses common two-df exposure splines and a reef random effect;',
            'curves are standardized over the 178 observed outcomes'
        )
    ) +
    theme_bw(base_size = 10) +
    theme(
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 35, hjust = 1)
    )

ggsave(
    file.path(
        output_dir,
        'cause_adjusted_salinity_brt_inla_partial_dependence.png'
    ),
    comparison_plot, width = 13, height = 8, dpi = 180
)

cat(
    'Wrote high-SSS audit and cause-adjusted INLA curves using',
    nrow(data), 'programme-depth outcomes from',
    n_distinct(data$ReefID), 'reefs.\n'
)
print(salinity_band_summary)

