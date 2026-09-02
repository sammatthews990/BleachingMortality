# Full-data retrospective model: event-specific deviations in each DHW basis.
# It is explanatory-only and never contributes to operational selection.

suppressPackageStartupMessages({
    library(dplyr)
    library(INLA)
    library(readr)
    library(gbm)
    library(ggplot2)
})
source('src/models/fit_inla_spatiotemporal_screen.R')
source('src/lib/model_diagnostics.R')

root <- project_root()
output_dir <- 'output/explanatory_event_dhw'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
candidate <- candidates |>
    filter(candidate ==
        'persistent_rw1_local_dhw_decomposed_hazards_freshwater_partial_pool')
stopifnot(nrow(candidate) == 1L)

base_make_fixed_rows <- make_fixed_rows
make_fixed_rows <- function(rows, component, observed_event_keys,
                            observed_event_indices) {
    base_make_fixed_rows(
        rows, component, observed_event_keys, observed_event_indices
    ) |>
        mutate(
            event_dhw_linear_index = event_dhw_index,
            event_dhw_linear_weight = ann_maxdhw_z,
            event_dhw_excess4_index = event_dhw_index,
            event_dhw_excess4_weight = dhw_excess4_z,
            event_dhw_excess8_index = event_dhw_index,
            event_dhw_excess8_weight = dhw_excess8_z
        )
}

make_formula <- function(candidate) {
    rhs <- c(
        paste0('layer_', seq_along(families)),
        candidate_shared_terms(candidate),
        'f(reef_event_index, model=\'iid\', constr=TRUE, hyper=pc_prec)',
        'f(event_effect_index, model=\'iid\', constr=TRUE, hyper=pc_prec)',
        'f(spatial_field, model=spde)',
        paste0(
            'f(event_dhw_linear_index, event_dhw_linear_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_prec)'
        ),
        paste0(
            'f(event_dhw_excess4_index, event_dhw_excess4_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_prec)'
        ),
        paste0(
            'f(event_dhw_excess8_index, event_dhw_excess8_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_prec)'
        ),
        paste0(
            'f(pool_acropora_index, pool_acropora_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_pool_prec)'
        ),
        paste0(
            'f(pool_novelty_index, pool_novelty_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_pool_prec)'
        ),
        paste0(
            'f(pool_freshwater_index, pool_freshwater_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_pool_prec)'
        ),
        paste0(
            'f(pool_cloud_index, pool_cloud_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_pool_prec)'
        ),
        paste0(
            'f(pool_current_index, pool_current_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_pool_prec)'
        ),
        paste0(
            'f(pool_cyclone_index, pool_cyclone_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_pool_prec)'
        ),
        paste0(
            'f(pool_cots_index, pool_cots_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_pool_prec)'
        ),
        paste0(
            'f(pool_freshwater_expanded_index, ',
            'pool_freshwater_expanded_weight, model=\'iid\', ',
            'constr=TRUE, hyper=pc_pool_prec)'
        )
    )
    as.formula(paste('response ~ -1 +', paste(rhs, collapse = ' + ')))
}

fit_path <- file.path(output_dir, 'full_fit.rds')
if (file.exists(fit_path)) {
    result <- readRDS(fit_path)
} else {
    result <- fit_candidate(
        data, data[1, , drop = FALSE], candidate,
        compute_criteria = TRUE, seed = 20260831L
    )
    saveRDS(result, fit_path)
}

corrected_rows <- apply_candidate_dhw(data, candidate)
prepared_rows <- prepare_fold(corrected_rows, corrected_rows)$analysis
fixed <- result$fit$summary.fixed |>
    as_tibble(rownames = 'term')
write_csv(fixed, file.path(output_dir, 'fixed_effects.csv'))
write_csv(
    result$fit$summary.hyperpar |>
        as_tibble(rownames = 'term'),
    file.path(output_dir, 'hyperparameters.csv')
)
write_csv(
    tibble(
        model_id = 'explanatory_event_dhw_interaction_full_data',
        scope = 'full_data_explanatory_only',
        waic = result$fit$waic$waic,
        dic = result$fit$dic$dic,
        mean_neg_log_cpo = mean(-log(
            result$fit$cpo$cpo[is.finite(result$fit$cpo$cpo)]
        ))
    ),
    file.path(output_dir, 'criteria.csv')
)
write_csv(
    prepared_rows |>
        group_by(event_year) |>
        summarise(
            observations = n(), reefs = n_distinct(ReefID),
            dhw_min = min(ann_maxdhw), dhw_max = max(ann_maxdhw),
            severe_losses = sum(mortality_prop >= 0.50), .groups = 'drop'
        ),
    file.path(output_dir, 'event_support.csv')
)

event_terms <- bind_rows(
    result$fit$summary.random$event_dhw_linear_index |>
        as_tibble() |>
        mutate(event_year = event_years[ID], thermal_basis = 'DHW'),
    result$fit$summary.random$event_dhw_excess4_index |>
        as_tibble() |>
        mutate(event_year = event_years[ID], thermal_basis = 'DHW_above_4'),
    result$fit$summary.random$event_dhw_excess8_index |>
        as_tibble() |>
        mutate(event_year = event_years[ID], thermal_basis = 'DHW_above_8')
) |>
    mutate(
        lower95 = .data[['0.025quant']],
        median = .data[['0.5quant']],
        upper95 = .data[['0.975quant']]
    ) |>
    select(event_year, thermal_basis, mean, sd, lower95, median, upper95)
write_csv(event_terms, file.path(output_dir, 'event_thermal_deviations.csv'))

coefficient <- function(term) {
    row <- fixed |> filter(.data$term == .env$term)
    if (nrow(row) != 1L) stop('Missing fixed effect: ', term)
    row
}
event_component <- function(name, year) {
    row <- event_terms |>
        filter(thermal_basis == .env$name, event_year == .env$year)
    if (nrow(row) != 1L) stop('Missing event effect')
    row
}
d_grid <- seq(0, 20, length.out = 161)
d_mean <- mean(prepared_rows$ann_maxdhw)
d_sd <- sd(prepared_rows$ann_maxdhw)
h4_mean <- mean(prepared_rows$dhw_excess4)
h4_sd <- sd(prepared_rows$dhw_excess4)
h8_mean <- mean(prepared_rows$dhw_excess8)
h8_sd <- sd(prepared_rows$dhw_excess8)
base_d <- coefficient('ann_maxdhw_z')
base_h4 <- coefficient('dhw_excess4_z')
base_h8 <- coefficient('dhw_excess8_z')
layer_occ <- coefficient('layer_1')
layer_mag <- coefficient('layer_4')

set.seed(20260831L)
n_draws <- 2000L
draw_response <- function(year) {
    event_d <- event_component('DHW', year)
    event_h4 <- event_component('DHW_above_4', year)
    event_h8 <- event_component('DHW_above_8', year)
    draws <- tibble(
        b_d = rnorm(n_draws, base_d$mean + event_d$mean,
                    sqrt(base_d$sd^2 + event_d$sd^2)),
        b_h4 = rnorm(n_draws, base_h4$mean + event_h4$mean,
                     sqrt(base_h4$sd^2 + event_h4$sd^2)),
        b_h8 = rnorm(n_draws, base_h8$mean + event_h8$mean,
                     sqrt(base_h8$sd^2 + event_h8$sd^2)),
        occ = rnorm(n_draws, layer_occ$mean, layer_occ$sd),
        mag = rnorm(n_draws, layer_mag$mean, layer_mag$sd)
    )
    purrr::map_dfr(d_grid, function(dhw) {
        dz <- (dhw - d_mean) / d_sd
        h4z <- (max(dhw - 4, 0) - h4_mean) / h4_sd
        h8z <- (max(dhw - 8, 0) - h8_mean) / h8_sd
        response <- plogis(draws$occ + draws$b_d * dz +
            draws$b_h4 * h4z + draws$b_h8 * h8z) *
            plogis(draws$mag + draws$b_d * dz +
                draws$b_h4 * h4z + draws$b_h8 * h8z)
        tibble(
            event_year = year, dhw = dhw,
            median = median(response),
            lower95 = quantile(response, 0.025),
            upper95 = quantile(response, 0.975)
        )
    })
}
event_curve <- purrr::map_dfr(event_years, draw_response)
write_csv(event_curve, file.path(output_dir, 'dose_response_event.csv'))
event_plot <- ggplot(
    event_curve,
    aes(dhw, median, colour = factor(event_year), fill = factor(event_year))
) +
    geom_ribbon(aes(ymin = lower95, ymax = upper95), alpha = 0.12,
                colour = NA) +
    geom_line(linewidth = 1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
        x = 'Local-first annual maximum DHW',
        y = 'Population mean relative mortality',
        colour = 'Event', fill = 'Event',
        title = 'Explanatory full-data event-specific DHW response',
        subtitle = 'Marginal posterior approximation; explanatory only'
    ) +
    theme_bw(base_size = 11)
ggsave(file.path(output_dir, 'event_dhw_response.png'), event_plot,
       width = 10, height = 6.5, dpi = 300)
ggsave(file.path(output_dir, 'event_dhw_response.pdf'), event_plot,
       width = 10, height = 6.5)
save_figure_bundle(
    event_plot, 'Fig-INLA-04_event_specific_dhw_response', event_curve,
    'Event-specific DHW response curves from the full-data explanatory INLA model.',
    'Shows how the fitted thermal response differs among observed bleaching events while holding the other fixed effects at their reference values.',
    'These curves are explanatory only: event effects are estimated using all data and do not represent prospective forecast skill.',
    'explanatory_event_dhw_interaction_full_data', 'INLA', 'event_specific_dose_response',
    'full_data_explanatory', root, TRUE,
    code_source = 'src/models/fit_inla_explanatory_event_dhw.R'
)

brt_rows <- prepared_rows |> mutate(event_factor = factor(event_year))
write_csv(brt_rows, file.path(output_dir, 'event_dhw_brt_data.csv'))
brt_predictors <- c(candidate_shared_terms(candidate), 'event_factor')
brt_formula <- as.formula(paste(
    'mortality_prop ~ event_factor * (ann_maxdhw_z + dhw_excess4_z + dhw_excess8_z)',
    '+', paste(setdiff(brt_predictors, c(
        'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z', 'event_factor'
    )), collapse = ' + ')
))
set.seed(20260831L)
brt_fit <- gbm(
    brt_formula, data = brt_rows, distribution = 'gaussian',
    n.trees = 1200, interaction.depth = 2, shrinkage = 0.01,
    n.minobsinnode = 12, bag.fraction = 0.7, train.fraction = 1,
    keep.data = TRUE, verbose = FALSE
)
saveRDS(brt_fit, file.path(output_dir, 'event_dhw_brt.rds'))
brt_influence <- summary(brt_fit, plotit = FALSE) |>
    as_tibble() |>
    rename(variable = var, relative_influence = rel.inf)
write_csv(brt_influence, file.path(output_dir, 'event_dhw_brt_influence.csv'))

brt_importance_plot <- brt_influence |>
    arrange(relative_influence) |>
    mutate(variable = factor(variable, levels = variable)) |>
    ggplot(aes(relative_influence, variable)) +
    geom_col(fill = '#2a9d8f', alpha = 0.9) +
    geom_text(aes(label = sprintf('%.1f%%', relative_influence)), hjust = -0.1, size = 3.2) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(
        title = 'Event-interaction BRT: variable importance',
        subtitle = 'Full-data explanatory model; relative influence is not a causal effect size',
        x = 'Relative influence (%)', y = NULL
    ) +
    theme_minimal(base_size = 12)
save_figure_bundle(
    brt_importance_plot, 'Fig-BRT-01_event_dhw_variable_importance', brt_influence,
    'Relative influence from the full-data explanatory BRT with event-by-thermal terms.',
    'Ranks predictor contribution within this fitted BRT; it is descriptive and not used for operational model selection.',
    'Relative influence is not a causal effect size and can be shared among correlated predictors.',
    'explanatory_event_dhw_interaction_full_data', 'BRT', 'variable_importance',
    'full_data_explanatory', root, TRUE,
    code_source = 'src/models/fit_inla_explanatory_event_dhw.R'
)

numeric_terms <- brt_influence |>
    arrange(desc(relative_influence)) |>
    pull(variable) |>
    intersect(names(brt_rows)) |>
    keep(~ is.numeric(brt_rows[[.x]])) |>
    head(6)
pdp_data <- map_dfr(numeric_terms, function(variable) {
    variable_index <- match(variable, brt_fit$var.names)
    grid <- gbm::plot.gbm(brt_fit, i.var = variable_index, return.grid = TRUE)
    tibble(
        variable = variable,
        predictor_value = grid[[1]],
        partial_dependence = grid[[2]]
    )
})
write_csv(pdp_data, file.path(output_dir, 'event_dhw_brt_partial_dependence.csv'))
brt_pdp_plot <- pdp_data |>
    ggplot(aes(predictor_value, partial_dependence)) +
    geom_line(linewidth = 0.75, colour = '#0072B2') +
    facet_wrap(~ variable, scales = 'free_x', ncol = 3) +
    labs(
        title = 'Event-interaction BRT: one-dimensional partial dependence',
        subtitle = 'Full-data explanatory model; other predictors are averaged over observed rows',
        x = 'Predictor value', y = 'Partial dependence (link scale)'
    ) +
    theme_minimal(base_size = 12)
save_figure_bundle(
    brt_pdp_plot, 'Fig-BRT-02_event_dhw_partial_dependence', pdp_data,
    'One-dimensional BRT partial-dependence curves for the most influential numeric predictors.',
    'Describes marginal patterns after averaging the other predictors across observed rows.',
    'Partial-dependence curves can be unreliable in sparsely supported predictor combinations; read with the event-specific INLA dose-response curves.',
    'explanatory_event_dhw_interaction_full_data', 'BRT', 'partial_dependence',
    'full_data_explanatory', root, TRUE,
    code_source = 'src/models/fit_inla_explanatory_event_dhw.R'
)
