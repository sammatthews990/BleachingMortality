# Reusable visual diagnostics for the model ledger.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(tidyr)
    library(purrr)
})
if (!exists('project_root', mode = 'function')) {
    source('scripts/model_registry.R')
}

figure_root <- function(root = project_root()) {
    path <- file.path(root, 'output', 'fig')
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

save_figure_bundle <- function(plot, figure_id, data, caption, interpretation,
                               caveats, model_id, framework, figure_type,
                               analysis_role, root = project_root(),
                               manuscript_candidate = FALSE,
                               version = 'v001',
                               code_source = 'scripts/model_diagnostics.R',
                               width = 10, height = 6.5) {
    dir <- figure_root(root)
    stem <- paste0(figure_id, '_', version)
    png_path <- file.path(dir, paste0(stem, '.png'))
    pdf_path <- file.path(dir, paste0(stem, '.pdf'))
    csv_path <- file.path(dir, paste0(stem, '_data.csv'))
    ggsave(png_path, plot, width = width, height = height, dpi = 300)
    ggsave(pdf_path, plot, width = width, height = height)
    write_csv(data, csv_path)

    manifest_path <- file.path(dir, 'manifest.csv')
    record <- tibble(
        figure_id = figure_id,
        relative_path = file.path('output', 'fig', basename(png_path)),
        model_id = model_id,
        framework = framework,
        figure_type = figure_type,
        analysis_role = analysis_role,
        selection_status = if_else(
            manuscript_candidate, 'manuscript_candidate', 'diagnostic'
        ),
        data_scope = 'relative bleaching mortality',
        validation_scheme = 'see model registry',
        created_utc = format(
            as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'
        ),
        code_source = code_source,
        input_hash = unname(tools::md5sum(csv_path)),
        caption = caption,
        interpretation = interpretation,
        caveats = caveats,
        manuscript_candidate = manuscript_candidate,
        version = version
    )
    existing <- if (file.exists(manifest_path)) {
        read_csv(manifest_path, show_col_types = FALSE) |>
            mutate(created_utc = as.character(created_utc))
    } else tibble()
    existing |>
        filter(!(figure_id == record$figure_id & version == record$version)) |>
        bind_rows(record) |>
        write_csv(manifest_path)
    invisible(list(png = png_path, pdf = pdf_path, data = csv_path))
}

write_figure_readme <- function(root = project_root()) {
    manifest_path <- file.path(figure_root(root), 'manifest.csv')
    if (!file.exists(manifest_path)) return(invisible(NULL))
    manifest <- read_csv(manifest_path, show_col_types = FALSE)
    lines <- c('# Figure registry', '')
    for (i in seq_len(nrow(manifest))) {
        row <- manifest[i, ]
        lines <- c(
            lines,
            paste0('## ', row$figure_id, ' — ', row$figure_type),
            '',
            paste0('Caption: ', row$caption),
            paste0('Interpretation: ', row$interpretation),
            paste0('Caveat: ', row$caveats),
            paste0('Model: ', row$model_id, ' (', row$framework, ').'),
            paste0('File: ', row$relative_path),
            ''
        )
    }
    writeLines(lines, file.path(figure_root(root), 'README.md'))
}

plot_fixed_effect_posteriors <- function(model_id, root = project_root()) {
    is_composite <- model_id %in% c(
        'cause_aware_exposure_gated_competing_hazards',
        'baseline_cots_logistic20_5_cyclone',
        'cots_raw_interval_logistic20_5_cyclone',
        'prospective_cots_log_peak_recency',
        'operational_rrn_raw_plus_manta_state'
    )
    if (is_composite) {
        thermal <- read_csv(
            file.path(root, 'output', 'cause_aware_competing_hazards',
                      'thermal_component_fixed_effects.csv'),
            show_col_types = FALSE
        ) |>
            mutate(component = 'thermal/freshwater')
        cause_path <- if (model_id %in% c(
            'prospective_cots_log_peak_recency',
            'operational_rrn_raw_plus_manta_state'
        )) file.path(root, 'output', 'prospective_cots_nowcast',
                     'selected_cause_fixed_effects.csv') else if (identical(
            model_id, 'cots_raw_interval_logistic20_5_cyclone'
        )) file.path(root, 'output', 'cots_raw_enso_occurrence',
                     'selected_cause_fixed_effects.csv') else
            file.path(root, 'output', 'cause_aware_competing_hazards',
                      'cause_fixed_effects.csv')
        causes <- read_csv(
            cause_path,
            show_col_types = FALSE
        ) |>
            filter(term != '(Intercept)') |>
            mutate(term = paste(component, term, sep = ': '))
        fixed <- bind_rows(thermal, causes) |>
            filter(!grepl('^layer_', term)) |>
            mutate(
                lower95 = .data[['0.025quant']],
                lower80 = mean - qnorm(0.90) * sd,
                upper80 = mean + qnorm(0.90) * sd,
                upper95 = .data[['0.975quant']],
                label = term,
                mechanism = component,
                panel = if_else(
                    component == 'thermal/freshwater',
                    'Thermal/freshwater fixed effects', 'Cause hazards'
                )
            ) |>
            select(term, mean, lower95, lower80, upper80, upper95,
                   label, mechanism, panel)
    } else {
        fixed <- read_csv(
            file.path(root, 'output', 'inla_spatiotemporal', 'full_fixed_effects.csv'),
            show_col_types = FALSE
        ) |>
            filter(candidate == .env$model_id) |>
        mutate(
            lower95 = .data[['0.025quant']],
            lower80 = mean - qnorm(0.90) * sd,
            upper80 = mean + qnorm(0.90) * sd,
            upper95 = .data[['0.975quant']]
        ) |>
            select(term, mean, lower95, lower80, upper80, upper95)
        terms <- read_model_terms(root)
        fixed <- fixed |>
            left_join(
                terms |>
                    transmute(term = id, label, mechanism = class),
                by = 'term'
            ) |>
            mutate(
                label = coalesce(label, term),
                mechanism = coalesce(mechanism, 'interaction_or_observation'),
                panel = case_when(
                    grepl('^dhw_x_', term) ~ 'Thermal interactions',
                    grepl('^layer_', term) ~ 'Observation layers',
                    TRUE ~ 'Fixed effects'
                )
            )
    }
    plot <- fixed |>
        arrange(mean) |>
        mutate(label = factor(label, levels = label)) |>
        ggplot(aes(mean, label, colour = mechanism)) +
        geom_vline(xintercept = 0, linetype = 2, colour = 'grey45') +
        geom_errorbar(aes(xmin = lower95, xmax = upper95), width = 0,
                      orientation = 'y', alpha = 0.35) +
        geom_errorbar(aes(xmin = lower80, xmax = upper80), width = 0,
                      orientation = 'y', linewidth = 1) +
        geom_point(size = 2) +
        facet_wrap(~ panel, scales = 'free_y') +
        labs(
            x = 'Posterior estimate on link scale', y = NULL,
            title = 'INLA fixed-effect posterior intervals',
            subtitle = paste('Model:', model_id)
        ) +
        theme_bw(base_size = 11) +
        theme(legend.position = 'bottom')
    save_figure_bundle(
        plot, 'Fig-INLA-01_fixed_effect_posteriors', fixed,
        'Posterior means with 80% and 95% credible intervals for selected INLA fixed effects.',
        'Shows direction and uncertainty of population-level link-scale effects.',
        'Coefficient magnitudes are not directly comparable across differently transformed predictors.',
        model_id, 'INLA', 'fixed_effect_posterior', 'operational_diagnostic',
        root, TRUE
    )
    list(plot = plot, data = fixed)
}

plot_cv_calibration <- function(model_id, root = project_root()) {
    prediction_path <- if (model_id %in% c(
        'prospective_cots_log_peak_recency',
        'operational_rrn_raw_plus_manta_state'
    )) file.path(root, 'output', 'prospective_cots_nowcast',
                 'cv_predictions.csv') else if (identical(
        model_id, 'cots_raw_interval_logistic20_5_cyclone'
    )) file.path(root, 'output', 'cots_raw_enso_occurrence',
                 'cots_cv_predictions.csv') else if (identical(
        model_id, 'baseline_cots_logistic20_5_cyclone'
    )) file.path(root, 'output', 'cots_timing_cover_cyclone_soft_gate',
                 'cv_predictions.csv') else if (identical(
        model_id, 'cause_aware_exposure_gated_competing_hazards'
    )) file.path(root, 'output', 'cause_aware_competing_hazards',
                 'gated_cv_predictions.csv') else
        file.path(root, 'output', 'inla_spatiotemporal', 'cv_predictions.csv')
    prediction_id <- if (identical(
        model_id, 'cots_raw_interval_logistic20_5_cyclone'
    )) 'cots_raw_interval_relative' else model_id
    rows <- read_csv(prediction_path, show_col_types = FALSE) |>
        filter(candidate == .env$prediction_id, scheme == 'leave_one_event_out')
    plot <- ggplot(
        rows,
        aes(predicted_mortality, observed_mortality, colour = factor(event_year))
    ) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, colour = 'grey35') +
        geom_point(alpha = 0.55, size = 1.8) +
        facet_wrap(~ programme_key) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
        labs(
            x = 'Held-out predicted mortality', y = 'Observed mortality',
            colour = 'Event',
            title = 'Leave-one-event-out calibration'
        ) +
        theme_bw(base_size = 11)
    save_figure_bundle(
        plot, 'Fig-INLA-02_heldout_calibration', rows,
        'Held-out observed versus predicted relative mortality by programme and bleaching event.',
        'Assesses calibration of the selected operational model under event transfer.',
        'Predictions are cross-validated; residual structure may reflect unresolved observation and disturbance processes.',
        model_id, 'INLA', 'heldout_calibration', 'operational_validation',
        root, TRUE
    )
    list(plot = plot, data = rows)
}

plot_empirical_dhw_response <- function(model_id, root = project_root()) {
    prediction_path <- if (model_id %in% c(
        'prospective_cots_log_peak_recency',
        'operational_rrn_raw_plus_manta_state'
    )) file.path(root, 'output', 'prospective_cots_nowcast',
                 'cv_predictions.csv') else if (identical(
        model_id, 'cots_raw_interval_logistic20_5_cyclone'
    )) file.path(root, 'output', 'cots_raw_enso_occurrence',
                 'cots_cv_predictions.csv') else if (identical(
        model_id, 'baseline_cots_logistic20_5_cyclone'
    )) file.path(root, 'output', 'cots_timing_cover_cyclone_soft_gate',
                 'cv_predictions.csv') else if (identical(
        model_id, 'cause_aware_exposure_gated_competing_hazards'
    )) file.path(root, 'output', 'cause_aware_competing_hazards',
                 'gated_cv_predictions.csv') else
        file.path(root, 'output', 'inla_spatiotemporal', 'cv_predictions.csv')
    prediction_id <- if (identical(
        model_id, 'cots_raw_interval_logistic20_5_cyclone'
    )) 'cots_raw_interval_relative' else model_id
    rows <- read_csv(prediction_path, show_col_types = FALSE) |>
        filter(candidate == .env$prediction_id, scheme == 'leave_one_event_out') |>
        mutate(dhw_bin = cut(
            ann_maxdhw, breaks = seq(0, 22, by = 1), include.lowest = TRUE
        )) |>
        group_by(programme_key, event_year, dhw_bin) |>
        summarise(
            dhw = mean(ann_maxdhw),
            observed = mean(observed_mortality),
            predicted = mean(predicted_mortality),
            observations = n(),
            .groups = 'drop'
        ) |>
        filter(is.finite(dhw))
    long <- rows |>
        pivot_longer(c(observed, predicted), names_to = 'series',
                     values_to = 'mortality')
    plot <- ggplot(long, aes(dhw, mortality, colour = series)) +
        geom_line(linewidth = 0.8) +
        geom_point(aes(size = observations), alpha = 0.7) +
        facet_grid(programme_key ~ event_year) +
        scale_colour_manual(values = c(observed = '#222222', predicted = '#0072B2')) +
        labs(
            x = 'Local-first annual maximum DHW',
            y = 'Mean mortality within DHW bin',
            colour = NULL, size = 'Observations',
            title = 'Empirical held-out DHW response by event'
        ) +
        theme_bw(base_size = 10)
    save_figure_bundle(
        plot, 'Fig-INLA-03_empirical_event_dhw_response', rows,
        'Binned held-out observed and predicted mortality across DHW by programme and event.',
        'Visualises event-specific response patterns that aggregate performance metrics can hide.',
        'This is an empirical diagnostic, not a covariate-adjusted partial dependence curve.',
        model_id, 'INLA', 'empirical_dhw_response', 'operational_validation',
        root, TRUE
    )
    list(plot = plot, data = rows)
}

plot_predictor_correlation <- function(model_id, root = project_root()) {
    fixed_model_id <- if (model_id %in% c(
        'cause_aware_exposure_gated_competing_hazards',
        'baseline_cots_logistic20_5_cyclone',
        'cots_raw_interval_logistic20_5_cyclone',
        'prospective_cots_log_peak_recency',
        'operational_rrn_raw_plus_manta_state'
    )) 'persistent_rw1_local_dhw_decomposed_hazards_freshwater_partial_pool' else
        model_id
    model_rows_path <- file.path(
        root, 'output', 'explanatory_event_dhw', 'event_dhw_brt_data.csv'
    )
    fixed_path <- file.path(
        root, 'output', 'inla_spatiotemporal', 'full_fixed_effects.csv'
    )
    if (!file.exists(model_rows_path) || !file.exists(fixed_path)) {
        stop('Correlation inputs are missing; run the explanatory and INLA pipelines')
    }
    model_rows <- read_csv(model_rows_path, show_col_types = FALSE)
    fixed_terms <- read_csv(fixed_path, show_col_types = FALSE) |>
        filter(candidate == .env$fixed_model_id) |>
        pull(term) |>
        unique()
    predictor_terms <- intersect(fixed_terms, names(model_rows))
    predictor_terms <- predictor_terms[vapply(
        model_rows[predictor_terms], is.numeric, logical(1)
    )]
    nonconstant <- vapply(
        model_rows[predictor_terms],
        function(x) sd(x, na.rm = TRUE) > 0,
        logical(1)
    )
    predictor_terms <- predictor_terms[nonconstant]
    correlation <- cor(
        model_rows[predictor_terms], use = 'pairwise.complete.obs',
        method = 'pearson'
    )
    labels <- read_model_terms(root) |>
        select(id, label)
    label_lookup <- setNames(labels$label, labels$id)
    display_labels <- ifelse(
        predictor_terms %in% names(label_lookup),
        label_lookup[predictor_terms], predictor_terms
    )
    names(display_labels) <- predictor_terms
    long <- as.data.frame(as.table(correlation), stringsAsFactors = FALSE) |>
        as_tibble() |>
        rename(term_x = Var1, term_y = Var2, correlation = Freq) |>
        mutate(
            label_x = unname(display_labels[term_x]),
            label_y = unname(display_labels[term_y]),
            label_x = factor(label_x, levels = display_labels),
            label_y = factor(label_y, levels = rev(display_labels)),
            label = if_else(
                abs(correlation) >= 0.65 | term_x == term_y,
                sprintf('%.2f', correlation), ''
            )
        )
    plot <- ggplot(long, aes(label_x, label_y, fill = correlation)) +
        geom_tile(colour = 'white', linewidth = 0.12) +
        geom_text(aes(label = label), size = 2.1) +
        scale_fill_gradient2(
            low = '#2166AC', mid = 'white', high = '#B2182B',
            midpoint = 0, limits = c(-1, 1)
        ) +
        coord_equal() +
        labs(
            title = 'Correlation among selected-model fixed-effect predictors',
            subtitle = 'Pearson correlations on full-data, fold-standardised inputs; labels shown for |r| >= 0.65',
            x = NULL, y = NULL, fill = 'Correlation'
        ) +
        theme_minimal(base_size = 9) +
        theme(
            axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
            panel.grid = element_blank(), legend.position = 'bottom'
        )
    save_figure_bundle(
        plot, 'Fig-DATA-01_predictor_correlation', long,
        'Pairwise Pearson correlation matrix for every numeric fixed-effect predictor in the selected model.',
        'Highlights structural redundancy among DHW spline bases, derived interactions and environmental predictors before interpreting individual coefficients.',
        'High correlation does not by itself justify removing an ecological mechanism; derived hinge and interaction terms are expected to be correlated and candidate removal must be judged by held-out prediction.',
        model_id, 'data', 'predictor_correlation', 'model_audit',
        root, TRUE, width = 14, height = 12
    )
    pairs <- long |>
        mutate(x_index = match(term_x, predictor_terms),
               y_index = match(term_y, predictor_terms)) |>
        filter(x_index < y_index) |>
        arrange(desc(abs(correlation))) |>
        select(term_x, term_y, label_x, label_y, correlation)
    list(plot = plot, data = long, pairs = pairs)
}

run_standard_inla_diagnostics <- function(model_id, root = project_root()) {
    outputs <- list(
        fixed_effects = plot_fixed_effect_posteriors(model_id, root),
        calibration = plot_cv_calibration(model_id, root),
        dhw_response = plot_empirical_dhw_response(model_id, root),
        predictor_correlation = plot_predictor_correlation(model_id, root)
    )
    write_figure_readme(root)
    outputs
}
