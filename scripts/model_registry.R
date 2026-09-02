# Shared model registry and report-contract utilities.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(yaml)
    library(jsonlite)
    library(tibble)
})

project_root <- function(start = getwd()) {
    candidate <- normalizePath(start, winslash = '/', mustWork = TRUE)
    while (TRUE) {
        if (file.exists(file.path(candidate, 'config', 'model_registry.yml'))) {
            return(candidate)
        }
        parent <- dirname(candidate)
        if (identical(parent, candidate)) {
            stop('Could not locate project root from: ', start)
        }
        candidate <- parent
    }
}

read_model_registry <- function(root = project_root()) {
    yaml::read_yaml(file.path(root, 'config', 'model_registry.yml'))
}

read_model_terms <- function(root = project_root()) {
    bind_rows(lapply(
        yaml::read_yaml(file.path(root, 'config', 'model_terms.yml'))$terms,
        function(x) tibble::as_tibble_row(x) |>
            mutate(operational = as.character(operational))
    ))
}

candidate_table <- function(registry = read_model_registry()) {
    bind_rows(lapply(registry$candidates, tibble::as_tibble_row))
}

resolve_fit_artifact <- function(model_id, root = project_root()) {
    registry <- read_model_registry(root)
    if (identical(model_id, registry$best_model$id) &&
        !is.null(registry$best_model$fit_artifact)) {
        artifact <- file.path(root, registry$best_model$fit_artifact)
        if (!file.exists(artifact)) stop('Registered fit artefact is missing: ', artifact)
        return(artifact)
    }
    candidates <- list.files(
        file.path(root, 'output', 'inla_spatiotemporal'),
        pattern = paste0('^full_.*', model_id, '\\.rds$'),
        full.names = TRUE
    )
    if (length(candidates) == 0L) {
        stop('No full INLA fit artefact found for model: ', model_id)
    }
    candidates[order(file.info(candidates)$mtime, decreasing = TRUE)][[1]]
}

model_metrics <- function(root = project_root()) {
    metrics <- read_csv(
        file.path(root, 'output', 'inla_spatiotemporal', 'model_comparison.csv'),
        show_col_types = FALSE
    )
    cause_path <- file.path(
        root, 'output', 'cause_aware_competing_hazards',
        'model_comparison.csv'
    )
    if (file.exists(cause_path)) {
        metrics <- bind_rows(
            metrics,
            read_csv(cause_path, show_col_types = FALSE) |>
                filter(candidate != 'selected_current')
        )
    }
    timing_path <- file.path(
        root, 'output', 'cots_timing_cover_cyclone_soft_gate',
        'model_comparison.csv'
    )
    if (file.exists(timing_path)) {
        metrics <- bind_rows(
            metrics,
            read_csv(timing_path, show_col_types = FALSE) |>
                filter(candidate != 'cause_aware_exposure_gated_competing_hazards')
        )
    }
    raw_cots_path <- file.path(
        root, 'output', 'cots_raw_enso_occurrence',
        'cots_model_comparison.csv'
    )
    if (file.exists(raw_cots_path)) {
        metrics <- bind_rows(
            metrics,
            read_csv(raw_cots_path, show_col_types = FALSE) |>
                filter(candidate != 'current_log_interval_selected') |>
                mutate(candidate = recode(
                    candidate,
                    cots_raw_interval_relative =
                        'cots_raw_interval_logistic20_5_cyclone'
                ))
        )
    }
    prospective_path <- file.path(
        root, 'output', 'prospective_cots_nowcast', 'model_comparison.csv'
    )
    if (file.exists(prospective_path)) {
        metrics <- bind_rows(
            metrics,
            read_csv(prospective_path, show_col_types = FALSE) |>
                filter(candidate != 'selected_retrospective_raw_interval')
        )
    }
    ensemble_path <- file.path(
        root, 'output', 'inla_brt_residual_ensemble', 'model_comparison.csv'
    )
    if (file.exists(ensemble_path)) {
        metrics <- bind_rows(
            metrics,
            read_csv(ensemble_path, show_col_types = FALSE) |>
                filter(grepl('residual BRT', candidate)) |>
                mutate(candidate = recode(
                    candidate,
                    `INLA + balanced residual BRT` =
                        'inla_brt_balanced_residual',
                    `INLA + severe-weighted residual BRT` =
                        'inla_brt_severe_weighted_residual'
                ))
        )
    }
    full_brt_path <- file.path(
        root, 'output', 'standalone_brt_envelope_ensemble',
        'model_comparison.csv'
    )
    if (file.exists(full_brt_path)) {
        metrics <- bind_rows(
            metrics,
            read_csv(full_brt_path, show_col_types = FALSE) |>
                filter(candidate %in% c(
                    'Standalone direct BRT',
                    'Standalone two-part BRT',
                    'Direct smooth 50% ensemble',
                    'Direct applicability champion'
                )) |>
                mutate(candidate = recode(
                    candidate,
                    `Standalone direct BRT` =
                        'standalone_direct_brt_full_predictors',
                    `Standalone two-part BRT` =
                        'standalone_two_part_brt_full_predictors',
                    `Direct smooth 50% ensemble` =
                        'direct_brt_smooth_environmental_envelope',
                    `Direct applicability champion` =
                        'direct_brt_applicability_environmental_envelope'
                ))
        )
    }
    registry <- read_model_registry(root)
    candidate_table(registry) |>
        left_join(metrics, by = c('id' = 'candidate')) |>
        arrange(factor(
            status,
            levels = c(
                'selected', 'comparator', 'under_test',
                'ensemble_or_uncertainty_only', 'explanatory_only', 'rejected'
            )
        ))
}

selected_model_snapshot <- function(root = project_root()) {
    registry <- read_model_registry(root)
    model_id <- registry$best_model$id
    fit_path <- resolve_fit_artifact(model_id, root)
    comparison_path <- if (!is.null(registry$best_model$comparison_artifact)) {
        file.path(root, registry$best_model$comparison_artifact)
    } else file.path(root, 'output', 'inla_spatiotemporal', 'model_comparison.csv')
    list(
        registry_version = registry$registry_version,
        analysis = registry$analysis,
        best_model = registry$best_model,
        resolved_fit_artifact = normalizePath(fit_path, winslash = '/'),
        fit_md5 = unname(tools::md5sum(fit_path)),
        comparison_md5 = unname(tools::md5sum(comparison_path)),
        generated_utc = format(
            as.POSIXct(Sys.time(), tz = 'UTC'), '%Y-%m-%dT%H:%M:%SZ'
        )
    )
}

write_model_snapshot <- function(root = project_root()) {
    output_dir <- file.path(root, 'output', 'model_registry')
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    snapshot <- selected_model_snapshot(root)
    write_json(
        snapshot, file.path(output_dir, 'current_model.json'),
        auto_unbox = TRUE, pretty = TRUE
    )
    write_csv(
        read_model_terms(root),
        file.path(output_dir, 'current_model_terms.csv')
    )
    write_csv(
        model_metrics(root),
        file.path(output_dir, 'candidate_status.csv')
    )
    invisible(snapshot)
}

current_model_latex <- function(model_id) {
    cause_model_ids <- c(
        'cause_aware_exposure_gated_competing_hazards',
        'baseline_cots_logistic20_5_cyclone',
        'cots_raw_interval_logistic20_5_cyclone',
        'prospective_cots_log_peak_recency',
        'operational_rrn_raw_plus_manta_state'
    )
    if (model_id %in% cause_model_ids) {
        cyclone_activation <- if (model_id %in% c(
            'baseline_cots_logistic20_5_cyclone',
            'cots_raw_interval_logistic20_5_cyclone',
            'prospective_cots_log_peak_recency',
            'operational_rrn_raw_plus_manta_state'
        )) '\\operatorname{logit}^{-1}[(W_i-20)/5]' else
            '\\mathbb{1}(W_i>20)'
        cots_intensity <- if (identical(
            model_id, 'operational_rrn_raw_plus_manta_state'
        )) paste0(
            '\\gamma^{(k)}_2(I^{\\rm RRN}_i-0.22)_+',
            '+\\gamma^{(k)}_3(I^{\\rm Manta}_i-0.22)_+',
            '+\\gamma^{(k)}_4(P^{\\rm Manta}_i-0.22)_+',
            '+\\gamma^{(k)}_5\\tau^{\\rm Manta}_i'
        ) else if (identical(
            model_id, 'prospective_cots_log_peak_recency'
        )) paste0(
            '\\gamma^{(k)}_2\\log[1+(L_i-0.22)_+]',
            '+\\gamma^{(k)}_3\\log[1+(P_i-0.22)_+]',
            '+\\gamma^{(k)}_4T_i'
        ) else if (identical(
            model_id, 'cots_raw_interval_logistic20_5_cyclone'
        )) '\\gamma^{(k)}_2(I_i-0.22)_+' else
            '\\gamma^{(k)}_2p_i\\log[1+(I_i-0.22)_+]'
        return(paste(
            '\\begin{aligned}',
            '\\eta^{(k)}_{T,i} ={}& \\alpha^{(k)}_{L_i}',
            '+ \\beta^{(k)}_1D_{zi}+\\beta^{(k)}_2H_{4,zi}',
            '+\\beta^{(k)}_3H_{8,zi} \\\\',
            '&+\\boldsymbol\\beta^{(k)\\top}_{\\rm composition}(A_i,C_i)',
            '+\\boldsymbol\\beta^{(k)\\top}_{\\rm history}(Q_i,N_i,E_i,R_i) \\\\',
            '&+\\boldsymbol\\beta^{(k)\\top}_{\\rm environment}',
            '(S_i,\\mathrm{Cloud}_i,J_i,\\mathrm{SSTskew}_i,',
            '\\mathrm{SSTkurt}_i,\\mathrm{Chl}_i) \\\\',
            '&+\\boldsymbol\\beta^{(k)\\top}_{\\rm freshwater}(G_i,W_i,P_i,U_i)',
            '+\\boldsymbol\\beta^{(k)\\top}_{\\rm DHW\\ interactions}\\mathbf h_i \\\\',
            '&+u^{(k)}_{\\mathrm{reef-event}(i)}+v^{(k)}_{\\mathrm{event}(i)}',
            '+s^{(k)}(\\mathbf x_i)+q^{(k)}_{\\mathrm{event}(i)}D_{zi}',
            '+\\delta^{(k)}_{L_i}, \\\\',
            'T_i &= \\operatorname{logit}^{-1}(\\eta^{(0)}_{T,i})',
            '\\operatorname{logit}^{-1}(\\eta^{(+)}_{T,i}), \\\\',
            '\\eta^{(k)}_{C,i} ={}& \\alpha^{(k)}_C+',
            paste0('\\gamma^{(k)}_1p_i+', cots_intensity),
            '+\\boldsymbol\\gamma^{(k)\\top}_3\\mathbf z_i, \\\\',
            if (identical(model_id, 'operational_rrn_raw_plus_manta_state'))
                'C_i &= a_i\\pi_{C,i}\\mu_{C,i},\\quad a_i=\\max\\{p_i,1-\\exp[-(I^{\\rm RRN}_i-0.22)_+/0.5]\\},' else if (identical(model_id, 'prospective_cots_log_peak_recency'))
                'C_i &= a_i\\pi_{C,i}\\mu_{C,i},\\quad a_i=\\max\\{p_i,1-\\exp[-(L_i-0.22)_+/0.5]\\},' else
                'C_i &= \\mathbb{1}(I_i>0.22)\\pi_{C,i}\\mu_{C,i},',
            '\\quad (\\pi_{C,i},\\mu_{C,i})=',
            '\\operatorname{logit}^{-1}(\\eta^{(0)}_{C,i},\\eta^{(+)}_{C,i}), \\\\',
            '\\eta^{(k)}_{S,i} ={}& \\alpha^{(k)}_S+',
            '\\theta^{(k)}_1\\log(1+W_i)+\\theta^{(k)}_2V_i',
            '+\\theta^{(k)}_3R_i+\\boldsymbol\\theta^{(k)\\top}_4\\mathbf z_i, \\\\',
            paste0('S_i &= ', cyclone_activation,
                   '\\pi_{S,i}\\mu_{S,i},'),
            '\\quad (\\pi_{S,i},\\mu_{S,i})=',
            '\\operatorname{logit}^{-1}(\\eta^{(0)}_{S,i},\\eta^{(+)}_{S,i}), \\\\',
            '\\widehat{M}_i &= 1-(1-T_i)(1-C_i)(1-S_i).',
            '\\end{aligned}', sep = '\n'
        ))
    }
    if (!identical(
        model_id,
        'persistent_rw1_local_dhw_decomposed_hazards_freshwater_partial_pool'
    )) {
        stop('No LaTex template defined for selected model: ', model_id)
    }
    paste(
        '\\begin{aligned}',
        '\\eta_i ={}& \\alpha_{L_i} + \\beta_1D_{zi} + \\beta_2H_{4,zi}',
        '+ \\beta_3H_{8,zi} \\\\',
        '&+ \\boldsymbol{\\beta}_{\\rm composition}^{\\top}(A_i,C_i)',
        '+ \\boldsymbol{\\beta}_{\\rm history}^{\\top}(Q_i,N_i,E_i,R_i) \\\\',
        '&+ \\boldsymbol{\\beta}_{\\rm environment}^{\\top}',
        '(S_i,\\mathrm{Cloud}_i,J_i,\\mathrm{SSTskew}_i,',
        '\\mathrm{SSTkurt}_i,\\mathrm{Chl}_i) \\\\',
        '&+ \\boldsymbol{\\beta}_{\\rm freshwater}^{\\top}(G_i,W_i,P_i,U_i)',
        '+ \\boldsymbol{\\beta}_{\\rm cyclone}^{\\top}(X_i,T_i,V_i)',
        '+ \\boldsymbol{\\beta}_{\\rm COTS}^{\\top}(B_i,O_i) \\\\',
        '&+ \\beta_{DA}H_{4,zi}A_i + \\beta_{DN}H_{4,zi}N_i',
        '+ \\beta_{DC}H_{4,zi}\\mathrm{Cloud}_i + \\beta_{DJ}D_{zi}J_i \\\\',
        '&+ \\beta_{DE}D_{zi}E_i + \\beta_{DR}D_{zi}R_i',
        '+ \\beta_{DG}H_{4,zi}G_i + \\beta_{DW}H_{4,zi}W_i',
        '+ \\beta_{DP}H_{4,zi}P_i + \\beta_{DU}H_{4,zi}U_i \\\\',
        '&+ u_{\\mathrm{reef-event}(i)} + v_{\\mathrm{event}(i)}',
        '+ s(\\mathbf{x}_i) + q_{\\mathrm{event}(i)}D_{zi}',
        '+ \\delta_{L_i}.',
        '\\end{aligned}',
        sep = '\n'
    )
}

render_model_contract <- function(root = project_root()) {
    registry <- read_model_registry(root)
    best <- registry$best_model
    metrics <- model_metrics(root) |>
        filter(id == best$id, scheme == best$validation_primary)
    if (nrow(metrics) != 1L) {
        stop('Expected exactly one primary metric row for selected model')
    }
    write_model_snapshot(root)
    metric_text <- paste0(
        'RMSE = ', sprintf('%.3f', metrics$rmse),
        '; predictive R2 = ', sprintf('%.3f', metrics$predictive_r2),
        '; severe RMSE = ', sprintf('%.3f', metrics$severe_rmse),
        '; false-extreme rate = ', scales::percent(
            metrics$false_extreme_rate, accuracy = 0.1
        ), '.'
    )
    candidates <- model_metrics(root) |>
        transmute(
            ID = id, Role = role, Status = status,
            Validation = coalesce(scheme, 'not fitted in this comparison'),
            RMSE = if_else(is.na(rmse), NA_character_, sprintf('%.3f', rmse)),
            Predictive_R2 = if_else(
                is.na(predictive_r2), NA_character_,
                sprintf('%.3f', predictive_r2)
            )
        ) |>
        knitr::kable(format = 'pipe') |>
        paste(collapse = '\n')
    if (identical(best$framework, 'INLA composite')) {
        terms <- tibble(
            Mechanism = c(
                'Thermal/freshwater', 'COTS occurrence', 'COTS severity',
                'Cyclone occurrence/severity', 'Cause activation',
                'Observation process'
            ),
            Predictor = c(
                'Local-first DHW spline, composition, prior exposure, WQC, rainfall, cloud, current, SST shape and chlorophyll',
                if (identical(best$id, 'operational_rrn_raw_plus_manta_state')) 'Event-year RRN seasonal density, event-year hindcast probability, latest and prior-three-year Manta pressure and time since peak' else if (identical(best$id, 'prospective_cots_log_peak_recency')) 'Event-start hindcast probability, latest manta density, prior three-year peak density and time since peak' else if (identical(best$id, 'cots_raw_interval_logistic20_5_cyclone')) 'Hindcast outbreak probability and raw RRN density excess above 0.22 COTS/tow' else 'Hindcast outbreak probability and probability-weighted log RRN excess density',
                if (identical(best$id, 'operational_rrn_raw_plus_manta_state')) 'Raw excess above 0.22 COTS/tow for RRN, latest Manta and prior Manta peak, with cover, composition and geography' else if (identical(best$id, 'prospective_cots_log_peak_recency')) 'Log excess above 0.22 COTS/tow for latest and prior peak density, with cover, composition and geography' else if (identical(best$id, 'cots_raw_interval_logistic20_5_cyclone')) 'Raw RRN excess density, cover, composition and geography' else 'Probability-weighted log RRN excess density, cover, composition and geography',
                'Damaging-wave hours, cyclone wind-distance index, rainfall, cover, composition and geography',
                if (best$id %in% c('baseline_cots_logistic20_5_cyclone', 'cots_raw_interval_logistic20_5_cyclone', 'prospective_cots_log_peak_recency', 'operational_rrn_raw_plus_manta_state')) {
                    'RRN COTS density >0.22 COTS/tow; continuous cyclone activation centred at 20 damaging-wave hours'
                } else {
                    'RRN COTS density >0.22 COTS/tow; damaging-wave exposure >20 hours'
                },
                'Programme-specific occurrence and positive-magnitude layers in the thermal response'
            ),
            Transformation = c(
                'Fold-standardised; DHW hinges at 4 and 8',
                if (identical(best$id, 'operational_rrn_raw_plus_manta_state')) 'Raw RRN event excess; dated Manta through end of February; raw latest/peak excess plus years since peak' else if (identical(best$id, 'prospective_cots_log_peak_recency')) 'Available through end of February; log[1 + max(density - 0.22, 0)] plus years since peak' else 'p x log[1 + max(I - 0.22, 0)]',
                'Two-part Bernoulli-beta annual-transition model',
                'log(1 + wave hours), scaled wind-distance and log rainfall',
                if (identical(best$id, 'operational_rrn_raw_plus_manta_state')) {
                    'Soft activation from hindcast probability or RRN seasonal density; logistic cyclone weight logit^-1[(hours - 20)/5]'
                } else if (identical(best$id, 'prospective_cots_log_peak_recency')) {
                    'Soft event-start activation from hindcast probability or latest density; logistic cyclone weight logit^-1[(hours - 20)/5]'
                } else if (identical(best$id, 'baseline_cots_logistic20_5_cyclone')) {
                    'Hard COTS gate; logistic cyclone weight logit^-1[(hours - 20)/5]'
                } else 'Exposure gate applied after cause-hazard prediction',
                'Bernoulli occurrence times conditional beta magnitude'
            ),
            Operational = c('TRUE', 'TRUE', 'TRUE', 'TRUE', 'TRUE', 'TRUE')
        ) |>
            knitr::kable(format = 'pipe') |>
            paste(collapse = '\n')
    } else {
        included_terms <- readr::read_csv(
            file.path(root, 'output', 'inla_spatiotemporal', 'full_fixed_effects.csv'),
            show_col_types = FALSE
        ) |>
            filter(candidate == best$id) |>
            pull(term) |>
            unique()
        terms <- read_model_terms(root) |>
            filter(id %in% included_terms) |>
            transmute(
                Mechanism = class, Predictor = label,
                Transformation = transformation, Operational = operational
            ) |>
            knitr::kable(format = 'pipe') |>
            paste(collapse = '\n')
    }
    caveats <- paste0('- ', unlist(best$caveats), collapse = '\n')

    knitr::asis_output(paste(
        '## Current model contract',
        '',
        paste0('Selected model: ', best$id, ' (', best$framework, '; ',
               best$purpose, ').'),
        '',
        paste0('Primary validation: ', best$validation_primary, '; ',
               metric_text),
        '',
        '### Response and observation model',
        '',
        '$$Z_i=\\mathbb{1}(Y_i>0),\\quad Z_i\\sim\\mathrm{Bernoulli}(\\pi_i),$$',
        '$$Y_i\\mid Z_i=1\\sim\\mathrm{Beta}(\\mu_i\\phi_{p(i)},',
        '(1-\\mu_i)\\phi_{p(i)}),$$',
        '$$\\mathrm{logit}(\\pi_i)=\\eta_i^{(0)},\\quad',
        '\\mathrm{logit}(\\mu_i)=\\eta_i^{(+)}.$$',
        '',
        'Six likelihood layers represent occurrence and positive conditional',
        'magnitude for LTMP, manta and MMP in the thermal/freshwater component.',
        'COTS- and cyclone-labelled annual cover transitions train two separate',
        'Bernoulli-beta hazards; held-out event intervals are excluded to avoid leakage.',
        '',
        '### Ecological linear predictor',
        '',
        '$$', current_model_latex(best$id), '$$',
        '',
        'DHW is local-first calibrated before its hinges and novelty are',
        'recomputed. The three expected hazards combine on the mortality scale,',
        'so explicitly labelled COTS and storm rows do not flatten the thermal',
        'dose-response. Disturbance labels train the cause components but are',
        'never required as GBR-wide prediction inputs.',
        '',
        '### Candidate ledger', '', candidates,
        '',
        '### Predictor dictionary', '', terms,
        '',
        '### Current caveats', '', caveats,
        sep = '\n'
    ))
}
