# Joint Bernoulli--Beta INLA screen with shared ecological effects and a small
# set of spatial and spatio-temporal alternatives. The same ecological
# coefficients and latent fields feed six observation likelihoods: occurrence
# and conditional magnitude for LTMP, manta tow and MMP.

suppressPackageStartupMessages({
    library(dplyr)
    library(INLA)
    library(Matrix)
    library(readr)
    library(sf)
    library(tidyr)
})
source('scripts/joint_compound_model_helpers.R')

set.seed(20260827L)
output_dir <- 'output/inla_spatiotemporal'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

programmes <- c('ltmp', 'manta', 'mmp')
event_years <- c(2016L, 2017L, 2020L, 2022L, 2024L)
families <- c(rep('binomial', 3L), rep('beta', 3L))
layer_lookup <- tibble(
    ID = seq_along(families),
    component = rep(c('occurrence', 'magnitude'), each = length(programmes)),
    programme_key = rep(programmes, times = 2L)
)

candidates <- tribble(
    ~candidate, ~spatial_structure, ~temporal_dhw_slope, ~feature_set,
    ~partial_pool, ~lizard_dhw_adjustment,
    'nonspatial', 'none', FALSE, 'core', FALSE, 'none',
    'persistent_spatial', 'persistent', FALSE, 'core', FALSE, 'none',
    'event_spatial_iid', 'event_iid', FALSE, 'core', FALSE, 'none',
    'event_spatial_ar1', 'event_ar1', FALSE, 'core', FALSE, 'none',
    'persistent_plus_dhw_rw1', 'persistent', TRUE, 'core', FALSE, 'none',
    'persistent_rw1_current', 'persistent', TRUE, 'current', FALSE, 'none',
    'persistent_rw1_current_sst_chla', 'persistent', TRUE,
    'current_sst_chla', FALSE, 'none',
    'persistent_rw1_full_partial_pool', 'persistent', TRUE,
    'current_sst_chla', TRUE, 'none',
    'persistent_rw1_cyclone_cots_partial_pool', 'persistent', TRUE,
    'current_sst_chla_pressure', TRUE, 'none',
    'persistent_rw1_cyclone_cots_disease_partial_pool', 'persistent', TRUE,
    'current_sst_chla_pressure_disease', TRUE, 'none',
    'persistent_rw1_cyclone_cots_lizard_adjusted_partial_pool',
    'persistent', TRUE, 'current_sst_chla_pressure', TRUE, 'flat_mean',
    'persistent_rw1_cyclone_cots_lizard_median_adjusted_partial_pool',
    'persistent', TRUE, 'current_sst_chla_pressure', TRUE, 'all_median',
    'persistent_rw1_cyclone_cots_lizard_slope_adjusted_partial_pool',
    'persistent', TRUE, 'current_sst_chla_pressure', TRUE, 'slope_mean',
    'persistent_rw1_cyclone_cots_rolling_dhw_partial_pool',
    'persistent', TRUE, 'current_sst_chla_pressure', TRUE, 'rolling_prior',
    'persistent_rw1_cyclone_cots_event_dhw_partial_pool',
    'persistent', TRUE, 'current_sst_chla_pressure', TRUE, 'within_event'
)

lizard_cluster_reef_ids <- c(
    '14-114', '14-116a', '14-118', '14-123', '14-126', '14-143'
)
lizard_logger_dhw <- read_csv(
    'data/processed/lizard_cluster_logger_dhw_uplift_2024.csv',
    show_col_types = FALSE
)
lizard_logger_uplifts <- c(
    none = 0,
    flat_mean = lizard_logger_dhw |>
        filter(grepl('FL', series)) |>
        summarise(x = mean(dhw_uplift_vs_official_noaa, na.rm = TRUE)) |>
        pull(x),
    all_median = median(
        lizard_logger_dhw$dhw_uplift_vs_official_noaa, na.rm = TRUE
    ),
    slope_mean = lizard_logger_dhw |>
        filter(grepl('SL', series)) |>
        summarise(x = mean(dhw_uplift_vs_official_noaa, na.rm = TRUE)) |>
        pull(x)
)
stopifnot(all(is.finite(lizard_logger_uplifts)))

raw_predictors <- c(
    'ann_maxdhw', 'dhw_excess4', 'dhw_excess8',
    'prop_acropora_pre', 'observed_pre_cover',
    'dhw10_load4', 'dhw_novelty10',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8',
    'secc3m_p10', 'cloudp_90',
    'log_coastal_rain30', 'wqc_prior10_percentile',
    'log1p_cyc_maxHrs4mw', 'log1p_cot_idwmeanpertow',
    'mcur_90', 'sst_summer_skewness',
    'sst_summer_excess_kurtosis', 'chla_wetseason_median',
    'tc_proximity100', 'tc_wind_distance_index',
    'cots_outbreak_probability',
    'disease_risk_max', 'disease_risk_days_ge1'
)

core_shared_terms <- c(
    paste0(raw_predictors[!raw_predictors %in% c(
        'mcur_90', 'sst_summer_skewness',
        'sst_summer_excess_kurtosis', 'chla_wetseason_median',
        'tc_proximity100', 'tc_wind_distance_index',
        'cots_outbreak_probability',
        'disease_risk_max', 'disease_risk_days_ge1'
    )], '_z'),
    'dhw_x_acropora', 'dhw_x_novelty', 'dhw_x_freshwater',
    'dhw_x_cloud', 'dhw_x_repeat_exposure', 'dhw_x_recovery_interval'
)
current_shared_terms <- c('mcur_90_z', 'dhw_x_current')
sst_chla_shared_terms <- c(
    'sst_summer_skewness_z', 'sst_summer_excess_kurtosis_z',
    'chla_wetseason_median_z'
)
pressure_shared_terms <- c(
    'tc_proximity100_z', 'tc_wind_distance_index_z',
    'cots_outbreak_probability_z'
)
disease_shared_terms <- c(
    'disease_risk_max_z', 'disease_risk_days_ge1_z'
)
shared_terms <- c(
    core_shared_terms, current_shared_terms, sst_chla_shared_terms,
    pressure_shared_terms, disease_shared_terms
)

candidate_shared_terms <- function(candidate) {
    switch(
        as.character(candidate$feature_set),
        core = core_shared_terms,
        current = c(core_shared_terms, current_shared_terms),
        current_sst_chla = c(
            core_shared_terms, current_shared_terms, sst_chla_shared_terms
        ),
        current_sst_chla_pressure = c(
            core_shared_terms, current_shared_terms, sst_chla_shared_terms,
            pressure_shared_terms
        ),
        current_sst_chla_pressure_disease = c(
            core_shared_terms, current_shared_terms, sst_chla_shared_terms,
            pressure_shared_terms, disease_shared_terms
        ),
        stop('Unknown feature set: ', candidate$feature_set)
    )
}

new_features <- read_csv(
    'data/processed/sst_chla_features_validation.csv',
    show_col_types = FALSE
) |>
    select(
        ReefID, year, sst_summer_skewness,
        sst_summer_excess_kurtosis, chla_wetseason_median
    )

cyclone <- read_csv(
    'data/processed/bom_cyclone_reef_year.csv', show_col_types = FALSE
) |>
    select(
        ReefID, event_year, tc_min_distance_km,
        tc_nearest_name, tc_nearest_max_wind_ms
    )

cots_hindcast <- read_csv(
    'data/gbrPredsAdj_20262408.csv', show_col_types = FALSE
) |>
    transmute(
        ReefName = reefName, event_year = as.integer(year),
        cots_outbreak_probability = outbrProb
    )

disease <- read_csv(
    'data/processed/noaa_disease_risk_validation.csv',
    show_col_types = FALSE
) |>
    select(
        ReefID, year, disease_risk_applicable,
        disease_risk_max, disease_risk_days_ge1,
        disease_risk_burden
    )

dhw_correction_path <- paste0(
    'data/processed/', 'noaa_dhw_correction_layer_validation.csv'
)
dhw_correction_hash <- substr(
    unname(tools::md5sum(dhw_correction_path)), 1, 8
)
candidate_cache_version <- function(candidate_name) {
    if (grepl('_(rolling|event)_dhw_', candidate_name)) {
        paste0('v3_', dhw_correction_hash)
    } else {
        'v2'
    }
}

dhw_correction <- read_csv(
    dhw_correction_path,
    show_col_types = FALSE
) |>
    select(
        ReefID, event_year, rolling_prior_correction,
        within_event_correction, nearest_current_logger_km,
        effective_current_loggers
    )

data <- load_joint_compound_rows() |>
    left_join(
        new_features,
        by = c('ReefID', 'event_year' = 'year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        cyclone, by = c('ReefID', 'event_year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        cots_hindcast, by = c('ReefName', 'event_year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        disease, by = c('ReefID', 'event_year' = 'year'),
        relationship = 'many-to-one'
    ) |>
    left_join(
        dhw_correction, by = c('ReefID', 'event_year'),
        relationship = 'many-to-one'
    ) |>
    mutate(
        tc_proximity100 = exp(-pmin(tc_min_distance_km, 1000) / 100),
        tc_wind_distance_index = pmax(
            tc_nearest_max_wind_ms - 17, 0
        ) * tc_proximity100,
        event_index = match(as.integer(event_year), event_years),
        reef_event_key = paste(ReefID, event_year, sep = '__')
    )

reef_locations <- data |>
    distinct(ReefID, lon, lat) |>
    st_as_sf(coords = c('lon', 'lat'), crs = 4326) |>
    st_transform(3577)
reef_xy <- st_coordinates(reef_locations) / 1000

hull <- inla.nonconvex.hull(reef_xy, convex = 80, concave = 160)
mesh <- inla.mesh.2d(
    boundary = hull, max.edge = c(80, 240), cutoff = 20,
    offset = c(80, 300)
)
spde <- inla.spde2.pcmatern(
    mesh = mesh, constr = TRUE,
    prior.range = c(200, 0.5),
    prior.sigma = c(0.7, 0.05)
)

location_lookup <- tibble(
    ReefID = reef_locations$ReefID,
    x_km = reef_xy[, 1], y_km = reef_xy[, 2]
)
data <- left_join(data, location_lookup, by = 'ReefID')

mesh_summary <- tibble(
    nodes = mesh$n,
    triangles = nrow(mesh$graph$tv),
    reefs = n_distinct(data$ReefID),
    reef_events = n_distinct(data$reef_event_key),
    events = n_distinct(data$event_year),
    crs = 'EPSG:3577; coordinates in km',
    max_edge_inner_km = 80,
    max_edge_outer_km = 240,
    pc_range_km = 200,
    pc_range_probability = 0.5,
    pc_sigma = 0.7,
    pc_sigma_probability = 0.05
)
write_csv(mesh_summary, file.path(output_dir, 'mesh_summary.csv'))
saveRDS(mesh, file.path(output_dir, 'mesh.rds'))

add_derived_terms <- function(rows) {
    rows |>
        mutate(
            freshwater_joint_z = (
                log_coastal_rain30_z + wqc_prior10_percentile_z
            ) / sqrt(2),
            dhw_x_acropora = dhw_excess4_z * prop_acropora_pre_z,
            dhw_x_novelty = dhw_excess4_z * dhw_novelty10_z,
            dhw_x_freshwater = dhw_excess4_z * freshwater_joint_z,
            dhw_x_cloud = dhw_excess4_z * cloudp_90_z,
            dhw_x_repeat_exposure = ann_maxdhw_z *
                dhw_events_since2016_n6_z,
            dhw_x_recovery_interval = ann_maxdhw_z *
                dhw_years_since_last_n6_capped8_z,
            dhw_x_current = ann_maxdhw_z * mcur_90_z,
            event_index = match(as.integer(event_year), event_years),
            reef_event_key = paste(ReefID, event_year, sep = '__')
        )
}

prepare_fold <- function(analysis, assessment) {
    prepared <- prepare_fold_predictors(
        analysis, assessment, raw_predictors
    )
    list(
        analysis = add_derived_terms(prepared$analysis),
        assessment = add_derived_terms(prepared$assessment),
        preprocessing = prepared$preprocessing
    )
}

apply_candidate_dhw <- function(rows, candidate) {
    adjustment <- as.character(candidate$lizard_dhw_adjustment)
    fixed_uplift <- if (adjustment %in% names(lizard_logger_uplifts)) {
        unname(lizard_logger_uplifts[[adjustment]])
    } else 0
    rows |>
        mutate(
            ann_maxdhw_original = ann_maxdhw,
            applied_dhw_uplift = case_when(
                adjustment == 'rolling_prior' ~
                    coalesce(rolling_prior_correction, 0),
                adjustment == 'within_event' ~
                    coalesce(within_event_correction, 0),
                fixed_uplift > 0 & event_year == 2024L &
                    ReefID %in% lizard_cluster_reef_ids ~ fixed_uplift,
                TRUE ~ 0
            ),
            ann_maxdhw = pmax(ann_maxdhw + applied_dhw_uplift, 0),
            dhw_excess4 = pmax(ann_maxdhw - 4, 0),
            dhw_excess8 = pmax(ann_maxdhw - 8, 0),
            # Novelty is current-event DHW minus the preceding 10-y maximum.
            dhw_novelty10 = dhw_novelty10 + applied_dhw_uplift
        )
}

likelihood_column <- function(programme, component) {
    base <- match(programme, programmes)
    if (component == 'magnitude') base + length(programmes) else base
}

make_response_matrix <- function(values, columns) {
    response <- matrix(
        NA_real_, nrow = length(values), ncol = length(families)
    )
    if (length(values) > 0L) {
        response[cbind(seq_along(values), columns)] <- values
    }
    response
}

make_fixed_rows <- function(rows, component, observed_event_keys,
                            observed_event_indices) {
    columns <- likelihood_column(rows$programme_key, component)
    intercepts <- matrix(
        0, nrow = nrow(rows), ncol = length(families),
        dimnames = list(NULL, paste0('layer_', seq_along(families)))
    )
    intercepts[cbind(seq_len(nrow(rows)), columns)] <- 1

    bind_cols(
        as_tibble(intercepts),
        rows |>
            transmute(
                across(all_of(shared_terms)),
                # Prediction-only levels must not enter constrained iid effects.
                # NA omits the term and gives the zero-centred new-level mean.
                reef_event_index = match(reef_event_key, observed_event_keys),
                event_effect_index = if_else(
                    event_index %in% observed_event_indices,
                    event_index, NA_integer_
                ),
                event_dhw_index = event_index,
                event_dhw_weight = ann_maxdhw_z,
                pool_acropora_index = columns,
                pool_acropora_weight = dhw_x_acropora,
                pool_novelty_index = columns,
                pool_novelty_weight = dhw_x_novelty,
                pool_freshwater_index = columns,
                pool_freshwater_weight = dhw_x_freshwater,
                pool_cloud_index = columns,
                pool_cloud_weight = dhw_x_cloud,
                pool_current_index = columns,
                pool_current_weight = dhw_x_current,
                likelihood_column = columns,
                x_km, y_km
            )
    )
}

make_stack_part <- function(rows, component, observed, tag,
                            observed_event_keys, observed_event_indices,
                            candidate) {
    fixed <- make_fixed_rows(
        rows, component, observed_event_keys, observed_event_indices
    )
    response <- make_response_matrix(observed, fixed$likelihood_column)
    locations <- as.matrix(select(fixed, x_km, y_km))

    A <- list(1)
    effects <- list(fixed)
    if (candidate$spatial_structure == 'persistent') {
        A[[length(A) + 1L]] <- inla.spde.make.A(mesh, loc = locations)
        effects[[length(effects) + 1L]] <- list(spatial_field = seq_len(mesh$n))
    }
    if (candidate$spatial_structure %in% c('event_iid', 'event_ar1')) {
        A[[length(A) + 1L]] <- inla.spde.make.A(
            mesh, loc = locations, group = rows$event_index,
            n.group = length(event_years)
        )
        effects[[length(effects) + 1L]] <- inla.spde.make.index(
            'event_spatial_field', n.spde = mesh$n,
            n.group = length(event_years)
        )
    }

    inla.stack(
        data = list(response = response), A = A, effects = effects,
        tag = tag, remove.unused = FALSE
    )
}

pc_prec <- list(
    prec = list(prior = 'pc.prec', param = c(0.7, 0.05))
)
pc_pool_prec <- list(
    prec = list(prior = 'pc.prec', param = c(0.5, 0.05))
)

make_formula <- function(candidate) {
    rhs <- c(
        paste0('layer_', seq_along(families)),
        candidate_shared_terms(candidate),
        'f(reef_event_index, model=\'iid\', constr=TRUE, hyper=pc_prec)',
        'f(event_effect_index, model=\'iid\', constr=TRUE, hyper=pc_prec)'
    )
    if (candidate$spatial_structure == 'persistent') {
        rhs <- c(rhs, 'f(spatial_field, model=spde)')
    }
    if (candidate$spatial_structure %in% c('event_iid', 'event_ar1')) {
        temporal_model <- if_else(
            candidate$spatial_structure == 'event_ar1', 'ar1', 'iid'
        )
        rhs <- c(rhs, paste0(
            'f(event_spatial_field, model=spde, ',
            'group=event_spatial_field.group, control.group=list(model=\'',
            temporal_model, '\'))'
        ))
    }
    if (candidate$temporal_dhw_slope) {
        rhs <- c(rhs, paste0(
            'f(event_dhw_index, event_dhw_weight, model=\'rw1\', ',
            'constr=TRUE, scale.model=TRUE, hyper=pc_prec)'
        ))
    }
    if (isTRUE(candidate$partial_pool)) {
        rhs <- c(
            rhs,
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
            )
        )
    }
    as.formula(paste('response ~ -1 +', paste(rhs, collapse = ' + ')))
}

fit_candidate <- function(analysis, assessment, candidate,
                          compute_criteria = FALSE, seed = 1L) {
    analysis <- apply_candidate_dhw(analysis, candidate)
    assessment <- apply_candidate_dhw(assessment, candidate)
    prepared <- prepare_fold(analysis, assessment)
    analysis <- prepared$analysis
    assessment <- prepared$assessment
    observed_event_keys <- unique(analysis$reef_event_key)
    observed_event_indices <- unique(analysis$event_index)

    positive <- analysis |>
        filter(mortality_prop > 0)
    upper_epsilon <- 0.5 / nrow(positive)
    magnitude_response <- pmin(
        positive$mortality_prop, 1 - upper_epsilon
    )

    stack <- inla.stack(
        make_stack_part(
            analysis, 'occurrence',
            as.numeric(analysis$mortality_prop > 0), 'analysis_occurrence',
            observed_event_keys, observed_event_indices, candidate
        ),
        make_stack_part(
            positive, 'magnitude', magnitude_response,
            'analysis_magnitude', observed_event_keys,
            observed_event_indices, candidate
        ),
        make_stack_part(
            assessment, 'occurrence', rep(NA_real_, nrow(assessment)),
            'prediction_occurrence', observed_event_keys,
            observed_event_indices, candidate
        ),
        make_stack_part(
            assessment, 'magnitude', rep(NA_real_, nrow(assessment)),
            'prediction_magnitude', observed_event_keys,
            observed_event_indices, candidate
        )
    )
    stack_data <- inla.stack.data(stack)
    link <- c(
        likelihood_column(analysis$programme_key, 'occurrence'),
        likelihood_column(positive$programme_key, 'magnitude'),
        likelihood_column(assessment$programme_key, 'occurrence'),
        likelihood_column(assessment$programme_key, 'magnitude')
    )

    beta_hyper <- list(
        theta = list(prior = 'loggamma', param = c(2, 0.1))
    )
    family_control <- c(
        rep(list(list()), 3L),
        rep(list(list(hyper = beta_hyper)), 3L)
    )

    set.seed(seed)
    started <- proc.time()[['elapsed']]
    fit <- inla(
        make_formula(candidate),
        data = stack_data,
        family = families,
        Ntrials = 1,
        control.family = family_control,
        control.fixed = list(mean = 0, prec = 4),
        control.predictor = list(
            A = inla.stack.A(stack), compute = TRUE, link = link
        ),
        control.compute = list(
            waic = compute_criteria, dic = compute_criteria,
            cpo = compute_criteria, config = TRUE
        ),
        verbose = FALSE
    )
    elapsed <- proc.time()[['elapsed']] - started

    occurrence_index <- inla.stack.index(
        stack, 'prediction_occurrence'
    )$data
    magnitude_index <- inla.stack.index(
        stack, 'prediction_magnitude'
    )$data
    occurrence <- fit$summary.fitted.values[occurrence_index, 'mean']
    magnitude <- fit$summary.fitted.values[magnitude_index, 'mean']

    predictions <- assessment |>
        transmute(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, ann_maxdhw_original, applied_dhw_uplift,
            ann_maxdhw,
            observed_mortality = mortality_prop,
            observed_occurrence = as.numeric(mortality_prop > 0),
            predicted_occurrence = occurrence,
            predicted_positive_mortality = magnitude,
            predicted_mortality = occurrence * magnitude
        )

    list(
        fit = fit, predictions = predictions,
        preprocessing = prepared$preprocessing,
        elapsed_seconds = elapsed,
        analysis_rows = nrow(analysis),
        assessment_rows = nrow(assessment)
    )
}

metric_summary <- function(rows) {
    rows |>
        summarise(
            n = n(), severe_n = sum(observed_mortality >= 0.50),
            rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
            mae = mean(abs(observed_mortality - predicted_mortality)),
            predictive_r2 = 1 -
                sum((observed_mortality - predicted_mortality)^2) /
                sum((observed_mortality - mean(observed_mortality))^2),
            bias = mean(predicted_mortality - observed_mortality),
            occurrence_brier = mean(
                (observed_occurrence - predicted_occurrence)^2
            ),
            severe_observed = if_else(
                severe_n > 0,
                mean(observed_mortality[observed_mortality >= 0.50]),
                NA_real_
            ),
            severe_predicted = if_else(
                severe_n > 0,
                mean(predicted_mortality[observed_mortality >= 0.50]),
                NA_real_
            ),
            severe_rmse = if_else(
                severe_n > 0,
                sqrt(mean((observed_mortality[observed_mortality >= 0.50] -
                               predicted_mortality[observed_mortality >= 0.50])^2)),
                NA_real_
            ),
            false_extreme_rate = mean(
                predicted_mortality[observed_mortality < 0.30] >= 0.30
            ),
            .groups = 'drop'
        )
}

extract_fit_summary <- function(result, candidate) {
    fit <- result$fit
    criteria <- tibble(
        candidate = candidate$candidate,
        waic = if (!is.null(fit$waic$waic)) fit$waic$waic else NA_real_,
        dic = if (!is.null(fit$dic$dic)) fit$dic$dic else NA_real_,
        log_marginal_likelihood = if (!is.null(fit$mlik)) {
            fit$mlik[1, 1]
        } else NA_real_,
        mean_neg_log_cpo = if (!is.null(fit$cpo$cpo)) {
            mean(-log(fit$cpo$cpo[is.finite(fit$cpo$cpo)]))
        } else NA_real_,
        elapsed_seconds = result$elapsed_seconds
    )
    fixed <- fit$summary.fixed |>
        as_tibble(rownames = 'term') |>
        mutate(candidate = candidate$candidate, .before = 1)
    hyper <- fit$summary.hyperpar |>
        as_tibble(rownames = 'term') |>
        mutate(candidate = candidate$candidate, .before = 1)
    list(criteria = criteria, fixed = fixed, hyper = hyper)
}

if (identical(Sys.getenv('INLA_ST_RUN'), '1')) {
    requested_candidates <- trimws(strsplit(
        Sys.getenv('INLA_ST_CANDIDATES', paste(candidates$candidate, collapse = ',')),
        ',', fixed = TRUE
    )[[1]])
    requested_schemes <- trimws(strsplit(
        Sys.getenv(
            'INLA_ST_SCHEMES',
            'reef_blocked_5fold,leave_one_event_out'
        ), ',', fixed = TRUE
    )[[1]])
    use_candidates <- candidates |>
        filter(candidate %in% requested_candidates)
    if (nrow(use_candidates) == 0L) stop('No requested INLA candidates')

    full_criteria <- tibble()
    full_fixed <- tibble()
    full_hyper <- tibble()
    temporal_slopes <- tibble()
    spatial_nodes <- tibble()
    partial_pool_deviations <- tibble()
    predictions <- tibble()
    runtime <- tibble()

    for (i in seq_len(nrow(use_candidates))) {
        candidate <- use_candidates[i, ]
        cache_file <- file.path(
            output_dir,
            paste0(
                'full_', candidate_cache_version(candidate$candidate), '_',
                candidate$candidate, '.rds'
            )
        )
        if (file.exists(cache_file)) {
            result <- readRDS(cache_file)
            message('Loaded full INLA: ', candidate$candidate)
        } else {
            result <- fit_candidate(
                data, data[1, , drop = FALSE], candidate,
                compute_criteria = TRUE, seed = 20260827L + i
            )
            saveRDS(result, cache_file)
            message('Fitted full INLA: ', candidate$candidate)
        }
        summary_parts <- extract_fit_summary(result, candidate)
        full_criteria <- bind_rows(full_criteria, summary_parts$criteria)
        full_fixed <- bind_rows(full_fixed, summary_parts$fixed)
        full_hyper <- bind_rows(full_hyper, summary_parts$hyper)
        runtime <- bind_rows(runtime, tibble(
            candidate = candidate$candidate, scheme = 'full', fold = 'all',
            analysis_rows = nrow(data), assessment_rows = 0L,
            elapsed_seconds = result$elapsed_seconds
        ))

        if ('event_dhw_index' %in% names(result$fit$summary.random)) {
            temporal_slopes <- bind_rows(
                temporal_slopes,
                result$fit$summary.random$event_dhw_index |>
                    as_tibble() |>
                    mutate(
                        candidate = candidate$candidate,
                        event_year = event_years[ID], .before = 1
                    )
            )
        }
        if ('spatial_field' %in% names(result$fit$summary.random)) {
            spatial_nodes <- bind_rows(
                spatial_nodes,
                result$fit$summary.random$spatial_field |>
                    as_tibble() |>
                    mutate(
                        candidate = candidate$candidate,
                        event_year = NA_integer_,
                        # INLA reports SPDE node IDs from zero on this build.
                        # Match by row order rather than treating IDs as R indices.
                        node = seq_len(n()),
                        x_km = mesh$loc[node, 1], y_km = mesh$loc[node, 2],
                        .before = 1
                    )
            )
        }
        if ('event_spatial_field' %in% names(result$fit$summary.random)) {
            field <- result$fit$summary.random$event_spatial_field |>
                as_tibble() |>
                mutate(
                    node = rep(seq_len(mesh$n), times = length(event_years)),
                    event_year = rep(event_years, each = mesh$n),
                    candidate = candidate$candidate,
                    x_km = mesh$loc[node, 1], y_km = mesh$loc[node, 2],
                    .before = 1
                )
            spatial_nodes <- bind_rows(spatial_nodes, field)
        }
        pooled_effects <- intersect(
            names(result$fit$summary.random),
            c(
                'pool_acropora_index', 'pool_novelty_index',
                'pool_freshwater_index', 'pool_cloud_index',
                'pool_current_index'
            )
        )
        for (pooled_effect in pooled_effects) {
            partial_pool_deviations <- bind_rows(
                partial_pool_deviations,
                result$fit$summary.random[[pooled_effect]] |>
                    as_tibble() |>
                    left_join(layer_lookup, by = 'ID') |>
                    mutate(
                        candidate = candidate$candidate,
                        interaction = gsub(
                            '^pool_|_index$', '', pooled_effect
                        ),
                        .before = 1
                    )
            )
        }
    }

    for (scheme in requested_schemes) {
        folds <- if (scheme == 'reef_blocked_5fold') {
            sort(unique(data$joint_reef_fold))
        } else if (scheme == 'leave_one_event_out') {
            event_years
        } else {
            stop('Unknown validation scheme: ', scheme)
        }
        for (i in seq_len(nrow(use_candidates))) {
            candidate <- use_candidates[i, ]
            for (fold in folds) {
                if (scheme == 'reef_blocked_5fold') {
                    held_reefs <- unique(
                        data$ReefID[data$joint_reef_fold == fold]
                    )
                    assessment <- data |>
                        filter(ReefID %in% held_reefs)
                    analysis <- data |>
                        filter(!ReefID %in% held_reefs)
                } else {
                    assessment <- data |>
                        filter(event_year == fold)
                    analysis <- data |>
                        filter(event_year != fold)
                }
                cache_file <- file.path(
                    output_dir,
                    paste0(
                        'cv_', candidate_cache_version(candidate$candidate),
                        '_', candidate$candidate, '_', scheme, '_',
                        fold, '.rds'
                    )
                )
                if (file.exists(cache_file)) {
                    cached <- readRDS(cache_file)
                    message(
                        'Loaded CV INLA: ', candidate$candidate, ' / ',
                        scheme, ' / ', fold
                    )
                } else {
                    result <- fit_candidate(
                        analysis, assessment, candidate,
                        compute_criteria = FALSE,
                        seed = 20260827L + i * 100L + as.integer(fold)
                    )
                    cached <- list(
                        predictions = result$predictions,
                        elapsed_seconds = result$elapsed_seconds,
                        analysis_rows = result$analysis_rows,
                        assessment_rows = result$assessment_rows
                    )
                    saveRDS(cached, cache_file)
                    message(
                        'Fitted CV INLA: ', candidate$candidate, ' / ',
                        scheme, ' / ', fold
                    )
                }
                predictions <- bind_rows(
                    predictions,
                    cached$predictions |>
                        mutate(
                            candidate = candidate$candidate,
                            scheme = scheme, fold = as.character(fold)
                        )
                )
                runtime <- bind_rows(runtime, tibble(
                    candidate = candidate$candidate, scheme = scheme,
                    fold = as.character(fold),
                    analysis_rows = cached$analysis_rows,
                    assessment_rows = cached$assessment_rows,
                    elapsed_seconds = cached$elapsed_seconds
                ))
            }
        }
    }

    metrics <- bind_rows(
        predictions |>
            group_by(candidate, scheme) |>
            metric_summary() |>
            mutate(programme_key = 'all'),
        predictions |>
            group_by(candidate, scheme, programme_key) |>
            metric_summary()
    )
    event_metrics <- predictions |>
        group_by(candidate, scheme, programme_key, event_year) |>
        metric_summary()

    adaptation_terms <- full_fixed |>
        filter(term %in% c(
            'dhw_events_since2016_n6_z',
            'dhw_years_since_last_n6_capped8_z',
            'dhw_x_repeat_exposure', 'dhw_x_recovery_interval',
            'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z'
        ))

    model_comparison <- metrics |>
        filter(programme_key == 'all') |>
        select(
            candidate, scheme, n, rmse, mae, predictive_r2, bias,
            occurrence_brier, severe_predicted, severe_observed,
            severe_rmse, false_extreme_rate
        ) |>
        left_join(full_criteria, by = 'candidate')

    expected_predictions <- nrow(data) * nrow(use_candidates) *
        length(requested_schemes)
    stopifnot(
        nrow(predictions) == expected_predictions,
        all(is.finite(predictions$predicted_mortality)),
        all(predictions$predicted_mortality >= 0 &
                predictions$predicted_mortality <= 1),
        all(predictions$predicted_occurrence >= 0 &
                predictions$predicted_occurrence <= 1),
        all(predictions$predicted_positive_mortality >= 0 &
                predictions$predicted_positive_mortality <= 1)
    )

    mesh_nodes <- tibble(
        node = seq_len(mesh$n),
        x_km = mesh$loc[, 1], y_km = mesh$loc[, 2]
    )
    write_csv(candidates, file.path(output_dir, 'candidate_definitions.csv'))
    write_csv(full_criteria, file.path(output_dir, 'full_fit_criteria.csv'))
    write_csv(full_fixed, file.path(output_dir, 'full_fixed_effects.csv'))
    write_csv(full_hyper, file.path(output_dir, 'full_hyperparameters.csv'))
    write_csv(adaptation_terms, file.path(output_dir, 'adaptation_terms.csv'))
    write_csv(temporal_slopes, file.path(output_dir, 'temporal_dhw_slopes.csv'))
    write_csv(spatial_nodes, file.path(output_dir, 'spatial_field_nodes.csv'))
    write_csv(
        partial_pool_deviations,
        file.path(output_dir, 'partial_pool_deviations.csv')
    )
    write_csv(mesh_nodes, file.path(output_dir, 'mesh_nodes.csv'))
    write_csv(predictions, file.path(output_dir, 'cv_predictions.csv'))
    write_csv(metrics, file.path(output_dir, 'cv_metrics.csv'))
    write_csv(event_metrics, file.path(output_dir, 'event_metrics.csv'))
    write_csv(model_comparison, file.path(output_dir, 'model_comparison.csv'))
    write_csv(runtime, file.path(output_dir, 'runtime.csv'))

    print(model_comparison)
}
