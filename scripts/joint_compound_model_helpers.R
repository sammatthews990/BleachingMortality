# Shared data and model definitions for the joint multi-program compound-stress
# candidate. Disturbance attribution is carried through unchanged and is not a
# predictor or an exclusion rule here.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(readr)
})
source('scripts/formal_model_helpers.R')

joint_programmes <- c('ltmp', 'manta', 'mmp')

joint_core_predictors <- c(
    'ann_maxdhw', 'dhw_excess4', 'dhw_excess8',
    'prop_acropora_pre', 'observed_pre_cover',
    'histmDHW6', 'yrsince6', 'dhw10_load4', 'dhw_novelty10',
    'secc3m_p10', 'cloudp_90', 'depth_within_programme'
)

joint_compound_predictors <- c(
    joint_core_predictors,
    'log_coastal_rain30', 'era5_wind_mean',
    'era5_wind_calm_fraction', 'era5_coastal_distance_km'
)

joint_rrn_wq_predictors <- c(joint_core_predictors, 'wqc_freqcc12')
joint_rrn_cyclone_predictors <- c(joint_core_predictors, 'cyc_maxHrs4mw')
joint_rrn_predictors <- c(
    joint_core_predictors, 'wqc_freqcc12', 'cyc_maxHrs4mw'
)
joint_rrn_relative_predictors <- c(
    joint_core_predictors,
    'wqc_prior10_percentile', 'log1p_cyc_maxHrs4mw'
)
joint_rrn_full_predictors <- c(
    joint_rrn_relative_predictors, 'log1p_cot_idwmeanpertow'
)
joint_compound_rrn_predictors <- c(
    joint_compound_predictors, 'wqc_freqcc12', 'cyc_maxHrs4mw'
)
joint_compound_rrn_relative_predictors <- c(
    joint_compound_predictors,
    'wqc_prior10_percentile', 'log1p_cyc_maxHrs4mw'
)
joint_compound_rrn_full_predictors <- c(
    joint_compound_rrn_relative_predictors,
    'log1p_cot_idwmeanpertow'
)
joint_all_predictors <- unique(c(
    joint_compound_rrn_predictors,
    joint_compound_rrn_full_predictors
))

joint_interactions <- c(
    'dhw_excess4_z:prop_acropora_pre_z',
    'dhw_excess8_z:prop_acropora_pre_z',
    'ann_maxdhw_z:log_coastal_rain30_z',
    'ann_maxdhw_z:era5_wind_calm_fraction_z'
)

assign_balanced_group_folds <- function(data, group, folds = 5L) {
    sizes <- data |>
        count(.data[[group]], name = 'observations') |>
        arrange(desc(observations), .data[[group]])
    number_of_folds <- min(as.integer(folds), nrow(sizes))
    fold_load <- rep(0L, number_of_folds)
    sizes$fold <- NA_integer_
    for (i in seq_len(nrow(sizes))) {
        chosen <- which.min(fold_load)
        sizes$fold[[i]] <- chosen
        fold_load[[chosen]] <- fold_load[[chosen]] + sizes$observations[[i]]
    }
    sizes$fold[match(data[[group]], sizes[[group]])]
}

load_joint_compound_rows <- function() {
    mortality <- bind_rows(lapply(joint_programmes, function(programme) {
        readRDS(validation_files[[programme]]) |>
            mutate(source_observation_id = as.character(source_observation_id))
    }))

    weather <- read_csv(
        'data/processed/era5_weather_reef_year.csv', show_col_types = FALSE
    ) |>
        group_by(year, grid_lat, grid_lon) |>
        summarise(
            log_coastal_rain30 = log1p(median(
                era5_coastal_rain_dec_mar_max_30day, na.rm = TRUE
            )),
            era5_wind_mean = median(
                era5_reef_wind_q1_mean, na.rm = TRUE
            ),
            era5_wind_calm_fraction = median(
                era5_reef_wind_q1_fraction_below_3, na.rm = TRUE
            ),
            era5_coastal_distance_km = median(
                coastal_grid_distance_km, na.rm = TRUE
            ),
            .groups = 'drop'
        )

    rrn <- read_csv(
        'data/processed/rrn_pressure_reef_year.csv', show_col_types = FALSE
    ) |>
        select(
            LABEL_ID, event_year, source_summer,
            wqc_freqcc12, wqc_prior10_n,
            wqc_prior10_percentile, wqc_prior10_delta,
            cyc_maxHrs4mw, log1p_cyc_maxHrs4mw,
            cot_meanpertow, cot_idwmeanpertow,
            log1p_cot_idwmeanpertow
        )

    joined <- mortality |>
        mutate(
            era5_grid_lat = floor(lat / 0.25 + 0.5) * 0.25,
            era5_grid_lon = floor(lon / 0.25 + 0.5) * 0.25
        ) |>
        left_join(
            weather,
            by = c(
                'event_year' = 'year',
                'era5_grid_lat' = 'grid_lat',
                'era5_grid_lon' = 'grid_lon'
            ),
            relationship = 'many-to-one'
        ) |>
        left_join(
            rrn,
            by = c('ReefID' = 'LABEL_ID', 'event_year'),
            relationship = 'many-to-one'
        ) |>
        group_by(programme_key) |>
        mutate(
            depth_within_programme = depth - median(depth, na.rm = TRUE)
        ) |>
        ungroup() |>
        mutate(
            dhw_excess4 = pmax(ann_maxdhw - 4, 0),
            dhw_excess8 = pmax(ann_maxdhw - 8, 0),
            programme_factor = factor(
                programme_key, levels = joint_programmes
            ),
            reef_effect = factor(ReefID),
            reef_event_effect = factor(paste(ReefID, event_year, sep = '__'))
        ) |>
        arrange(event_year, ReefID, programme_factor, source_observation_id)

    if (nrow(joined) != nrow(mortality)) {
        stop('ERA5 join changed the mortality row count')
    }
    missing_weather <- !complete.cases(
        joined[, joint_compound_predictors, drop = FALSE]
    )
    if (any(missing_weather)) {
        stop(sum(missing_weather), ' joint rows have incomplete predictors')
    }
    missing_rrn <- !complete.cases(
        joined[, c(
            'wqc_freqcc12', 'wqc_prior10_percentile',
            'cyc_maxHrs4mw', 'log1p_cyc_maxHrs4mw'
        ), drop = FALSE]
    )
    if (any(missing_rrn)) {
        stop(sum(missing_rrn), ' joint rows have incomplete RRN pressures')
    }
    missing_cots <- sum(!is.finite(joined$log1p_cot_idwmeanpertow))
    if (missing_cots > 0L) {
        message(
            missing_cots,
            ' joint rows lack modelled COTS pressure; ',
            'these are imputed within each analysis fold'
        )
    }
    joined$joint_reef_fold <- assign_balanced_group_folds(
        joined, 'ReefID', folds = 5L
    )
    joined
}

prepare_joint_predictors <- function(analysis, assessment) {
    prepare_fold_predictors(
        analysis, assessment, joint_all_predictors
    )
}

joint_mean_rhs <- function() {
    paste(
        c(
            '0 + programme_factor',
            paste0(joint_compound_predictors, '_z'),
            joint_interactions,
            '(1 | reef_event_effect)'
        ),
        collapse = ' + '
    )
}

joint_boundary_rhs <- function() {
    paste(
        c(
            '0 + programme_factor',
            'ann_maxdhw_z', 'dhw_excess4_z', 'dhw_excess8_z',
            'prop_acropora_pre_z', 'observed_pre_cover_z',
            'log_coastal_rain30_z', 'era5_wind_mean_z',
            'era5_wind_calm_fraction_z',
            'ann_maxdhw_z:log_coastal_rain30_z',
            'ann_maxdhw_z:era5_wind_calm_fraction_z'
        ),
        collapse = ' + '
    )
}

make_joint_compound_formula <- function() {
    bf(
        as.formula(paste('mortality_prop ~', joint_mean_rhs())),
        phi = ~ 0 + programme_factor,
        zoi = as.formula(paste('~', joint_boundary_rhs())),
        coi = ~ 0 + programme_factor
    )
}

joint_compound_priors <- function() {
    programme_coefficients <- paste0('programme_factor', joint_programmes)
    priors <- c(
        prior(normal(0, 0.5), class = 'b'),
        prior(exponential(1), class = 'sd'),
        prior(normal(0, 0.7), class = 'b', dpar = 'zoi'),
        prior(normal(2, 0.8), class = 'b', dpar = 'phi'),
        prior(normal(-4, 1.2), class = 'b', dpar = 'coi')
    )
    for (coefficient in programme_coefficients) {
        priors <- c(
            priors,
            set_prior(
                'normal(-2, 1)', class = 'b', coef = coefficient
            )
        )
    }
    c(
        priors,
        prior(normal(0.15, 0.2), class = 'b', coef = 'ann_maxdhw_z'),
        prior(normal(0.2, 0.25), class = 'b', coef = 'dhw_excess4_z'),
        prior(normal(0.35, 0.3), class = 'b', coef = 'dhw_excess8_z'),
        prior(
            normal(0.2, 0.25), class = 'b',
            coef = 'ann_maxdhw_z:log_coastal_rain30_z'
        ),
        prior(
            normal(0.1, 0.25), class = 'b',
            coef = 'ann_maxdhw_z:era5_wind_calm_fraction_z'
        ),
        prior(
            normal(-0.15, 0.25), class = 'b', coef = 'ann_maxdhw_z',
            dpar = 'zoi'
        ),
        prior(
            normal(-0.2, 0.3), class = 'b',
            coef = 'ann_maxdhw_z:log_coastal_rain30_z', dpar = 'zoi'
        ),
        prior(
            normal(-0.1, 0.3), class = 'b',
            coef = 'ann_maxdhw_z:era5_wind_calm_fraction_z', dpar = 'zoi'
        )
    )
}

joint_prediction_components <- function(fit, assessment) {
    draws <- posterior_epred(
        fit, newdata = assessment, re_formula = NA,
        allow_new_levels = TRUE
    )
    zoi <- posterior_linpred(
        fit, newdata = assessment, dpar = 'zoi', transform = TRUE,
        re_formula = NA, allow_new_levels = TRUE
    )
    coi <- posterior_linpred(
        fit, newdata = assessment, dpar = 'coi', transform = TRUE,
        re_formula = NA, allow_new_levels = TRUE
    )
    occurrence <- 1 - zoi * (1 - coi)
    positive <- draws / pmax(occurrence, 1e-8)
    tibble(
        predicted_occurrence = colMeans(occurrence),
        predicted_positive_loss = pmin(pmax(colMeans(positive), 0), 1),
        predicted_mortality = colMeans(draws),
        prediction_q05 = apply(draws, 2, quantile, 0.05),
        prediction_q50 = apply(draws, 2, quantile, 0.50),
        prediction_q95 = apply(draws, 2, quantile, 0.95)
    )
}

joint_assessment_index <- function(data, scheme, fold) {
    if (scheme == 'leave_one_event_out') {
        return(data$event_year == as.integer(fold))
    }
    if (scheme == 'reef_blocked_2024') {
        return(data$event_year == 2024L & data$joint_reef_fold == as.integer(fold))
    }
    stop('Unknown joint validation scheme: ', scheme)
}

joint_analysis_index <- function(data, scheme, fold) {
    assessment <- joint_assessment_index(data, scheme, fold)
    if (scheme == 'leave_one_event_out') return(!assessment)
    held_out_reefs <- unique(data$ReefID[assessment])
    !data$ReefID %in% held_out_reefs
}
