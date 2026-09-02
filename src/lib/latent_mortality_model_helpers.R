# Shared latent ecological mortality model.
#
# One ecological state is estimated for each reef-event. Environmental and
# community predictors enter that state once. Manta tow is the reference
# observation scale; LTMP and MMP receive separate bias, precision and zero/one
# observation parameters. A common response loading is used because MMP has
# insufficient paired overlap to identify its own ecological scale reliably.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
})
source('src/lib/joint_compound_model_helpers.R')

latent_mortality_predictors <- c(
    'ann_maxdhw', 'dhw_excess4', 'dhw_excess8',
    'prop_acropora_pre', 'observed_pre_cover',
    'histmDHW6', 'yrsince6', 'dhw10_load4', 'dhw_novelty10',
    'dhw_events_since2016_n6', 'dhw_years_since_last_n6_capped8',
    'dhw_no_prior_n6',
    'secc3m_p10', 'cloudp_90',
    'log_coastal_rain30', 'era5_wind_mean',
    'era5_wind_calm_fraction', 'era5_coastal_distance_km',
    'wqc_freqcc12', 'log1p_cyc_maxHrs4mw'
)

median_finite <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else median(x)
}

make_latent_event_rows <- function(observations) {
    observations |>
        group_by(ReefID, event_year) |>
        summarise(
            ReefName = first(ReefName),
            region_block = first(region_block),
            across(all_of(latent_mortality_predictors), median_finite),
            programmes_observed = paste(
                sort(unique(programme_key)), collapse = '+'
            ),
            observations = n(),
            .groups = 'drop'
        ) |>
        mutate(reef_event_effect = factor(paste(ReefID, event_year, sep = '__')))
}

prepare_latent_fold <- function(analysis_observations, assessment_observations) {
    analysis_events <- make_latent_event_rows(analysis_observations)
    assessment_events <- make_latent_event_rows(assessment_observations)
    prepared <- prepare_fold_predictors(
        analysis_events, assessment_events, latent_mortality_predictors
    )
    event_columns <- c(
        'ReefID', 'event_year', 'reef_event_effect',
        paste0(latent_mortality_predictors, '_z')
    )

    attach_events <- function(observations, events) {
        observations |>
            select(-any_of(c(
                'reef_event_effect',
                paste0(latent_mortality_predictors, '_z')
            ))) |>
            left_join(
                select(events, all_of(event_columns)),
                by = c('ReefID', 'event_year'),
                relationship = 'many-to-one'
            ) |>
            mutate(
                programme_factor = factor(
                    programme_key, levels = joint_programmes
                ),
                is_ltmp = as.numeric(programme_key == 'ltmp'),
                is_mmp = as.numeric(programme_key == 'mmp'),
                reef_event_effect = factor(reef_event_effect)
            )
    }

    list(
        analysis = attach_events(analysis_observations, prepared$analysis),
        assessment = attach_events(
            assessment_observations, prepared$assessment
        ),
        analysis_events = prepared$analysis,
        assessment_events = prepared$assessment,
        preprocessing = prepared$preprocessing
    )
}

latent_ecological_rhs <- function() {
    interactions <- c(
        'dhw_excess4_z:prop_acropora_pre_z',
        'dhw_excess8_z:prop_acropora_pre_z',
        'ann_maxdhw_z:dhw_events_since2016_n6_z',
        'ann_maxdhw_z:dhw_years_since_last_n6_capped8_z',
        'ann_maxdhw_z:log_coastal_rain30_z',
        'ann_maxdhw_z:wqc_freqcc12_z',
        'ann_maxdhw_z:era5_wind_calm_fraction_z'
    )
    paste(
        c(
            '1', paste0(latent_mortality_predictors, '_z'),
            interactions, '(1 | reef_event_effect)'
        ),
        collapse = ' + '
    )
}

make_latent_mortality_formula <- function() {
    bf(
        # The beta-family logit link transforms this predictor to 0--1.
        # Do not apply inv_logit here or the mean is transformed twice.
        mortality_prop ~
            biasltmp * is_ltmp + biasmmp * is_mmp +
                latenteta,
        latenteta = as.formula(paste('~', latent_ecological_rhs())),
        biasltmp ~ 1,
        biasmmp ~ 1,
        phi ~ 0 + programme_factor,
        zoi ~ 0 + programme_factor,
        coi ~ 0 + programme_factor,
        nl = TRUE
    )
}

latent_mortality_priors <- function() {
    priors <- c(
        prior(normal(0, 0.5), class = 'b', nlpar = 'latenteta'),
        prior(
            normal(-2, 0.8), class = 'b', coef = 'Intercept',
            nlpar = 'latenteta'
        ),
        prior(exponential(1), class = 'sd', nlpar = 'latenteta'),
        prior(
            normal(0, 0.35), class = 'b', coef = 'Intercept',
            nlpar = 'biasltmp'
        ),
        prior(
            normal(0, 0.50), class = 'b', coef = 'Intercept',
            nlpar = 'biasmmp'
        ),
        prior(normal(2, 0.8), class = 'b', dpar = 'phi'),
        prior(normal(-2, 1), class = 'b', dpar = 'zoi'),
        prior(normal(-4, 1.2), class = 'b', dpar = 'coi')
    )
    priors
}

latent_prediction_components <- function(fit, assessment) {
    # Held-out reef-events have genuinely new random-effect levels. Sampling
    # those levels marginalises process heterogeneity for operational use.
    expected <- posterior_epred(
        fit, newdata = assessment, re_formula = NULL,
        allow_new_levels = TRUE, sample_new_levels = 'gaussian'
    )
    latent_eta <- posterior_linpred(
        fit, newdata = assessment, nlpar = 'latenteta',
        re_formula = NULL, allow_new_levels = TRUE,
        sample_new_levels = 'gaussian', transform = FALSE
    )
    latent_mortality <- plogis(latent_eta)

    tibble(
        predicted_mortality = colMeans(expected),
        prediction_q05 = apply(expected, 2, quantile, 0.05),
        prediction_q50 = apply(expected, 2, quantile, 0.50),
        prediction_q95 = apply(expected, 2, quantile, 0.95),
        latent_ecological_mortality = colMeans(latent_mortality),
        latent_q05 = apply(latent_mortality, 2, quantile, 0.05),
        latent_q50 = apply(latent_mortality, 2, quantile, 0.50),
        latent_q95 = apply(latent_mortality, 2, quantile, 0.95)
    )
}

predict_operational_latent_mortality <- function(
    fit, new_events, preprocessing, include_process_variation = TRUE
) {
    missing_predictors <- setdiff(
        latent_mortality_predictors, names(new_events)
    )
    if (length(missing_predictors) > 0L) {
        stop(
            'Operational rows are missing: ',
            paste(missing_predictors, collapse = ', ')
        )
    }
    if (!all(c('ReefID', 'event_year') %in% names(new_events))) {
        stop('Operational rows require ReefID and event_year')
    }

    prepared <- new_events
    for (predictor in latent_mortality_predictors) {
        rule <- preprocessing[preprocessing$predictor == predictor, ]
        if (nrow(rule) != 1L) stop('Missing preprocessing rule for ', predictor)
        values <- prepared[[predictor]]
        values[!is.finite(values)] <- rule$median[[1]]
        prepared[[paste0(predictor, '_z')]] <-
            (values - rule$mean[[1]]) / rule$sd[[1]]
    }
    prepared <- prepared |>
        mutate(
            reef_event_effect = factor(paste(ReefID, event_year, sep = '__')),
            programme_factor = factor('manta', levels = joint_programmes),
            is_ltmp = 0,
            is_mmp = 0
        )

    fixed_eta <- posterior_linpred(
        fit, newdata = prepared, nlpar = 'latenteta',
        re_formula = NA, allow_new_levels = TRUE, transform = FALSE
    )
    process_eta <- if (include_process_variation) {
        posterior_linpred(
            fit, newdata = prepared, nlpar = 'latenteta',
            re_formula = NULL, allow_new_levels = TRUE,
            sample_new_levels = 'gaussian', transform = FALSE
        )
    } else {
        fixed_eta
    }
    fixed_mortality <- plogis(fixed_eta)
    process_mortality <- plogis(process_eta)

    bind_cols(
        select(new_events, any_of(c('ReefID', 'ReefName', 'event_year'))),
        tibble(
            latent_population_mean = colMeans(fixed_mortality),
            latent_process_mean = colMeans(process_mortality),
            latent_process_q05 = apply(
                process_mortality, 2, quantile, 0.05
            ),
            latent_process_q50 = apply(
                process_mortality, 2, quantile, 0.50
            ),
            latent_process_q95 = apply(
                process_mortality, 2, quantile, 0.95
            )
        )
    )
}
