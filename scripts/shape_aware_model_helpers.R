# Shared definitions for the shape-aware mortality candidates. These models are
# deliberately separate from the formal benchmark so the existing results and
# caches remain reproducible.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
})
source('scripts/formal_model_helpers.R')

shape_environmental_predictors <- c(
    'prop_acropora_pre', 'observed_pre_cover',
    'dhw_excess4', 'dhw_excess8',
    'histmDHW6', 'yrsince6', 'histmDHW4', 'yrsince4',
    'ann_maxsst', 'winyear_mean', 'winyear_sd', 'mcur_90',
    'secc3m', 'cloudp_90', 'dhw10_load4', 'dhw_novelty10',
    'secc3m_p10', 'dhw_events_since2016_n6',
    'dhw_years_since_last_n6_capped8', 'dhw_no_prior_n6'
)

shape_predictors <- function(programme_key) {
    if (programme_key == 'mmp') {
        c(shape_environmental_predictors, 'depth')
    } else {
        shape_environmental_predictors
    }
}

add_shape_features <- function(data) {
    data |>
        mutate(
            dhw_excess4 = pmax(ann_maxdhw - 4, 0),
            dhw_excess8 = pmax(ann_maxdhw - 8, 0)
        )
}

load_shape_rows <- function(programme_key) {
    load_programme_rows(programme_key) |>
        add_shape_features()
}

shape_mean_rhs <- function(predictors) {
    scaled_main_effects <- paste0(predictors, '_z')
    thermal_composition <- c(
        'dhw_excess4_z:prop_acropora_pre',
        'dhw_excess8_z:prop_acropora_pre',
        'dhw_excess4_z:dhw_events_since2016_n6_z',
        'dhw_excess8_z:dhw_events_since2016_n6_z',
        'dhw_excess4_z:dhw_years_since_last_n6_capped8_z',
        'dhw_excess8_z:dhw_years_since_last_n6_capped8_z'
    )
    paste(
        c(
            scaled_main_effects,
            thermal_composition,
            '(1 | reef_effect)',
            '(1 | region_effect)'
        ),
        collapse = ' + '
    )
}

shape_boundary_rhs <- function() {
    paste(
        c(
            'dhw_excess4_z', 'dhw_excess8_z',
            'prop_acropora_pre_z', 'observed_pre_cover_z',
            'cloudp_90_z', 'secc3m_z', 'histmDHW6_z',
            'dhw10_load4_z', 'dhw_novelty10_z', 'secc3m_p10_z',
            'dhw_events_since2016_n6_z',
            'dhw_years_since_last_n6_capped8_z', 'dhw_no_prior_n6_z',
            'dhw_excess4_z:prop_acropora_pre',
            'dhw_excess8_z:prop_acropora_pre',
            'dhw_excess4_z:dhw_events_since2016_n6_z',
            'dhw_excess8_z:dhw_events_since2016_n6_z',
            'dhw_excess4_z:dhw_years_since_last_n6_capped8_z',
            'dhw_excess8_z:dhw_years_since_last_n6_capped8_z',
            '(1 | region_effect)'
        ),
        collapse = ' + '
    )
}

make_shape_formula <- function(programme_key, predictors) {
    mean_formula <- as.formula(
        paste('mortality_prop ~', shape_mean_rhs(predictors))
    )
    boundary_formula <- as.formula(paste('~', shape_boundary_rhs()))

    if (programme_key == 'mmp') {
        bf(mean_formula, zoi = boundary_formula, coi = ~ 1)
    } else {
        bf(mean_formula, zi = boundary_formula)
    }
}

shape_family <- function(programme_key) {
    if (programme_key == 'mmp') {
        zero_one_inflated_beta()
    } else {
        zero_inflated_beta()
    }
}

shape_priors <- function(programme_key) {
    boundary_parameter <- if (programme_key == 'mmp') 'zoi' else 'zi'
    boundary_intercept <- if (programme_key == 'mmp') -0.5 else 0
    priors <- c(
        prior(normal(0, 0.5), class = 'b'),
        prior(normal(-1.5, 1), class = 'Intercept'),
        prior(exponential(1), class = 'sd'),
        prior(gamma(2, 0.1), class = 'phi'),

        # Ecological shape priors strongly favour a rising response and a
        # steeper slope above 8 DHW without hard-coding a 50% outcome.
        prior(
            normal(0.3, 0.25), class = 'b', coef = 'dhw_excess4_z'
        ),
        prior(
            normal(0.6, 0.35), class = 'b', coef = 'dhw_excess8_z'
        ),
        prior(
            normal(0.25, 0.25), class = 'b',
            coef = 'dhw_excess4_z:prop_acropora_pre'
        ),
        prior(
            normal(0.6, 0.4), class = 'b',
            coef = 'dhw_excess8_z:prop_acropora_pre'
        ),

        set_prior(
            'normal(0, 0.7)', class = 'b', dpar = boundary_parameter
        ),
        set_prior(
            paste0('normal(', boundary_intercept, ', 1)'),
            class = 'Intercept', dpar = boundary_parameter
        ),

        # In the zero-inflation parameterisation, negative coefficients make
        # mortality occurrence increasingly likely as heat increases.
        set_prior(
            'normal(-0.3, 0.25)', class = 'b', coef = 'dhw_excess4_z',
            dpar = boundary_parameter
        ),
        set_prior(
            'normal(-0.6, 0.35)', class = 'b', coef = 'dhw_excess8_z',
            dpar = boundary_parameter
        ),
        set_prior(
            'normal(-0.25, 0.25)', class = 'b',
            coef = 'dhw_excess4_z:prop_acropora_pre',
            dpar = boundary_parameter
        ),
        set_prior(
            'normal(-0.6, 0.4)', class = 'b',
            coef = 'dhw_excess8_z:prop_acropora_pre',
            dpar = boundary_parameter
        )
    )

    if (programme_key == 'mmp') {
        c(
            priors,
            prior(normal(-2.7, 1), class = 'Intercept', dpar = 'coi')
        )
    } else {
        priors
    }
}

shape_prediction_components <- function(fit, assessment, programme_key) {
    draws <- posterior_epred(
        fit, newdata = assessment, re_formula = NA,
        allow_new_levels = TRUE
    )
    if (programme_key == 'mmp') {
        boundary <- posterior_linpred(
            fit, newdata = assessment, dpar = 'zoi', transform = TRUE,
            re_formula = NA, allow_new_levels = TRUE
        )
        complete <- posterior_linpred(
            fit, newdata = assessment, dpar = 'coi', transform = TRUE,
            re_formula = NA, allow_new_levels = TRUE
        )
        occurrence <- 1 - boundary * (1 - complete)
    } else {
        zero_probability <- posterior_linpred(
            fit, newdata = assessment, dpar = 'zi', transform = TRUE,
            re_formula = NA, allow_new_levels = TRUE
        )
        occurrence <- 1 - zero_probability
    }
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
