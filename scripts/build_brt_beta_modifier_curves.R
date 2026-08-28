# Compare the production BRT and Bayesian beta/Bernoulli-boundary models on a
# common set of population-level DHW profiles. Each ecological modifier is
# shown at its programme-specific 10th, 50th and 90th percentile while all
# other predictors are held at their programme median.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(gbm)
    library(readr)
    library(tidyr)
})
source('scripts/formal_model_helpers.R')

set.seed(20260826L)
model_dir <- 'output/models/formal'
output_dir <- 'output/zero_mortality_sensitivity'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

programmes <- c('manta', 'ltmp', 'mmp')
dhw_grid <- seq(0, 16, by = 0.25)

# These three scenarios are the strongest non-redundant, operationally
# available modifiers shared by the two production learners. Acropora is
# represented coherently by varying both its proportion and absolute cover.
modifier_definitions <- tibble(
    modifier = c('acropora_composition', 'dhw_novelty10', 'cloudp_90'),
    modifier_label = c(
        'Pre-event Acropora composition',
        'Thermal novelty relative to prior decade',
        'Peak-season cloud fraction'
    )
)

apply_preprocessing <- function(rows, rules) {
    for (i in seq_len(nrow(rules))) {
        name <- rules$predictor[[i]]
        rows[[name]][!is.finite(rows[[name]])] <- rules$median[[i]]
        rows[[paste0(name, '_z')]] <- (
            rows[[name]] - rules$mean[[i]]
        ) / rules$sd[[i]]
    }
    rows
}

predict_brt_component <- function(component, newdata, n_trees, type) {
    if (component$type == 'constant') {
        return(rep(component$value, nrow(newdata)))
    }
    predict(
        component$model, newdata = newdata,
        n.trees = n_trees, type = type
    )
}

predict_brt <- function(cache, newdata) {
    trees <- as.integer(cache$tuning$best$n_trees[[1]])
    occurrence <- predict_brt_component(
        cache$model$occurrence, newdata, trees, 'response'
    )
    magnitude <- predict_brt_component(
        cache$model$magnitude, newdata, trees, 'response'
    )
    if (cache$model$magnitude_mode == 'logit_gaussian') {
        magnitude <- plogis(magnitude)
    }
    occurrence <- pmin(pmax(occurrence, 0), 1)
    magnitude <- pmin(pmax(magnitude, 0), 1)
    tibble(
        occurrence = occurrence,
        positive_mortality = magnitude,
        mortality = occurrence * magnitude,
        q05 = NA_real_, q95 = NA_real_
    )
}

predict_brms <- function(cache, newdata, programme) {
    fit <- cache$fit
    expected <- posterior_epred(
        fit, newdata = newdata, re_formula = NA,
        allow_new_levels = TRUE
    )
    if (programme == 'mmp') {
        boundary <- posterior_linpred(
            fit, newdata = newdata, dpar = 'zoi', transform = TRUE,
            re_formula = NA, allow_new_levels = TRUE
        )
        complete <- posterior_linpred(
            fit, newdata = newdata, dpar = 'coi', transform = TRUE,
            re_formula = NA, allow_new_levels = TRUE
        )
        occurrence_draws <- 1 - boundary * (1 - complete)
    } else {
        zero_draws <- posterior_linpred(
            fit, newdata = newdata, dpar = 'zi', transform = TRUE,
            re_formula = NA, allow_new_levels = TRUE
        )
        occurrence_draws <- 1 - zero_draws
    }
    positive_draws <- expected / pmax(occurrence_draws, 1e-8)
    tibble(
        occurrence = colMeans(occurrence_draws),
        positive_mortality = pmin(pmax(colMeans(positive_draws), 0), 1),
        mortality = colMeans(expected),
        q05 = apply(expected, 2, quantile, 0.05),
        q95 = apply(expected, 2, quantile, 0.95)
    )
}

make_profile_grid <- function(data, preprocessing, modifier, programme) {
    probabilities <- c(0.10, 0.50, 0.90)
    level_labels <- c('Low (10th percentile)', 'Median (50th percentile)',
                      'High (90th percentile)')
    if (modifier == 'acropora_composition') {
        values <- quantile(
            data$prop_acropora_pre[is.finite(data$prop_acropora_pre)],
            probabilities, names = FALSE
        )
        units <- 'proportion of pre-event hard coral'
    } else {
        values <- quantile(
            data[[modifier]][is.finite(data[[modifier]])],
            probabilities, names = FALSE
        )
        units <- if (modifier == 'cloudp_90') 'fraction' else 'DHW anomaly'
    }

    profile <- expand_grid(
        ann_maxdhw = dhw_grid,
        modifier_level = factor(level_labels, levels = level_labels)
    ) |>
        mutate(
            modifier_value = values[match(modifier_level, level_labels)],
            modifier_units = units
        )
    rows <- data[rep(1, nrow(profile)), , drop = FALSE]
    for (predictor in preprocessing$predictor) {
        rows[[predictor]] <- preprocessing$median[
            match(predictor, preprocessing$predictor)
        ]
    }
    rows$ann_maxdhw <- profile$ann_maxdhw
    if (modifier == 'acropora_composition') {
        rows$prop_acropora_pre <- profile$modifier_value
        rows$acropora_cover_pre <- profile$modifier_value *
            median(data$observed_pre_cover, na.rm = TRUE)
    } else {
        rows[[modifier]] <- profile$modifier_value
    }
    rows <- apply_preprocessing(rows, preprocessing)
    bind_cols(
        profile,
        rows |>
            select(-ann_maxdhw) |>
            mutate(
                programme_key = .env$programme,
                modifier = .env$modifier
            )
    )
}

all_curves <- tibble()
scenario_values <- tibble()
for (programme in programmes) {
    data <- load_programme_rows(programme)
    brt_cache <- readRDS(file.path(
        model_dir, paste0('brt_', programme, '_production.rds')
    ))
    brms_cache <- readRDS(file.path(
        model_dir, paste0('brms_', programme, '_production.rds')
    ))
    observed_min <- min(data$ann_maxdhw)
    observed_max <- max(data$ann_maxdhw)

    for (modifier in modifier_definitions$modifier) {
        brt_grid <- make_profile_grid(
            data, brt_cache$preprocessing, modifier, programme
        )
        brms_grid <- make_profile_grid(
            data, brms_cache$preprocessing, modifier, programme
        )
        profile_columns <- c(
            'programme_key', 'modifier', 'modifier_level',
            'modifier_value', 'modifier_units', 'ann_maxdhw'
        )

        brt_prediction <- predict_brt(brt_cache, brt_grid)
        brms_prediction <- predict_brms(brms_cache, brms_grid, programme)
        new_curves <- bind_rows(
            brt_grid |>
                select(all_of(profile_columns)) |>
                bind_cols(brt_prediction) |>
                mutate(learner = 'BRT'),
            brms_grid |>
                select(all_of(profile_columns)) |>
                bind_cols(brms_prediction) |>
                mutate(learner = 'Bayesian beta + Bernoulli boundary')
        ) |>
            mutate(
                observed_dhw_min = observed_min,
                observed_dhw_max = observed_max,
                outside_observed_dhw = ann_maxdhw < observed_min |
                    ann_maxdhw > observed_max
            )
        all_curves <- bind_rows(all_curves, new_curves)

        scenario_values <- bind_rows(
            scenario_values,
            brt_grid |>
                distinct(
                    programme_key, modifier, modifier_level,
                    modifier_value, modifier_units
                )
        )
    }
    message('Completed modifier profiles for ', programme)
}

all_curves <- all_curves |>
    left_join(modifier_definitions, by = 'modifier') |>
    relocate(modifier_label, .after = modifier)

# Rank candidate modifiers using both BRT relative influence and the largest
# standardised Bayesian fixed effect containing that predictor. The Acropora
# scenario combines its proportion and absolute-cover representations.
brt_influence <- read_csv(
    'output/formal_models/brt_production_influence.csv',
    show_col_types = FALSE
)
brms_effects <- read_csv(
    'output/formal_models/brms_production_coefficients.csv',
    show_col_types = FALSE
)
candidate_groups <- list(
    acropora_composition = c('prop_acropora_pre', 'acropora_cover_pre'),
    dhw_novelty10 = 'dhw_novelty10',
    cloudp_90 = 'cloudp_90',
    mcur_90 = 'mcur_90',
    histmDHW6 = 'histmDHW6',
    secc3m = 'secc3m',
    winyear_mean = 'winyear_mean',
    winyear_sd = 'winyear_sd'
)
selection_audit <- bind_rows(lapply(names(candidate_groups), function(group) {
    terms <- candidate_groups[[group]]
    brt_value <- brt_influence |>
        filter(predictor %in% terms) |>
        group_by(programme_key, component) |>
        summarise(value = sum(relative_influence), .groups = 'drop') |>
        summarise(value = mean(value)) |>
        pull(value)
    pattern <- paste(terms, collapse = '|')
    brms_value <- brms_effects |>
        filter(grepl(pattern, term)) |>
        group_by(programme_key) |>
        summarise(value = max(abs(estimate)), .groups = 'drop') |>
        summarise(value = mean(value)) |>
        pull(value)
    tibble(
        modifier = group,
        mean_brt_relative_influence = brt_value,
        mean_max_abs_standardised_beta = brms_value
    )
})) |>
    mutate(
        brt_rank = min_rank(desc(mean_brt_relative_influence)),
        bayesian_rank = min_rank(desc(mean_max_abs_standardised_beta)),
        combined_rank_score = brt_rank + bayesian_rank,
        selected = modifier %in% modifier_definitions$modifier
    ) |>
    arrange(combined_rank_score, brt_rank)

curve_landmarks <- all_curves |>
    filter(ann_maxdhw %in% c(4, 8, 12)) |>
    select(
        programme_key, learner, modifier, modifier_label,
        modifier_level, modifier_value, ann_maxdhw,
        occurrence, positive_mortality, mortality,
        outside_observed_dhw
    )

write_csv(all_curves, file.path(output_dir, 'modifier_dhw_curves.csv'), na = '')
write_csv(scenario_values, file.path(output_dir, 'modifier_scenario_values.csv'), na = '')
write_csv(selection_audit, file.path(output_dir, 'modifier_selection_audit.csv'), na = '')
write_csv(curve_landmarks, file.path(output_dir, 'modifier_curve_landmarks.csv'), na = '')

print(selection_audit)
print(curve_landmarks)
