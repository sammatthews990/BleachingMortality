# Compare the shape-aware candidates with the untouched formal benchmarks and
# generate composition-specific production dose-response curves.

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(gbm)
    library(readr)
    library(tidyr)
})
source('scripts/shape_aware_model_helpers.R')

output_dir <- 'output/shape_aware_models'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required <- c(
    file.path(output_dir, 'shape_brms_event_predictions.csv'),
    file.path(output_dir, 'shape_brt_event_predictions.csv'),
    'output/formal_models/model_comparison_predictions.csv'
)
if (any(!file.exists(required))) {
    stop('Run both shape-aware fitting scripts before comparison.')
}

convex_weight <- function(observed, first_prediction, second_prediction) {
    difference <- first_prediction - second_prediction
    denominator <- sum(difference^2)
    if (!is.finite(denominator) || denominator < 1e-12) return(0.5)
    min(max(
        sum(difference * (observed - second_prediction)) / denominator,
        0
    ), 1)
}

shape_brms <- read_csv(required[[1]], show_col_types = FALSE) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )
shape_brt <- read_csv(required[[2]], show_col_types = FALSE) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    )
keys <- c(
    'programme_key', 'source_observation_id', 'ReefID', 'ReefName',
    'event_year', 'region_block', 'observed_mortality', 'scheme', 'fold',
    'ann_maxdhw', 'prop_acropora_pre', 'observed_pre_cover'
)
matched <- inner_join(
    shape_brms |>
        select(
            all_of(keys), brms_prediction = predicted_mortality,
            brms_occurrence = predicted_occurrence,
            brms_positive = predicted_positive_loss
        ),
    shape_brt |>
        select(
            all_of(keys), brt_prediction = predicted_mortality,
            brt_occurrence = predicted_occurrence,
            brt_positive = predicted_positive_loss
        ),
    by = keys
)
if (nrow(matched) != nrow(shape_brms) || nrow(matched) != nrow(shape_brt)) {
    stop('Shape-aware BRMS and BRT predictions do not match one-to-one.')
}

shape_ensemble <- tibble()
shape_fold_weights <- tibble()
for (programme in unique(matched$programme_key)) {
    programme_rows <- matched[matched$programme_key == programme, ]
    for (held_out_event in unique(programme_rows$fold)) {
        meta_rows <- programme_rows$fold != held_out_event
        assessment_rows <- !meta_rows
        weight_brms <- convex_weight(
            programme_rows$observed_mortality[meta_rows],
            programme_rows$brms_prediction[meta_rows],
            programme_rows$brt_prediction[meta_rows]
        )
        shape_ensemble <- bind_rows(
            shape_ensemble,
            programme_rows[assessment_rows, ] |>
                transmute(
                    across(all_of(keys)),
                    learner = 'shape_ensemble',
                    predicted_occurrence = weight_brms * brms_occurrence +
                        (1 - weight_brms) * brt_occurrence,
                    predicted_positive_loss = weight_brms * brms_positive +
                        (1 - weight_brms) * brt_positive,
                    predicted_mortality = weight_brms * brms_prediction +
                        (1 - weight_brms) * brt_prediction
                )
        )
        shape_fold_weights <- bind_rows(
            shape_fold_weights,
            tibble(
                programme_key = programme, fold = held_out_event,
                weight_brms, weight_brt = 1 - weight_brms,
                meta_analysis_rows = sum(meta_rows)
            )
        )
    }
}

shape_production_weights <- matched |>
    group_by(programme_key) |>
    group_modify(~ {
        weight_brms <- convex_weight(
            .x$observed_mortality, .x$brms_prediction, .x$brt_prediction
        )
        tibble(weight_brms, weight_brt = 1 - weight_brms, n = nrow(.x))
    }) |>
    ungroup()

context <- bind_rows(lapply(names(validation_files), function(programme) {
    readRDS(validation_files[[programme]]) |>
        transmute(
            programme_key,
            source_observation_id = as.character(source_observation_id),
            ann_maxdhw, prop_acropora_pre, observed_pre_cover
        )
}))
baseline <- read_csv(required[[3]], show_col_types = FALSE) |>
    mutate(
        source_observation_id = as.character(source_observation_id),
        fold = as.character(fold)
    ) |>
    filter(
        scheme == 'leave_one_event_out',
        learner %in% c('brms', 'brt', 'formal_ensemble')
    ) |>
    left_join(context, by = c('programme_key', 'source_observation_id')) |>
    mutate(learner = recode(
        learner,
        brms = 'benchmark_brms', brt = 'benchmark_brt',
        formal_ensemble = 'benchmark_ensemble'
    ))

comparison_rows <- bind_rows(
    baseline |>
        select(
            programme_key, source_observation_id, ReefID, ReefName,
            event_year, region_block, observed_mortality, scheme, fold,
            ann_maxdhw, prop_acropora_pre, observed_pre_cover,
            learner, predicted_occurrence, predicted_positive_loss,
            predicted_mortality
        ),
    shape_brms,
    shape_brt,
    shape_ensemble
) |>
    mutate(
        residual = observed_mortality - predicted_mortality,
        absolute_error = abs(residual),
        severe = observed_mortality >= 0.2,
        high_dhw = ann_maxdhw >= 10,
        high_dhw_high_acropora = high_dhw & prop_acropora_pre >= 0.7,
        low_mortality_event = event_year %in% c(2020L, 2022L)
    )

recommended_routes <- tibble(
    programme_key = c('ltmp', 'manta', 'mmp'),
    within_support_model = c(
        'shape_ensemble', 'benchmark_ensemble', 'shape_brms'
    ),
    outside_support_model = c(
        'shape_brms', 'benchmark_brms', 'shape_brms'
    ),
    reason = c(
        paste(
            'Shape ensemble improves overall and severe-tail validation;',
            'use BRMS alone beyond tree training support.'
        ),
        paste(
            'Shape candidate worsens event-held-out calibration;',
            'retain the benchmark and drop its BRT outside support.'
        ),
        paste(
            'Shape BRMS has the best overall, severe-tail, and 2024',
            'event-held-out performance.'
        )
    )
)
recommended_rows <- bind_rows(lapply(seq_len(nrow(recommended_routes)), function(i) {
    route <- recommended_routes[i, ]
    comparison_rows |>
        filter(
            programme_key == route$programme_key,
            learner == route$within_support_model
        ) |>
        mutate(learner = 'recommended_route')
}))
comparison_rows <- bind_rows(comparison_rows, recommended_rows)

overall_metrics <- comparison_rows |>
    group_by(programme_key, learner) |>
    summarise(
        n = n(),
        rmse = sqrt(mean(residual^2)),
        mae = mean(absolute_error),
        predictive_r2 = 1 - sum(residual^2) /
            sum((observed_mortality - mean(observed_mortality))^2),
        bias_observed_minus_predicted = mean(residual),
        severe_mae = mean(absolute_error[severe]),
        high_dhw_mae = mean(absolute_error[high_dhw]),
        high_dhw_bias = mean(residual[high_dhw]),
        high_dhw_high_acropora_mae = mean(
            absolute_error[high_dhw_high_acropora]
        ),
        high_dhw_high_acropora_bias = mean(
            residual[high_dhw_high_acropora]
        ),
        low_event_mae = mean(absolute_error[low_mortality_event]),
        low_event_bias = mean(residual[low_mortality_event]),
        .groups = 'drop'
    )

event_metrics <- comparison_rows |>
    group_by(programme_key, learner, event_year) |>
    summarise(
        n = n(), observed_mean = mean(observed_mortality),
        predicted_mean = mean(predicted_mortality),
        bias_observed_minus_predicted = mean(residual),
        rmse = sqrt(mean(residual^2)),
        severe_mae = mean(absolute_error[severe]),
        .groups = 'drop'
    )

brms_components <- function(cache, new_data, programme_key) {
    draws <- posterior_epred(cache$fit, newdata = new_data, re_formula = NA)
    if (programme_key == 'mmp') {
        boundary <- posterior_linpred(
            cache$fit, newdata = new_data, dpar = 'zoi',
            transform = TRUE, re_formula = NA
        )
        complete <- posterior_linpred(
            cache$fit, newdata = new_data, dpar = 'coi',
            transform = TRUE, re_formula = NA
        )
        occurrence <- 1 - boundary * (1 - complete)
    } else {
        zero_probability <- posterior_linpred(
            cache$fit, newdata = new_data, dpar = 'zi',
            transform = TRUE, re_formula = NA
        )
        occurrence <- 1 - zero_probability
    }
    positive <- draws / pmax(occurrence, 1e-8)
    tibble(
        occurrence = colMeans(occurrence),
        positive = pmin(pmax(colMeans(positive), 0), 1),
        mortality = colMeans(draws),
        q05 = apply(draws, 2, quantile, 0.05),
        q95 = apply(draws, 2, quantile, 0.95)
    )
}

brt_components <- function(cache, new_data) {
    trees <- as.integer(cache$tuning$best$n_trees[[1]])
    occurrence <- if (cache$model$occurrence$type == 'constant') {
        rep(cache$model$occurrence$value, nrow(new_data))
    } else {
        predict(
            cache$model$occurrence$model, newdata = new_data,
            n.trees = trees, type = 'response'
        )
    }
    positive <- if (cache$model$magnitude$type == 'constant') {
        rep(cache$model$magnitude$value, nrow(new_data))
    } else {
        predict(
            cache$model$magnitude$model, newdata = new_data,
            n.trees = trees, type = 'response'
        )
    }
    if (cache$model$magnitude_mode == 'logit_gaussian') positive <- plogis(positive)
    occurrence <- pmin(pmax(occurrence, 0), 1)
    positive <- pmin(pmax(positive, 0), 1)
    tibble(
        occurrence, positive, mortality = occurrence * positive,
        q05 = NA_real_, q95 = NA_real_
    )
}

dose_curves <- tibble()
baseline_weights <- read_csv(
    'output/formal_models/production_ensemble_weights.csv',
    show_col_types = FALSE
) |>
    filter(learner == 'formal_ensemble')

for (programme in names(validation_files)) {
    data <- load_programme_rows(programme)
    scenarios <- tibble(
        composition = c('Low Acropora', 'Typical', 'Acropora dominated'),
        prop_acropora_pre = c(0.1, median(data$prop_acropora_pre), 0.9),
        observed_pre_cover = c(
            median(data$observed_pre_cover),
            median(data$observed_pre_cover), 0.55
        )
    )
    grid_key <- crossing(
        scenarios,
        ann_maxdhw = seq(0, 14, by = 0.25)
    )

    caches <- list(
        benchmark_brms = readRDS(file.path(
            'output/models/formal', paste0('brms_', programme, '_production.rds')
        )),
        benchmark_brt = readRDS(file.path(
            'output/models/formal', paste0('brt_', programme, '_production.rds')
        )),
        shape_brms = readRDS(file.path(
            'output/models/shape_aware',
            paste0('shape_brms_', programme, '_production.rds')
        )),
        shape_brt = readRDS(file.path(
            'output/models/shape_aware',
            paste0('shape_brt_', programme, '_production.rds')
        ))
    )

    model_predictions <- list()
    for (model_name in names(caches)) {
        cache <- caches[[model_name]]
        new_data <- as_tibble(as.list(setNames(
            cache$preprocessing$median, cache$preprocessing$predictor
        )))
        new_data <- new_data[rep(1, nrow(grid_key)), , drop = FALSE]
        new_data$ann_maxdhw <- grid_key$ann_maxdhw
        new_data$prop_acropora_pre <- grid_key$prop_acropora_pre
        if ('observed_pre_cover' %in% names(new_data)) {
            new_data$observed_pre_cover <- grid_key$observed_pre_cover
        }
        if ('acropora_cover_pre' %in% names(new_data)) {
            new_data$acropora_cover_pre <-
                grid_key$prop_acropora_pre * grid_key$observed_pre_cover
        }
        new_data <- add_shape_features(new_data)
        new_data$reef_effect <- factor(
            levels(data$reef_effect)[1], levels = levels(data$reef_effect)
        )
        new_data$event_effect <- factor(
            levels(data$event_effect)[1], levels = levels(data$event_effect)
        )
        new_data$region_effect <- factor(
            levels(data$region_effect)[1], levels = levels(data$region_effect)
        )
        for (i in seq_len(nrow(cache$preprocessing))) {
            predictor <- cache$preprocessing$predictor[[i]]
            new_data[[paste0(predictor, '_z')]] <- (
                new_data[[predictor]] - cache$preprocessing$mean[[i]]
            ) / cache$preprocessing$sd[[i]]
        }
        model_predictions[[model_name]] <- if (grepl('brms$', model_name)) {
            brms_components(cache, new_data, programme)
        } else {
            brt_components(cache, new_data)
        }
        dose_curves <- bind_rows(
            dose_curves,
            bind_cols(
                tibble(programme_key = programme, learner = model_name),
                grid_key, model_predictions[[model_name]]
            )
        )
    }

    benchmark_weight <- baseline_weights |>
        filter(programme_key == programme) |>
        pull(weight_brms)
    shape_weight <- shape_production_weights |>
        filter(programme_key == programme) |>
        pull(weight_brms)
    for (ensemble_name in c('benchmark_ensemble', 'shape_ensemble')) {
        prefix <- sub('_ensemble$', '', ensemble_name)
        weight <- if (prefix == 'benchmark') benchmark_weight else shape_weight
        first <- model_predictions[[paste0(prefix, '_brms')]]
        second <- model_predictions[[paste0(prefix, '_brt')]]
        dose_curves <- bind_rows(
            dose_curves,
            bind_cols(
                tibble(programme_key = programme, learner = ensemble_name),
                grid_key,
                tibble(
                    occurrence = weight * first$occurrence +
                        (1 - weight) * second$occurrence,
                    positive = weight * first$positive +
                        (1 - weight) * second$positive,
                    mortality = weight * first$mortality +
                        (1 - weight) * second$mortality,
                    q05 = NA_real_, q95 = NA_real_
                )
            )
        )
    }

    route <- recommended_routes |>
        filter(programme_key == programme)
    within <- model_predictions[[route$within_support_model]]
    if (is.null(within)) {
        prefix <- sub('_ensemble$', '', route$within_support_model)
        weight <- if (prefix == 'benchmark') benchmark_weight else shape_weight
        first <- model_predictions[[paste0(prefix, '_brms')]]
        second <- model_predictions[[paste0(prefix, '_brt')]]
        within <- tibble(
            occurrence = weight * first$occurrence +
                (1 - weight) * second$occurrence,
            positive = weight * first$positive +
                (1 - weight) * second$positive,
            mortality = weight * first$mortality +
                (1 - weight) * second$mortality
        )
    }
    outside <- model_predictions[[route$outside_support_model]]
    observed_dhw_max <- max(data$ann_maxdhw, na.rm = TRUE)
    outside_support <- grid_key$ann_maxdhw > observed_dhw_max
    dose_curves <- bind_rows(
        dose_curves,
        bind_cols(
            tibble(programme_key = programme, learner = 'recommended_route'),
            grid_key,
            tibble(
                occurrence = if_else(
                    outside_support, outside$occurrence, within$occurrence
                ),
                positive = if_else(
                    outside_support, outside$positive, within$positive
                ),
                mortality = if_else(
                    outside_support, outside$mortality, within$mortality
                ),
                q05 = NA_real_, q95 = NA_real_,
                outside_tree_support = outside_support,
                observed_dhw_max = observed_dhw_max
            )
        )
    )
}

write_csv(comparison_rows, file.path(output_dir, 'comparison_rows.csv'))
write_csv(overall_metrics, file.path(output_dir, 'comparison_metrics.csv'))
write_csv(event_metrics, file.path(output_dir, 'event_metrics.csv'))
write_csv(shape_fold_weights, file.path(output_dir, 'shape_fold_weights.csv'))
write_csv(
    shape_production_weights,
    file.path(output_dir, 'shape_production_weights.csv')
)
write_csv(
    recommended_routes,
    file.path(output_dir, 'recommended_model_routes.csv')
)
write_csv(dose_curves, file.path(output_dir, 'production_dose_curves.csv'))

print(overall_metrics)
print(shape_production_weights)
