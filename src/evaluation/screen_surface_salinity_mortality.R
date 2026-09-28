# Screen two-dimensional eReefs surface-salinity summaries against the
# bleaching-mortality outcome. Only the season-ending 2024 layer currently
# overlaps a modelled bleaching event, so reef-blocked validation is the
# defensible predictive test; leave-one-event-out promotion is not possible.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(ggplot2)
    library(INLA)
    library(readr)
    library(splines)
    library(tidyr)
})

source('src/lib/surface_salinity_cause_helpers.R')
validate_surface_salinity_cause_classifier()

set.seed(20260915L)

salinity_file <- 'data/processed/ereefs_surface_salinity_reef_event.csv'
model_rows_file <- 'output/explanatory_event_dhw/event_dhw_brt_data.csv'
output_dir <- 'output/surface_salinity_mortality'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
effects_only <- '--effects-only' %in% commandArgs(trailingOnly = TRUE)

if (!all(file.exists(c(salinity_file, model_rows_file)))) {
    stop(
        'Missing inputs. Run src/data/build_reef_surface_salinity_metrics.R ',
        'and the explanatory-event-DHW build first.'
    )
}

salinity <- read_csv(salinity_file, show_col_types = FALSE)
if (anyDuplicated(salinity[c('ReefID', 'event_year')])) {
    stop('Salinity data have duplicate ReefID-event_year keys')
}

rows_all <- read_csv(model_rows_file, show_col_types = FALSE) |>
    mutate(
        ReefID = toupper(trimws(ReefID)),
        mortality_row_id = as.character(.mortality_row_id),
        is_manta = as.numeric(programme_key == 'manta'),
        is_mmp = as.numeric(programme_key == 'mmp'),
        depth_within_programme = as.numeric(depth_within_programme),
        disturbance_has_bleaching = coalesce(
            disturbance_has_bleaching, FALSE
        ),
        disturbance_has_cots = coalesce(disturbance_has_cots, FALSE),
        disturbance_has_cyclone = coalesce(
            disturbance_has_cyclone, FALSE
        ),
        disturbance_has_flood = coalesce(disturbance_has_flood, FALSE)
    ) |>
    left_join(
        salinity,
        by = c('ReefID', 'event_year'),
        relationship = 'many-to-one'
    ) |>
    classify_surface_salinity_causes()

if (anyDuplicated(rows_all$mortality_row_id)) {
    stop('Mortality row IDs are not unique after joining salinity')
}

overlap <- rows_all |>
    group_by(event_year) |>
    summarise(
        mortality_rows = n(),
        mortality_reefs = n_distinct(ReefID),
        salinity_rows = sum(is.finite(sss_min_area_mean_psu)),
        salinity_reefs = n_distinct(
            ReefID[is.finite(sss_min_area_mean_psu)]
        ),
        .groups = 'drop'
    ) |>
    full_join(
        salinity |>
            group_by(event_year, season) |>
            summarise(
                salinity_layer_reefs = sum(
                    is.finite(sss_min_area_mean_psu)
                ),
                .groups = 'drop'
            ),
        by = 'event_year'
    ) |>
    arrange(event_year)
write_csv(overlap, file.path(output_dir, 'event_overlap.csv'), na = '')

# Direct salinity is never fold-imputed. Restrict every candidate to the same
# observed rows so metric deltas cannot arise from changing the test sample.
direct_metrics <- c(
    'sss_min_area_mean_psu', 'sss_min_worst_part_psu',
    'log1p_hours_below30_area_mean',
    'log1p_hours_below30_worst_part'
)
rows_2024 <- rows_all |>
    filter(
        event_year == 2024L,
        if_all(all_of(direct_metrics), is.finite)
    )
if (nrow(rows_2024) == 0L) stop('No 2024 mortality rows match salinity')

cohort_rows <- list(
    freshwater_inclusive = rows_2024 |>
        filter(!disturbance_has_cots),
    bleaching_compatible = rows_2024 |>
        filter(!disturbance_has_cots, !disturbance_has_cyclone),
    target_descriptions_strict = rows_2024 |>
        filter(target_description_strict_reef),
    target_descriptions_flood_inclusive = rows_2024 |>
        filter(target_description_flood_reef)
)

write_csv(
    rows_2024 |>
        select(
            mortality_row_id, ReefID, ReefName, event_year, programme_key,
            depth, mortality_prop, joint_reef_fold, DISTURBANCE_TYPE,
            description, tooltip, disturbance_text,
            disturbance_has_bleaching, disturbance_has_cyclone,
            disturbance_has_flood, disturbance_has_cots,
            disturbance_has_disease, target_description_class,
            target_description_strict_row,
            target_description_flood_row,
            target_description_strict_reef,
            target_description_flood_reef,
            reef_has_flood_description,
            reef_has_competing_description,
            sss_min_area_mean_psu, sss_min_worst_part_psu,
            hours_below30_area_mean, hours_below30_worst_part
        ),
    file.path(output_dir, 'cause_classification_audit.csv'), na = ''
)
write_csv(
    cohort_rows$target_descriptions_flood_inclusive |>
        mutate(
            included_in_strict_target_set =
                target_description_strict_reef
        ) |>
        select(
            mortality_row_id, ReefID, ReefName, event_year, programme_key,
            depth, mortality_prop, joint_reef_fold,
            target_description_class, included_in_strict_target_set,
            reef_has_flood_description, DISTURBANCE_TYPE, description,
            tooltip, disturbance_text, everything()
        ),
    file.path(output_dir, 'cause_focused_reef_dataset.csv'), na = ''
)

cause_support <- rows_2024 |>
    group_by(programme_key, depth, target_description_class) |>
    summarise(
        rows = n(), reefs = n_distinct(ReefID),
        mean_mortality = mean(mortality_prop),
        median_mortality = median(mortality_prop),
        minimum_sss_area_mean_psu = min(sss_min_area_mean_psu),
        reefs_area_mean_below30 = n_distinct(
            ReefID[sss_min_area_mean_psu < 30]
        ),
        maximum_hours_below30_area_mean = max(
            hours_below30_area_mean
        ),
        strict_target_reefs = n_distinct(
            ReefID[target_description_strict_reef]
        ),
        flood_inclusive_target_reefs = n_distinct(
            ReefID[target_description_flood_reef]
        ),
        .groups = 'drop'
    )
write_csv(
    cause_support,
    file.path(output_dir, 'cause_description_support.csv'), na = ''
)

mmp_paired_depth <- cohort_rows$target_descriptions_flood_inclusive |>
    filter(programme_key == 'mmp', depth %in% c(2, 5)) |>
    group_by(
        ReefID, ReefName, depth, sss_min_area_mean_psu,
        hours_below30_area_mean, reef_has_flood_description
    ) |>
    summarise(
        rows = n(),
        mortality = mean(mortality_prop),
        description_classes = paste(
            sort(unique(target_description_class)), collapse = '; '
        ),
        .groups = 'drop'
    ) |>
    pivot_wider(
        names_from = depth,
        values_from = c(rows, mortality, description_classes),
        names_glue = '{.value}_{depth}m'
    ) |>
    filter(is.finite(mortality_2m), is.finite(mortality_5m)) |>
    mutate(shallow_minus_deep_mortality = mortality_2m - mortality_5m)
write_csv(
    mmp_paired_depth,
    file.path(output_dir, 'mmp_cause_focused_paired_depth.csv'), na = ''
)

bootstrap_paired_association <- function(x, y, metric, outcome,
                                         replicates = 2000L) {
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(x) < 5L || length(unique(x)) < 2L) {
        return(tibble(
            metric, outcome, paired_reefs = length(x), spearman_rho = NA_real_,
            rho_q025 = NA_real_, rho_q975 = NA_real_
        ))
    }
    set.seed(20261315L + match(
        paste(metric, outcome),
        c(
            'sss_min mortality_2m', 'sss_min mortality_5m',
            'sss_min shallow_minus_deep',
            'hours_below30 mortality_2m', 'hours_below30 mortality_5m',
            'hours_below30 shallow_minus_deep'
        )
    ))
    boot <- replicate(replicates, {
        index <- sample(seq_along(x), length(x), replace = TRUE)
        suppressWarnings(cor(x[index], y[index], method = 'spearman'))
    })
    boot <- boot[is.finite(boot)]
    tibble(
        metric, outcome, paired_reefs = length(x),
        spearman_rho = suppressWarnings(cor(x, y, method = 'spearman')),
        rho_q025 = if (length(boot)) quantile(boot, 0.025) else NA_real_,
        rho_q975 = if (length(boot)) quantile(boot, 0.975) else NA_real_
    )
}

mmp_paired_associations <- bind_rows(
    bootstrap_paired_association(
        mmp_paired_depth$sss_min_area_mean_psu,
        mmp_paired_depth$mortality_2m,
        'sss_min', 'mortality_2m'
    ),
    bootstrap_paired_association(
        mmp_paired_depth$sss_min_area_mean_psu,
        mmp_paired_depth$mortality_5m,
        'sss_min', 'mortality_5m'
    ),
    bootstrap_paired_association(
        mmp_paired_depth$sss_min_area_mean_psu,
        mmp_paired_depth$shallow_minus_deep_mortality,
        'sss_min', 'shallow_minus_deep'
    ),
    bootstrap_paired_association(
        mmp_paired_depth$hours_below30_area_mean,
        mmp_paired_depth$mortality_2m,
        'hours_below30', 'mortality_2m'
    ),
    bootstrap_paired_association(
        mmp_paired_depth$hours_below30_area_mean,
        mmp_paired_depth$mortality_5m,
        'hours_below30', 'mortality_5m'
    ),
    bootstrap_paired_association(
        mmp_paired_depth$hours_below30_area_mean,
        mmp_paired_depth$shallow_minus_deep_mortality,
        'hours_below30', 'shallow_minus_deep'
    )
)
write_csv(
    mmp_paired_associations,
    file.path(output_dir, 'mmp_cause_focused_paired_associations.csv'),
    na = ''
)

core_predictors <- c(
    'ann_maxdhw', 'observed_pre_cover', 'prop_acropora_pre',
    'log1p_cyc_interval_maxHrs4mw',
    'depth_within_programme', 'is_manta', 'is_mmp'
)
proxy_predictors <- c(
    core_predictors,
    'log_coastal_rain30', 'wqc_freqcc12',
    'wqc_prior10_percentile', 'wqc_10yr_sum'
)
candidate_predictors <- list(
    context_only = core_predictors,
    proxy_context = proxy_predictors,
    proxy_plus_sss_mean = c(proxy_predictors, 'sss_min_area_mean_psu'),
    proxy_plus_hours_mean = c(
        proxy_predictors, 'log1p_hours_below30_area_mean'
    ),
    proxy_plus_both_mean = c(
        proxy_predictors, 'sss_min_area_mean_psu',
        'log1p_hours_below30_area_mean'
    ),
    proxy_plus_both_tail = c(
        proxy_predictors, 'sss_min_worst_part_psu',
        'log1p_hours_below30_worst_part'
    )
)

prepare_fold <- function(analysis, assessment, predictors, scale = FALSE) {
    rules <- tibble(
        predictor = predictors,
        median = NA_real_, mean = NA_real_, sd = NA_real_,
        zero_variance = FALSE
    )
    active <- character()
    for (i in seq_along(predictors)) {
        predictor <- predictors[[i]]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        if (length(observed) == 0L) next
        fill <- median(observed)
        analysis[[predictor]][!is.finite(analysis[[predictor]])] <- fill
        assessment[[predictor]][!is.finite(assessment[[predictor]])] <- fill
        centre <- mean(analysis[[predictor]])
        spread <- sd(analysis[[predictor]])
        zero <- !is.finite(spread) || spread == 0
        rules$median[[i]] <- fill
        rules$mean[[i]] <- centre
        rules$sd[[i]] <- ifelse(zero, 1, spread)
        rules$zero_variance[[i]] <- zero
        if (!zero) {
            active <- c(active, predictor)
            if (scale) {
                z_name <- paste0(predictor, '_z')
                analysis[[z_name]] <-
                    (analysis[[predictor]] - centre) / spread
                assessment[[z_name]] <-
                    (assessment[[predictor]] - centre) / spread
            }
        }
    }
    list(
        analysis = analysis, assessment = assessment,
        active = active, rules = rules
    )
}

fit_gbm_component <- function(data, response, predictors, distribution,
                              seed, n_trees = 1200L) {
    outcome <- data[[response]]
    if (length(unique(outcome)) < 2L || !is.finite(sd(outcome)) ||
        sd(outcome) == 0) {
        return(list(type = 'constant', value = mean(outcome)))
    }
    set.seed(seed)
    model <- gbm(
        reformulate(predictors, response = response),
        data = data,
        distribution = distribution,
        n.trees = n_trees,
        interaction.depth = 2L,
        shrinkage = 0.02,
        n.minobsinnode = min(5L, max(2L, floor(nrow(data) / 20L))),
        bag.fraction = 0.7,
        train.fraction = 1,
        keep.data = FALSE,
        verbose = FALSE
    )
    list(type = 'gbm', model = model, n_trees = n_trees)
}

predict_gbm_component <- function(component, new_data, response_scale = TRUE) {
    if (component$type == 'constant') {
        return(rep(component$value, nrow(new_data)))
    }
    predict(
        component$model, new_data, n.trees = component$n_trees,
        type = ifelse(response_scale, 'response', 'link')
    )
}

fit_brt <- function(data, predictors, seed, n_trees = 1200L) {
    data$has_loss <- as.numeric(data$mortality_prop > 0)
    positive <- data |> filter(has_loss == 1)
    if (nrow(positive) == 0L) stop('No positive mortality rows for BRT')
    epsilon <- min(0.01, 0.5 / nrow(positive))
    positive$positive_logit <- qlogis(pmin(
        pmax(positive$mortality_prop, epsilon), 1 - epsilon
    ))
    list(
        occurrence = fit_gbm_component(
            data, 'has_loss', predictors, 'bernoulli', seed, n_trees
        ),
        magnitude = fit_gbm_component(
            positive, 'positive_logit', predictors, 'gaussian', seed + 1L,
            n_trees
        ),
        direct = fit_gbm_component(
            data, 'mortality_prop', predictors, 'gaussian', seed + 2L,
            n_trees
        )
    )
}

predict_brt <- function(model, assessment) {
    occurrence <- pmin(pmax(
        predict_gbm_component(model$occurrence, assessment), 0
    ), 1)
    magnitude_link <- predict_gbm_component(
        model$magnitude, assessment, response_scale = FALSE
    )
    magnitude <- if (model$magnitude$type == 'constant') {
        plogis(magnitude_link)
    } else {
        plogis(magnitude_link)
    }
    direct <- pmin(pmax(
        predict_gbm_component(model$direct, assessment), 0
    ), 1)
    tibble(
        predicted_occurrence = occurrence,
        predicted_positive_mortality = magnitude,
        predicted_two_part = occurrence * magnitude,
        predicted_direct = direct
    )
}

inla_random_term <- paste0(
    "f(reef_index, model='iid', ",
    "hyper=list(prec=list(prior='pc.prec', param=c(0.7,0.05))))"
)

fit_inla_two_part <- function(analysis, assessment, predictors,
                              compute_criteria = FALSE,
                              config = FALSE) {
    prepared <- prepare_fold(analysis, assessment, predictors, scale = TRUE)
    analysis <- prepared$analysis
    assessment <- prepared$assessment
    terms <- paste0(prepared$active, '_z')
    if (length(terms) == 0L) stop('No active predictors for INLA')

    all_levels <- bind_rows(analysis, assessment)
    reef_levels <- unique(all_levels$ReefID)
    index_rows <- function(x) {
        x |> mutate(reef_index = match(ReefID, reef_levels))
    }
    analysis <- index_rows(analysis)
    assessment <- index_rows(assessment)
    rhs <- paste(c(terms, inla_random_term), collapse = ' + ')

    occurrence_data <- bind_rows(
        analysis |> mutate(response = as.numeric(mortality_prop > 0)),
        assessment |> mutate(response = NA_real_)
    )
    occurrence <- inla(
        as.formula(paste('response ~', rhs)),
        family = 'binomial', Ntrials = 1, data = occurrence_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0,
            prec.intercept = 1
        ),
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(
            waic = compute_criteria, dic = compute_criteria,
            config = config
        ),
        verbose = FALSE
    )

    positive <- analysis |> filter(mortality_prop > 0)
    positive_n <- nrow(positive)
    if (positive_n < 5L) stop('Too few positive rows for INLA beta model')
    positive$response <- (
        positive$mortality_prop * (positive_n - 1) + 0.5
    ) / positive_n
    magnitude_data <- bind_rows(
        positive,
        assessment |> mutate(response = NA_real_)
    )
    magnitude <- inla(
        as.formula(paste('response ~', rhs)),
        family = 'beta', data = magnitude_data,
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
        control.compute = list(
            waic = compute_criteria, dic = compute_criteria,
            config = config
        ),
        verbose = FALSE
    )

    occurrence_index <- nrow(analysis) + seq_len(nrow(assessment))
    magnitude_index <- positive_n + seq_len(nrow(assessment))
    occurrence_summary <-
        occurrence$summary.fitted.values[occurrence_index, ]
    magnitude_summary <-
        magnitude$summary.fitted.values[magnitude_index, ]
    list(
        predictions = tibble(
            predicted_occurrence = occurrence_summary$mean,
            predicted_positive_mortality = magnitude_summary$mean,
            predicted_mortality = occurrence_summary$mean *
                magnitude_summary$mean
        ),
        occurrence = occurrence,
        magnitude = magnitude,
        preprocessing = prepared$rules,
        active_predictors = prepared$active
    )
}

metric_summary <- function(predictions) {
    predictions |>
        group_by(cohort, candidate, framework) |>
        summarise(
            rows = n(), reefs = n_distinct(ReefID),
            rmse = sqrt(mean(
                (observed_mortality - predicted_mortality)^2
            )),
            mae = mean(abs(observed_mortality - predicted_mortality)),
            predictive_r2 = 1 - sum(
                (observed_mortality - predicted_mortality)^2
            ) / sum(
                (observed_mortality - mean(observed_mortality))^2
            ),
            bias = mean(predicted_mortality - observed_mortality),
            occurrence_brier = if (all(is.finite(predicted_occurrence))) {
                mean((observed_occurrence - predicted_occurrence)^2)
            } else NA_real_,
            severe_rows = sum(observed_mortality >= 0.2),
            severe_rmse = if (severe_rows > 0L) sqrt(mean(
                (observed_mortality[observed_mortality >= 0.2] -
                    predicted_mortality[observed_mortality >= 0.2])^2
            )) else NA_real_,
            false_extreme_rate = mean(
                predicted_mortality >= 0.5 & observed_mortality < 0.2
            ),
            .groups = 'drop'
        )
}

if (!effects_only) {
predictions <- tibble()
preprocessing <- tibble()
for (cohort_name in names(cohort_rows)) {
    cohort <- cohort_rows[[cohort_name]]
    folds <- sort(unique(cohort$joint_reef_fold))
    for (candidate in names(candidate_predictors)) {
        predictors <- candidate_predictors[[candidate]]
        for (fold in folds) {
            held_out <- cohort$joint_reef_fold == fold
            analysis <- cohort[!held_out, , drop = FALSE]
            assessment <- cohort[held_out, , drop = FALSE]

            prepared_brt <- prepare_fold(
                analysis, assessment, predictors, scale = FALSE
            )
            brt <- fit_brt(
                prepared_brt$analysis, prepared_brt$active,
                20260915L + 1000L * match(cohort_name, names(cohort_rows)) +
                    100L * match(candidate, names(candidate_predictors)) +
                    fold
            )
            brt_prediction <- predict_brt(brt, prepared_brt$assessment)
            prediction_keys <- prepared_brt$assessment |>
                transmute(
                    mortality_row_id, ReefID, ReefName, programme_key,
                    event_year, depth, joint_reef_fold,
                    observed_mortality = mortality_prop,
                    observed_occurrence = as.numeric(mortality_prop > 0)
                )
            predictions <- bind_rows(
                predictions,
                bind_cols(prediction_keys, brt_prediction) |>
                    transmute(
                        across(everything()), cohort = cohort_name,
                        candidate, framework = 'BRT two-part',
                        predicted_mortality = predicted_two_part
                    ) |>
                    select(-predicted_two_part, -predicted_direct),
                bind_cols(prediction_keys, brt_prediction) |>
                    transmute(
                        across(everything()), cohort = cohort_name,
                        candidate, framework = 'BRT direct',
                        predicted_mortality = predicted_direct,
                        predicted_occurrence = NA_real_,
                        predicted_positive_mortality = NA_real_
                    ) |>
                    select(-predicted_two_part, -predicted_direct)
            )
            preprocessing <- bind_rows(
                preprocessing,
                prepared_brt$rules |>
                    mutate(
                        cohort = cohort_name, candidate,
                        framework = 'BRT', fold
                    )
            )

            fitted_inla <- fit_inla_two_part(
                analysis, assessment, predictors
            )
            predictions <- bind_rows(
                predictions,
                bind_cols(prediction_keys, fitted_inla$predictions) |>
                    mutate(
                        cohort = cohort_name, candidate,
                        framework = 'INLA two-part'
                    )
            )
            preprocessing <- bind_rows(
                preprocessing,
                fitted_inla$preprocessing |>
                    mutate(
                        cohort = cohort_name, candidate,
                        framework = 'INLA', fold
                    )
            )
            message(
                cohort_name, ' / ', candidate, ' / reef fold ', fold
            )
        }
    }
}

metrics <- metric_summary(predictions)
metrics_by_programme_depth <- predictions |>
    group_by(
        cohort, candidate, framework, programme_key, depth
    ) |>
    summarise(
        rows = n(), reefs = n_distinct(ReefID),
        rmse = sqrt(mean(
            (observed_mortality - predicted_mortality)^2
        )),
        mae = mean(abs(observed_mortality - predicted_mortality)),
        bias = mean(predicted_mortality - observed_mortality),
        severe_rows = sum(observed_mortality >= 0.2),
        severe_rmse = if (severe_rows > 0L) sqrt(mean(
            (observed_mortality[observed_mortality >= 0.2] -
                predicted_mortality[observed_mortality >= 0.2])^2
        )) else NA_real_,
        .groups = 'drop'
    )
write_csv(
    predictions,
    file.path(output_dir, 'reef_blocked_predictions.csv'), na = ''
)
write_csv(
    metrics,
    file.path(output_dir, 'reef_blocked_metrics.csv'), na = ''
)
write_csv(
    metrics_by_programme_depth,
    file.path(output_dir, 'reef_blocked_metrics_by_programme_depth.csv'),
    na = ''
)
write_csv(
    preprocessing,
    file.path(output_dir, 'fold_preprocessing.csv'), na = ''
)

bootstrap_deltas <- function(prediction_rows, replicates = 1000L) {
    baseline <- prediction_rows |>
        filter(candidate == 'proxy_context') |>
        select(
            cohort, framework, mortality_row_id, ReefID,
            observed_mortality, baseline_prediction = predicted_mortality
        )
    comparisons <- prediction_rows |>
        filter(candidate != 'proxy_context') |>
        inner_join(
            baseline,
            by = c(
                'cohort', 'framework', 'mortality_row_id', 'ReefID',
                'observed_mortality'
            ),
            relationship = 'many-to-one'
        )
    results <- tibble()
    groups <- comparisons |>
        distinct(cohort, framework, candidate)
    for (i in seq_len(nrow(groups))) {
        spec <- groups[i, ]
        x <- comparisons |>
            filter(
                cohort == spec$cohort,
                framework == spec$framework,
                candidate == spec$candidate
            )
        reef_ids <- unique(x$ReefID)
        set.seed(20261015L + i)
        delta <- replicate(replicates, {
            sampled <- sample(reef_ids, length(reef_ids), replace = TRUE)
            index <- unlist(lapply(sampled, function(id) which(x$ReefID == id)))
            sqrt(mean(
                (x$observed_mortality[index] -
                    x$predicted_mortality[index])^2
            )) - sqrt(mean(
                (x$observed_mortality[index] -
                    x$baseline_prediction[index])^2
            ))
        })
        observed_delta <- sqrt(mean(
            (x$observed_mortality - x$predicted_mortality)^2
        )) - sqrt(mean(
            (x$observed_mortality - x$baseline_prediction)^2
        ))
        results <- bind_rows(
            results,
            tibble(
                cohort = spec$cohort,
                framework = spec$framework,
                candidate = spec$candidate,
                rows = nrow(x), reefs = length(reef_ids),
                delta_rmse = observed_delta,
                delta_rmse_q025 = quantile(delta, 0.025),
                delta_rmse_q975 = quantile(delta, 0.975),
                bootstrap_probability_improved = mean(delta < 0)
            )
        )
    }
    results
}

metric_deltas <- bootstrap_deltas(predictions)
write_csv(
    metric_deltas,
    file.path(output_dir, 'reef_blocked_metric_deltas.csv'), na = ''
)

coverage <- bind_rows(lapply(names(cohort_rows), function(cohort_name) {
    cohort_rows[[cohort_name]] |>
        group_by(programme_key, depth) |>
        summarise(
            rows = n(), reefs = n_distinct(ReefID),
            positive_rows = sum(mortality_prop > 0),
            severe_rows = sum(mortality_prop >= 0.2),
            mean_mortality = mean(mortality_prop),
            reefs_area_mean_below30 = n_distinct(
                ReefID[sss_min_area_mean_psu < 30]
            ),
            reefs_worst_part_below30 = n_distinct(
                ReefID[sss_min_worst_part_psu < 30]
            ),
            minimum_sss_area_mean_psu = min(sss_min_area_mean_psu),
            maximum_hours_below30_area_mean = max(
                hours_below30_area_mean
            ),
            .groups = 'drop'
        ) |>
        mutate(cohort = cohort_name)
}))
write_csv(coverage, file.path(output_dir, 'analysis_coverage.csv'), na = '')

correlation_variables <- c(
    'mortality_prop', 'sss_min_area_mean_psu',
    'hours_below30_area_mean', 'ann_maxdhw',
    'log1p_cyc_interval_maxHrs4mw', 'log_coastal_rain30',
    'wqc_freqcc12', 'wqc_prior10_percentile', 'wqc_10yr_sum'
)
correlations <- bind_rows(lapply(names(cohort_rows), function(cohort_name) {
    x <- cohort_rows[[cohort_name]][correlation_variables]
    as.data.frame(cor(x, method = 'spearman', use = 'pairwise.complete.obs')) |>
        tibble::rownames_to_column('variable_1') |>
        pivot_longer(-variable_1, names_to = 'variable_2', values_to = 'rho') |>
        mutate(cohort = cohort_name)
}))
write_csv(correlations, file.path(output_dir, 'correlations.csv'), na = '')

# Full-data BRT importance remains descriptive. It is retained because the
# direct BRT is the project's nonlinear diagnostic, not a selected forecast.
importance <- tibble()
full_brt_models <- list()
for (cohort_name in names(cohort_rows)) {
    x <- cohort_rows[[cohort_name]]
    for (candidate in c('proxy_plus_both_mean', 'proxy_plus_both_tail')) {
        prepared <- prepare_fold(
            x, x[0, , drop = FALSE], candidate_predictors[[candidate]],
            scale = FALSE
        )
        model <- fit_brt(
            prepared$analysis, prepared$active,
            20261115L + match(cohort_name, names(cohort_rows)) * 10L +
                match(candidate, names(candidate_predictors))
        )
        full_brt_models[[paste(cohort_name, candidate, sep = '__')]] <- model
        for (component in names(model)) {
            if (model[[component]]$type != 'gbm') next
            importance <- bind_rows(
                importance,
                summary(model[[component]]$model, plotit = FALSE) |>
                    as_tibble() |>
                    transmute(
                        cohort = cohort_name, candidate, component,
                        predictor = var,
                        relative_influence = rel.inf,
                        rank = rank(-rel.inf, ties.method = 'min')
                    )
            )
        }
    }
}
write_csv(
    importance, file.path(output_dir, 'full_data_brt_importance.csv'), na = ''
)
saveRDS(
    full_brt_models, file.path(output_dir, 'full_data_brt_models.rds')
)
} else {
    predictions <- read_csv(
        file.path(output_dir, 'reef_blocked_predictions.csv'),
        show_col_types = FALSE
    )
    metrics <- read_csv(
        file.path(output_dir, 'reef_blocked_metrics.csv'),
        show_col_types = FALSE
    )
    metric_deltas <- read_csv(
        file.path(output_dir, 'reef_blocked_metric_deltas.csv'),
        show_col_types = FALSE
    )
}

# Programme-specific effect models. Only MMP estimates a within-programme
# contrast between 2 m and 5 m; Manta and LTMP are separate 9 m curves.
effect_base_predictors <- c(
    'ann_maxdhw', 'observed_pre_cover', 'prop_acropora_pre',
    'log1p_cyc_interval_maxHrs4mw', 'log_coastal_rain30',
    'wqc_freqcc12', 'wqc_prior10_percentile', 'wqc_10yr_sum'
)
effect_specs <- list(
    sss_min = list(
        raw = 'sss_min_area_mean_psu',
        model = 'sss_min_area_mean_psu',
        label = 'Area-mean minimum surface salinity (PSU)',
        low_is_stress = TRUE
    ),
    hours_below30 = list(
        raw = 'hours_below30_area_mean',
        model = 'log1p_hours_below30_area_mean',
        label = 'Area-mean hours below 30 PSU',
        low_is_stress = FALSE
    )
)

fit_effect_inla <- function(data, programme, spec, posterior_samples = 1000L) {
    data <- data |>
        filter(
            programme_key == .env$programme,
            is.finite(.data[[spec$raw]]),
            is.finite(.data[[spec$model]])
        )
    prepared <- prepare_fold(
        data, data[0, , drop = FALSE], effect_base_predictors,
        scale = TRUE
    )
    data <- prepared$analysis
    active_base <- paste0(prepared$active, '_z')
    model_x <- data[[spec$model]]
    boundary_knots <- range(model_x)
    interior_values <- model_x[
        model_x > boundary_knots[[1]] & model_x < boundary_knots[[2]]
    ]
    if (length(interior_values) == 0L) {
        stop('No interior exposure values for ', programme, ' / ', spec$raw)
    }
    interior_knot <- median(interior_values)
    basis <- ns(
        model_x, knots = interior_knot,
        Boundary.knots = boundary_knots
    )
    data$metric_ns1 <- basis[, 1]
    data$metric_ns2 <- basis[, 2]
    depth_terms <- character()
    if (programme == 'mmp') {
        data$depth5 <- as.numeric(data$depth == 5)
        data$metric_ns1_depth5 <- data$metric_ns1 * data$depth5
        data$metric_ns2_depth5 <- data$metric_ns2 * data$depth5
        depth_terms <- c(
            'depth5', 'metric_ns1_depth5', 'metric_ns2_depth5'
        )
    }
    fixed_terms <- c(active_base, 'metric_ns1', 'metric_ns2', depth_terms)
    data$reef_index <- match(data$ReefID, unique(data$ReefID))
    rhs <- paste(c(fixed_terms, inla_random_term), collapse = ' + ')

    occurrence_data <- data |>
        mutate(response = as.numeric(mortality_prop > 0))
    occurrence <- inla(
        as.formula(paste('response ~', rhs)),
        family = 'binomial', Ntrials = 1, data = occurrence_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0,
            prec.intercept = 1
        ),
        control.predictor = list(compute = TRUE, link = 1),
        control.compute = list(config = TRUE), verbose = FALSE
    )
    positive <- data |> filter(mortality_prop > 0)
    positive_n <- nrow(positive)
    positive$response <- (
        positive$mortality_prop * (positive_n - 1) + 0.5
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
        control.compute = list(config = TRUE), verbose = FALSE
    )

    raw_range <- quantile(data[[spec$raw]], c(0.02, 0.98), na.rm = TRUE)
    if (!spec$low_is_stress) raw_range[[1]] <- 0
    raw_grid <- seq(raw_range[[1]], raw_range[[2]], length.out = 80L)
    model_grid <- if (spec$model == spec$raw) raw_grid else log1p(raw_grid)
    grid_basis <- ns(
        model_grid,
        knots = attr(basis, 'knots'),
        Boundary.knots = attr(basis, 'Boundary.knots')
    )
    depths <- if (programme == 'mmp') c(2, 5) else 9
    grid <- tidyr::crossing(depth = depths, grid_index = seq_along(raw_grid)) |>
        mutate(
            metric_value = raw_grid[grid_index],
            metric_ns1 = grid_basis[grid_index, 1],
            metric_ns2 = grid_basis[grid_index, 2]
        )
    for (term in active_base) grid[[term]] <- 0
    if (programme == 'mmp') {
        grid$depth5 <- as.numeric(grid$depth == 5)
        grid$metric_ns1_depth5 <- grid$metric_ns1 * grid$depth5
        grid$metric_ns2_depth5 <- grid$metric_ns2 * grid$depth5
    }
    fixed_formula <- reformulate(fixed_terms)
    design <- model.matrix(fixed_formula, grid)

    fixed_draws <- function(fit, samples) {
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
    occurrence_beta <- fixed_draws(occurrence, posterior_samples)
    magnitude_beta <- fixed_draws(magnitude, posterior_samples)
    occurrence_probability <- plogis(design %*% occurrence_beta)
    positive_magnitude <- plogis(design %*% magnitude_beta)
    expected_mortality <- occurrence_probability * positive_magnitude

    reference_index <- vapply(depths, function(target_depth) {
        candidates <- which(grid$depth == target_depth)
        if (spec$low_is_stress) tail(candidates, 1) else candidates[[1]]
    }, integer(1))
    reference_for_row <- reference_index[match(grid$depth, depths)]
    delta <- expected_mortality - expected_mortality[reference_for_row, ]
    effect <- grid |>
        transmute(
            programme_key = programme, depth, metric_value,
            predicted_mortality = rowMeans(expected_mortality),
            lower = apply(expected_mortality, 1, quantile, 0.025),
            upper = apply(expected_mortality, 1, quantile, 0.975),
            delta_from_reference = rowMeans(delta),
            delta_lower = apply(delta, 1, quantile, 0.025),
            delta_upper = apply(delta, 1, quantile, 0.975),
            reference_value = raw_grid[
                ifelse(spec$low_is_stress, length(raw_grid), 1L)
            ],
            framework = 'INLA two-part'
        )
    list(
        effect = effect,
        occurrence = occurrence,
        magnitude = magnitude,
        preprocessing = prepared$rules,
        basis = list(
            knots = attr(basis, 'knots'),
            boundary_knots = attr(basis, 'Boundary.knots')
        )
    )
}

fit_effect_brt <- function(data, programme, spec, bootstraps = 60L,
                           seed = 1L) {
    data <- data |>
        filter(
            programme_key == .env$programme,
            is.finite(.data[[spec$raw]]),
            is.finite(.data[[spec$model]])
        )
    predictors <- c(effect_base_predictors, spec$model)
    if (programme == 'mmp') predictors <- c(predictors, 'depth')
    prepared <- prepare_fold(
        data, data[0, , drop = FALSE], predictors, scale = FALSE
    )
    data <- prepared$analysis
    active <- prepared$active
    raw_range <- quantile(data[[spec$raw]], c(0.02, 0.98), na.rm = TRUE)
    if (!spec$low_is_stress) raw_range[[1]] <- 0
    raw_grid <- seq(raw_range[[1]], raw_range[[2]], length.out = 80L)
    depths <- if (programme == 'mmp') c(2, 5) else 9

    partial_curve <- function(model, integration_data) {
        bind_rows(lapply(depths, function(target_depth) {
            rows_per_grid <- nrow(integration_data)
            new_data <- integration_data[
                rep(seq_len(rows_per_grid), times = length(raw_grid)),
                , drop = FALSE
            ]
            new_data[[spec$raw]] <- rep(raw_grid, each = rows_per_grid)
            if (spec$model != spec$raw) {
                new_data[[spec$model]] <- log1p(new_data[[spec$raw]])
            }
            if (programme == 'mmp') new_data$depth <- target_depth
            prediction <- predict_brt(model, new_data)
            tibble(
                depth = target_depth,
                metric_value = rep(raw_grid, each = rows_per_grid),
                predicted_mortality = prediction$predicted_direct
            ) |>
                group_by(depth, metric_value) |>
                summarise(
                    predicted_mortality = mean(predicted_mortality),
                    .groups = 'drop'
                )
        }))
    }

    model <- fit_brt(data, active, seed, n_trees = 500L)
    point <- partial_curve(model, data)
    reefs <- unique(data$ReefID)
    set.seed(seed + 10000L)
    bootstrap_curves <- tibble()
    for (b in seq_len(bootstraps)) {
        sampled <- sample(reefs, length(reefs), replace = TRUE)
        index <- unlist(lapply(sampled, function(id) which(data$ReefID == id)))
        boot_data <- data[index, , drop = FALSE]
        boot_prepared <- prepare_fold(
            boot_data, boot_data[0, , drop = FALSE], predictors,
            scale = FALSE
        )
        boot_model <- tryCatch(
            fit_brt(
                boot_prepared$analysis, boot_prepared$active, seed + b,
                n_trees = 500L
            ),
            error = function(e) NULL
        )
        if (is.null(boot_model)) next
        bootstrap_curves <- bind_rows(
            bootstrap_curves,
            partial_curve(boot_model, boot_prepared$analysis) |>
                mutate(bootstrap = b)
        )
    }
    if (n_distinct(bootstrap_curves$bootstrap) < 50L) {
        stop('Fewer than 50 successful BRT bootstrap fits for ', programme)
    }
    reference <- bootstrap_curves |>
        group_by(bootstrap, depth) |>
        slice(if (spec$low_is_stress) n() else 1L) |>
        ungroup() |>
        select(bootstrap, depth, reference_prediction = predicted_mortality)
    bootstrap_curves <- bootstrap_curves |>
        left_join(
            reference, by = c('bootstrap', 'depth'),
            relationship = 'many-to-one'
        ) |>
        mutate(delta = predicted_mortality - reference_prediction)
    interval <- bootstrap_curves |>
        group_by(depth, metric_value) |>
        summarise(
            lower = quantile(predicted_mortality, 0.025),
            upper = quantile(predicted_mortality, 0.975),
            delta_from_reference = mean(delta),
            delta_lower = quantile(delta, 0.025),
            delta_upper = quantile(delta, 0.975),
            .groups = 'drop'
        )
    point |>
        left_join(interval, by = c('depth', 'metric_value')) |>
        mutate(
            programme_key = programme,
            reference_value = if (spec$low_is_stress) {
                max(raw_grid)
            } else min(raw_grid),
            framework = 'BRT direct'
        )
}

effect_cohort_programmes <- list(
    freshwater_inclusive = c('mmp', 'manta', 'ltmp'),
    target_descriptions_flood_inclusive = 'mmp'
)
effect_curves <- tibble()
effect_models <- list()
for (effect_cohort_name in names(effect_cohort_programmes)) {
    primary_effect_rows <- cohort_rows[[effect_cohort_name]]
    programmes <- effect_cohort_programmes[[effect_cohort_name]]
    for (metric_name in names(effect_specs)) {
        spec <- effect_specs[[metric_name]]
        for (programme in programmes) {
            fitted_inla <- fit_effect_inla(
                primary_effect_rows, programme, spec
            )
            model_key <- paste(
                'inla', effect_cohort_name, metric_name, programme,
                sep = '__'
            )
            effect_models[[model_key]] <- fitted_inla[c(
                'occurrence', 'magnitude', 'preprocessing', 'basis'
            )]
            effect_curves <- bind_rows(
                effect_curves,
                fitted_inla$effect |>
                    mutate(
                        cohort = effect_cohort_name,
                        metric = metric_name, metric_label = spec$label
                    ),
                fit_effect_brt(
                    primary_effect_rows, programme, spec,
                    seed = 20261215L +
                        10000L * match(
                            effect_cohort_name,
                            names(effect_cohort_programmes)
                        ) +
                        100L * match(metric_name, names(effect_specs)) +
                        match(programme, programmes)
                ) |>
                    mutate(
                        cohort = effect_cohort_name,
                        metric = metric_name, metric_label = spec$label
                    )
            )
            message(
                'Effect curves: ', effect_cohort_name, ' / ',
                metric_name, ' / ', programme
            )
        }
    }
}
effect_curves <- effect_curves |>
    mutate(
        effect_cohort_label = recode(
            cohort,
            freshwater_inclusive = 'Broad no-COTS sensitivity',
            target_descriptions_flood_inclusive =
                'Target descriptions: flood/bleaching/none'
        ),
        programme_depth = case_when(
            programme_key == 'mmp' ~ paste0('MMP ', depth, ' m'),
            programme_key == 'manta' ~ 'Manta 9 m',
            TRUE ~ 'LTMP 9 m'
        )
    )
write_csv(
    effect_curves, file.path(output_dir, 'depth_effect_curves.csv'), na = ''
)
saveRDS(effect_models, file.path(output_dir, 'inla_depth_effect_models.rds'))

count_supporting_reefs <- function(cohort_value, programme_value,
                                   metric_value, threshold_value) {
    if (is.na(threshold_value)) return(NA_integer_)
    x <- cohort_rows[[cohort_value]] |>
        filter(programme_key == programme_value)
    if (metric_value == 'sss_min') {
        n_distinct(x$ReefID[x$sss_min_area_mean_psu <= threshold_value])
    } else {
        n_distinct(x$ReefID[x$hours_below30_area_mean >= threshold_value])
    }
}

thresholds <- effect_curves |>
    group_by(cohort, effect_cohort_label, framework, metric, metric_label,
             programme_key, depth, programme_depth, reference_value) |>
    summarise(
        threshold = {
            supported <- metric_value[delta_lower > 0]
            if (length(supported) == 0L) NA_real_ else if (
                first(metric) == 'sss_min'
            ) max(supported) else min(supported)
        },
        maximum_estimated_uplift = max(delta_from_reference),
        .groups = 'drop'
    ) |>
    rowwise() |>
    mutate(
        supporting_reefs = count_supporting_reefs(
            cohort, programme_key, metric, threshold
        ),
        threshold_supported_by_at_least_five_reefs =
            supporting_reefs >= 5L
    ) |>
    ungroup()
write_csv(thresholds, file.path(output_dir, 'effect_thresholds.csv'), na = '')

raw_plot_data <- bind_rows(lapply(
    names(effect_cohort_programmes),
    function(cohort_name) {
        cohort_rows[[cohort_name]] |>
            filter(programme_key %in% effect_cohort_programmes[[cohort_name]]) |>
            mutate(cohort = cohort_name)
    }
)) |>
    mutate(
        programme_depth = case_when(
            programme_key == 'mmp' ~ paste0('MMP ', depth, ' m'),
            programme_key == 'manta' ~ 'Manta 9 m',
            TRUE ~ 'LTMP 9 m'
        )
    )
write_csv(
    raw_plot_data |>
        select(
            cohort, mortality_row_id, ReefID, ReefName, programme_key, depth,
            programme_depth, mortality_prop, ann_maxdhw,
            target_description_class, reef_has_flood_description,
            sss_min_area_mean_psu, sss_min_worst_part_psu,
            hours_below30_area_mean, hours_below30_worst_part,
            disturbance_has_cyclone, disturbance_has_flood,
            cyc_interval_maxHrs4mw
        ),
    file.path(output_dir, 'effect_observations.csv'), na = ''
)

depth_plot <- ggplot(
    effect_curves,
    aes(metric_value, predicted_mortality, colour = programme_depth,
        fill = programme_depth)
) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.10,
                colour = NA) +
    geom_line(linewidth = 0.8) +
    facet_grid(effect_cohort_label + framework ~ metric_label,
               scales = 'free_x') +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
        x = NULL, y = 'Expected relative mortality',
        colour = 'Survey programme/depth', fill = 'Survey programme/depth',
        title = 'MMP target-description curves test the apparent salinity threshold',
        subtitle = paste(
            'Broad no-COTS curves retain all programmes;',
            'focused curves retain only flood, bleaching or blank annotations'
        )
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'bottom')
ggsave(
    file.path(output_dir, 'depth_effect_curves.png'), depth_plot,
    width = 12, height = 11, dpi = 180
)

delta_plot <- metric_deltas |>
    filter(candidate != 'context_only') |>
    ggplot(aes(delta_rmse, candidate, colour = framework)) +
    geom_vline(xintercept = 0, colour = 'grey50', linetype = 2) +
    geom_errorbar(
        aes(xmin = delta_rmse_q025, xmax = delta_rmse_q975),
        orientation = 'y', width = 0.2,
        position = position_dodge(width = 0.5)
    ) +
    geom_point(position = position_dodge(width = 0.5)) +
    facet_wrap(~ cohort, scales = 'free_y') +
    labs(
        x = 'Change in reef-blocked RMSE versus existing-proxy context',
        y = NULL, colour = 'Framework',
        title = 'Direct surface salinity must improve a fixed-row proxy baseline',
        subtitle = 'Negative values favour the candidate; intervals resample reefs'
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = 'bottom')
ggsave(
    file.path(output_dir, 'reef_blocked_rmse_deltas.png'), delta_plot,
    width = 11, height = 7, dpi = 180
)

cat('Wrote surface-salinity mortality screen to', output_dir, '\n')
print(metrics |> arrange(cohort, framework, rmse))
print(thresholds)
