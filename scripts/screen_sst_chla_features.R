# Blocked BRT and INLA screen for SST distribution shape and median chlorophyll.
# This is a feature screen, not the final model: BRT settings and INLA priors
# are held fixed so candidate predictor sets receive a fair, fast comparison.

suppressPackageStartupMessages({
    library(dplyr)
    library(gbm)
    library(INLA)
    library(readr)
    library(tidyr)
})
source('scripts/formal_model_helpers.R')

set.seed(20260826L)
output_dir <- 'output/sst_chla_screen'
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

new_features <- read_csv(
    'data/processed/sst_chla_features_validation.csv',
    show_col_types = FALSE
) |>
    select(
        ReefID, year,
        sst_summer_skewness, sst_summer_excess_kurtosis,
        sst_prevyear_skewness, sst_prevyear_excess_kurtosis,
        chla_wetseason_median, chla_wetseason_n,
        chla_q1_median, chla_q1_n
    )
if (anyDuplicated(new_features[c('ReefID', 'year')])) {
    stop('SST/chlorophyll feature table has duplicate reef-year keys')
}

load_screen_rows <- function(programme) {
    rows <- load_programme_rows(programme) |>
        left_join(
            new_features,
            by = c('ReefID', 'event_year' = 'year'),
            relationship = 'many-to-one'
        )
    if (any(!is.finite(rows$sst_summer_skewness))) {
        stop('Missing SST shape features after joining ', programme)
    }
    if (any(!is.finite(rows$chla_wetseason_median))) {
        stop('Missing wet-season chlorophyll after joining ', programme)
    }
    rows
}

# The earlier INLA history screen supported the newer exposure history for
# manta tow, but not LTMP or reef-blocked MMP. Hold that decision fixed here so
# this test isolates the new SST/chlorophyll information.
history_for_programme <- c(manta = 'paper', ltmp = 'legacy', mmp = 'legacy')
history_variables <- list(
    legacy = c('histmDHW6', 'yrsince6'),
    paper = c(
        'dhw_events_since2016_n6',
        'dhw_years_since_last_n6_capped8', 'dhw_no_prior_n6'
    )
)

anchors <- c('ann_maxdhw', 'prop_acropora_pre', 'observed_pre_cover')
existing_modifiers <- c(
    'dhw10_load4', 'dhw_novelty10', 'secc3m_p10', 'cloudp_90'
)
new_primary <- c(
    'sst_summer_skewness', 'sst_summer_excess_kurtosis',
    'chla_wetseason_median'
)
new_sensitivity <- c(
    'sst_prevyear_skewness', 'sst_prevyear_excess_kurtosis'
)
candidates <- c('core', 'sst_chla_raw', 'modifier_pca')
schemes <- c('reef_blocked_5fold', 'leave_one_event_out')

impute_from_analysis <- function(analysis, assessment, variables) {
    rules <- tibble(variable = variables, median = NA_real_)
    for (i in seq_along(variables)) {
        variable <- variables[[i]]
        observed <- analysis[[variable]][is.finite(analysis[[variable]])]
        if (length(observed) == 0) stop('No analysis values for ', variable)
        value <- median(observed)
        analysis[[variable]][!is.finite(analysis[[variable]])] <- value
        assessment[[variable]][!is.finite(assessment[[variable]])] <- value
        rules$median[[i]] <- value
    }
    list(analysis = analysis, assessment = assessment, rules = rules)
}

prepare_candidate <- function(analysis, assessment, programme, candidate) {
    history <- history_variables[[history_for_programme[[programme]]]]
    raw_variables <- unique(c(
        anchors, existing_modifiers, history, new_primary, new_sensitivity
    ))
    imputed <- impute_from_analysis(analysis, assessment, raw_variables)
    analysis <- imputed$analysis
    assessment <- imputed$assessment
    pca <- NULL

    if (candidate == 'core') {
        predictors <- c(anchors, existing_modifiers, history)
    } else if (candidate == 'sst_chla_raw') {
        predictors <- c(anchors, existing_modifiers, history, new_primary)
    } else {
        # PCA is restricted to environmental modifiers. Current DHW,
        # composition, cover and interpretable exposure history remain direct.
        pca_variables <- c(
            existing_modifiers, new_primary, new_sensitivity
        )
        pca <- prcomp(
            analysis[pca_variables], center = TRUE, scale. = TRUE
        )
        cumulative <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
        components <- min(which(cumulative >= 0.85))
        components <- min(components, 6L)
        analysis_scores <- predict(pca, analysis[pca_variables])[
            , seq_len(components), drop = FALSE
        ]
        pc_names <- paste0('modifier_pc', seq_len(components))
        assessment_scores <- if (nrow(assessment) > 0) {
            predict(pca, assessment[pca_variables])[
                , seq_len(components), drop = FALSE
            ]
        } else {
            matrix(
                numeric(), nrow = 0, ncol = components,
                dimnames = list(NULL, pc_names)
            )
        }
        colnames(analysis_scores) <- pc_names
        colnames(assessment_scores) <- pc_names
        analysis[pc_names] <- analysis_scores
        assessment[pc_names] <- assessment_scores
        predictors <- c(anchors, history, pc_names)
    }
    list(
        analysis = analysis, assessment = assessment,
        predictors = predictors, pca = pca
    )
}

fit_brt <- function(analysis, predictors) {
    occurrence_formula <- reformulate(predictors, response = 'has_loss')
    magnitude_formula <- reformulate(predictors, response = 'positive_logit')
    analysis$has_loss <- as.numeric(analysis$mortality_prop > 0)
    positive <- analysis |> filter(has_loss == 1)
    epsilon <- min(0.01, 0.5 / nrow(positive))
    positive$positive_logit <- qlogis(pmin(
        pmax(positive$mortality_prop, epsilon), 1 - epsilon
    ))
    settings <- list(
        n.trees = 1000L, interaction.depth = 2L, shrinkage = 0.02,
        n.minobsinnode = min(8L, max(2L, floor(nrow(analysis) / 12L))),
        bag.fraction = 0.7, train.fraction = 1, keep.data = FALSE,
        verbose = FALSE
    )
    list(
        occurrence = do.call(gbm, c(
            list(
                formula = occurrence_formula, data = analysis,
                distribution = 'bernoulli'
            ),
            settings
        )),
        magnitude = do.call(gbm, c(
            list(
                formula = magnitude_formula, data = positive,
                distribution = 'gaussian'
            ),
            settings
        ))
    )
}

predict_brt <- function(model, assessment) {
    occurrence <- predict(
        model$occurrence, assessment, n.trees = 1000, type = 'response'
    )
    magnitude <- plogis(predict(
        model$magnitude, assessment, n.trees = 1000, type = 'response'
    ))
    tibble(
        predicted_occurrence = pmin(pmax(occurrence, 0), 1),
        predicted_positive_mortality = pmin(pmax(magnitude, 0), 1),
        predicted_mortality = predicted_occurrence *
            predicted_positive_mortality
    )
}

standardise_for_inla <- function(analysis, assessment, predictors) {
    scaled <- character()
    for (predictor in predictors) {
        centre <- mean(analysis[[predictor]])
        spread <- sd(analysis[[predictor]])
        if (!is.finite(spread) || spread == 0) spread <- 1
        name <- paste0(predictor, '_z')
        analysis[[name]] <- (analysis[[predictor]] - centre) / spread
        assessment[[name]] <- (assessment[[predictor]] - centre) / spread
        scaled <- c(scaled, name)
    }
    list(analysis = analysis, assessment = assessment, predictors = scaled)
}

fit_predict_inla <- function(analysis, assessment, predictors) {
    prepared <- standardise_for_inla(analysis, assessment, predictors)
    analysis <- prepared$analysis
    assessment <- prepared$assessment
    predictors <- prepared$predictors

    levels <- bind_rows(analysis, assessment) |>
        transmute(ReefID, event_year, region_block)
    reef_levels <- unique(levels$ReefID)
    event_levels <- sort(unique(levels$event_year))
    region_levels <- unique(levels$region_block)
    add_indices <- function(rows) {
        rows |>
            mutate(
                reef_index = match(ReefID, reef_levels),
                event_index = match(event_year, event_levels),
                region_index = match(region_block, region_levels)
            )
    }
    analysis <- add_indices(analysis)
    assessment <- add_indices(assessment)
    random <- paste(
        "f(reef_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.7,0.05))))",
        "f(event_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.7,0.05))))",
        "f(region_index, model='iid', hyper=list(prec=list(prior='pc.prec', param=c(0.5,0.05))))",
        sep = ' + '
    )
    formula <- as.formula(paste(
        'response ~', paste(predictors, collapse = ' + '), '+', random
    ))
    occurrence_data <- bind_rows(
        analysis |> mutate(response = as.numeric(mortality_prop > 0)),
        assessment |> mutate(response = NA_real_)
    )
    occurrence <- inla(
        formula, family = 'binomial', Ntrials = 1,
        data = occurrence_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
        ),
        control.predictor = list(compute = TRUE, link = 1),
        verbose = FALSE
    )

    positive <- analysis |> filter(mortality_prop > 0)
    positive_n <- nrow(positive)
    positive$response <- (
        positive$mortality_prop * (positive_n - 1) + 0.5
    ) / positive_n
    magnitude_data <- bind_rows(
        positive, assessment |> mutate(response = NA_real_)
    )
    magnitude <- inla(
        formula, family = 'beta', data = magnitude_data,
        control.fixed = list(
            mean = 0, prec = 4, mean.intercept = 0, prec.intercept = 1
        ),
        control.family = list(
            hyper = list(theta = list(prior = 'loggamma', param = c(2, 0.1)))
        ),
        control.predictor = list(compute = TRUE, link = 1),
        verbose = FALSE
    )
    occurrence_rows <- nrow(analysis) + seq_len(nrow(assessment))
    magnitude_rows <- positive_n + seq_len(nrow(assessment))
    occurrence_mean <- occurrence$summary.fitted.values[
        occurrence_rows, 'mean'
    ]
    magnitude_mean <- magnitude$summary.fitted.values[
        magnitude_rows, 'mean'
    ]
    tibble(
        predicted_occurrence = occurrence_mean,
        predicted_positive_mortality = magnitude_mean,
        predicted_mortality = occurrence_mean * magnitude_mean
    )
}

predictions <- tibble()
runtime <- tibble()

for (programme in names(validation_files)) {
    rows <- load_screen_rows(programme)
    for (candidate in candidates) {
        for (scheme in schemes) {
            for (fold in fold_values(rows, scheme)) {
                assessment_rows_index <- assessment_rows(rows, scheme, fold)
                prepared <- prepare_candidate(
                    rows[!assessment_rows_index, , drop = FALSE],
                    rows[assessment_rows_index, , drop = FALSE],
                    programme, candidate
                )
                for (learner in c('brt', 'inla')) {
                    started <- proc.time()[['elapsed']]
                    predicted <- if (learner == 'brt') {
                        predict_brt(
                            fit_brt(prepared$analysis, prepared$predictors),
                            prepared$assessment
                        )
                    } else {
                        fit_predict_inla(
                            prepared$analysis, prepared$assessment,
                            prepared$predictors
                        )
                    }
                    elapsed <- proc.time()[['elapsed']] - started
                    predictions <- bind_rows(
                        predictions,
                        bind_cols(
                            prepared$assessment |>
                                transmute(
                                    programme_key = .env$programme,
                                    source_observation_id, ReefID, ReefName,
                                    event_year,
                                    observed_mortality = mortality_prop,
                                    observed_occurrence =
                                        as.numeric(mortality_prop > 0)
                                ),
                            predicted
                        ) |>
                            mutate(
                                candidate = candidate, learner = learner,
                                scheme = scheme, fold = as.character(fold)
                            )
                    )
                    runtime <- bind_rows(
                        runtime,
                        tibble(
                            programme_key = programme, candidate, learner,
                            scheme, fold = as.character(fold),
                            elapsed_seconds = elapsed
                        )
                    )
                }
                message(
                    programme, ' / ', candidate, ' / ', scheme, ' / ', fold
                )
            }
        }
    }
}

metrics <- predictions |>
    group_by(programme_key, candidate, learner, scheme) |>
    summarise(
        rows = n(),
        rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
        mae = mean(abs(observed_mortality - predicted_mortality)),
        predictive_r2 = 1 -
            sum((observed_mortality - predicted_mortality)^2) /
            sum((observed_mortality - mean(observed_mortality))^2),
        bias = mean(predicted_mortality - observed_mortality),
        occurrence_brier = mean(
            (observed_occurrence - predicted_occurrence)^2
        ),
        severe_rows = sum(observed_mortality >= 0.5),
        severe_bias = if_else(
            severe_rows > 0,
            mean(
                predicted_mortality[observed_mortality >= 0.5] -
                observed_mortality[observed_mortality >= 0.5]
            ),
            NA_real_
        ),
        .groups = 'drop'
    )

metrics_2024 <- predictions |>
    filter(scheme == 'reef_blocked_5fold', event_year == 2024) |>
    group_by(programme_key, candidate, learner) |>
    summarise(
        rows = n(),
        rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
        mae = mean(abs(observed_mortality - predicted_mortality)),
        predictive_r2 = 1 -
            sum((observed_mortality - predicted_mortality)^2) /
            sum((observed_mortality - mean(observed_mortality))^2),
        bias = mean(predicted_mortality - observed_mortality),
        severe_rows = sum(observed_mortality >= 0.5),
        severe_bias = if_else(
            severe_rows > 0,
            mean(
                predicted_mortality[observed_mortality >= 0.5] -
                observed_mortality[observed_mortality >= 0.5]
            ),
            NA_real_
        ),
        .groups = 'drop'
    )

# Full-data PCA audit. These are descriptive loadings; every predictive fold
# above learned its own imputation, scaling and rotation from analysis rows.
pca_loadings <- tibble()
pca_variance <- tibble()
for (programme in names(validation_files)) {
    rows <- load_screen_rows(programme)
    prepared <- prepare_candidate(
        rows, rows[0, , drop = FALSE], programme, 'modifier_pca'
    )
    pca <- prepared$pca
    variance <- pca$sdev^2 / sum(pca$sdev^2)
    retained <- grep('^modifier_pc', prepared$predictors, value = TRUE)
    number <- length(retained)
    pca_loadings <- bind_rows(
        pca_loadings,
        as.data.frame(pca$rotation[, seq_len(number), drop = FALSE]) |>
            tibble::rownames_to_column('feature') |>
            pivot_longer(-feature, names_to = 'component', values_to = 'loading') |>
            mutate(programme_key = programme)
    )
    pca_variance <- bind_rows(
        pca_variance,
        tibble(
            programme_key = programme,
            component = paste0('PC', seq_along(variance)),
            variance_explained = variance,
            cumulative_variance = cumsum(variance),
            retained = seq_along(variance) <= number
        )
    )
}

feature_rows <- bind_rows(lapply(names(validation_files), load_screen_rows))
correlation_variables <- c(
    'ann_maxdhw', existing_modifiers, new_primary, new_sensitivity
)
correlations <- cor(
    feature_rows[correlation_variables],
    use = 'pairwise.complete.obs', method = 'spearman'
) |>
    as.data.frame() |>
    tibble::rownames_to_column('feature_1') |>
    pivot_longer(-feature_1, names_to = 'feature_2', values_to = 'spearman_rho')

write_csv(predictions, file.path(output_dir, 'predictions.csv'), na = '')
write_csv(metrics, file.path(output_dir, 'metrics.csv'), na = '')
write_csv(metrics_2024, file.path(output_dir, 'metrics_2024.csv'), na = '')
write_csv(runtime, file.path(output_dir, 'runtime.csv'), na = '')
write_csv(pca_loadings, file.path(output_dir, 'pca_loadings.csv'), na = '')
write_csv(pca_variance, file.path(output_dir, 'pca_variance.csv'), na = '')
write_csv(correlations, file.path(output_dir, 'correlations.csv'), na = '')

print(metrics)
print(metrics_2024)
