# Shared, deliberately small helpers for the formal BRMS/BRT prediction models.
# The modelling scripts remain linear and readable; this file only centralises
# fold definitions, predictor preprocessing, and repeated metric calculations.

suppressPackageStartupMessages({
    library(dplyr)
    library(splines)
})

validation_files <- c(
    manta = "data/processed/validation_rows_manta.rds",
    ltmp = "data/processed/validation_rows_ltmp.rds",
    mmp = "data/processed/validation_rows_mmp.rds"
)

environmental_predictors <- c(
    'prop_acropora_pre', 'acropora_cover_pre',
    "ann_maxdhw",
    "histmDHW6", "yrsince6",
    "histmDHW4", "yrsince4",
    "ann_maxsst", "winyear_mean", "winyear_sd",
    "mcur_90", "dist_to_er_km",
    "secc3m", "cloudp_90"
)
environmental_predictors <- c(
    environmental_predictors,
    'dhw10_load4', 'dhw_novelty10', 'secc3m_p10'
)

selection_schemes <- c("leave_one_event_out", "reef_blocked_5fold")
all_validation_schemes <- c(selection_schemes, "region_blocked")

programme_predictors <- function(programme_key) {
    if (programme_key == "mmp") {
        c(environmental_predictors, "depth")
    } else {
        environmental_predictors
    }
}

load_programme_rows <- function(programme_key) {
    path <- validation_files[[programme_key]]
    if (is.null(path) || !file.exists(path)) {
        stop("Missing validation rows for ", programme_key,
             ". Run scripts/build_validation_splits.R first.")
    }
    readRDS(path) |>
        arrange(event_year, ReefID, source_observation_id) |>
        mutate(
            reef_effect = factor(ReefID),
            event_effect = factor(event_year),
            region_effect = factor(region_block)
        )
}

fold_values <- function(data, scheme) {
    switch(
        scheme,
        leave_one_event_out = sort(unique(data$event_year)),
        reef_blocked_5fold = sort(unique(data$reef_fold)),
        region_blocked = sort(unique(data$region_block)),
        stop("Unknown validation scheme: ", scheme)
    )
}

assessment_rows <- function(data, scheme, fold) {
    switch(
        scheme,
        leave_one_event_out = data$event_year == as.integer(fold),
        reef_blocked_5fold = data$reef_fold == as.integer(fold),
        region_blocked = data$region_block == as.character(fold),
        stop("Unknown validation scheme: ", scheme)
    )
}

prepare_fold_predictors <- function(analysis, assessment, predictors) {
    preprocessing <- tibble(
        predictor = predictors,
        median = NA_real_,
        mean = NA_real_,
        sd = NA_real_,
        zero_variance = FALSE
    )

    for (i in seq_along(predictors)) {
        predictor <- predictors[[i]]
        observed <- analysis[[predictor]][is.finite(analysis[[predictor]])]
        if (length(observed) == 0L) {
            stop("No finite analysis values for ", predictor)
        }

        impute_value <- median(observed)
        analysis[[predictor]][!is.finite(analysis[[predictor]])] <- impute_value
        assessment[[predictor]][!is.finite(assessment[[predictor]])] <- impute_value

        centre <- mean(analysis[[predictor]])
        scale <- sd(analysis[[predictor]])
        zero_variance <- !is.finite(scale) || scale == 0
        if (zero_variance) scale <- 1

        scaled_name <- paste0(predictor, "_z")
        analysis[[scaled_name]] <- (analysis[[predictor]] - centre) / scale
        assessment[[scaled_name]] <- (assessment[[predictor]] - centre) / scale

        preprocessing$median[[i]] <- impute_value
        preprocessing$mean[[i]] <- centre
        preprocessing$sd[[i]] <- scale
        preprocessing$zero_variance[[i]] <- zero_variance
    }

    list(
        analysis = analysis,
        assessment = assessment,
        preprocessing = preprocessing
    )
}

make_balanced_inner_folds <- function(data, folds = 4L) {
    reef_sizes <- data |>
        count(ReefID, name = "observations") |>
        arrange(desc(observations), ReefID)
    number_of_folds <- min(folds, nrow(reef_sizes))
    fold_load <- rep(0L, number_of_folds)
    reef_sizes$inner_fold <- NA_integer_

    for (i in seq_len(nrow(reef_sizes))) {
        chosen <- which.min(fold_load)
        reef_sizes$inner_fold[[i]] <- chosen
        fold_load[[chosen]] <- fold_load[[chosen]] + reef_sizes$observations[[i]]
    }

    reef_sizes$inner_fold[match(data$ReefID, reef_sizes$ReefID)]
}

brms_rhs <- function(predictors) {
    scaled <- paste0(predictors, "_z")
    main_effects <- c(
        "splines::ns(ann_maxdhw_z, df = 3)",
        setdiff(scaled, "ann_maxdhw_z")
    )
    interactions <- c(
        'ann_maxdhw_z:prop_acropora_pre_z',
        "ann_maxdhw_z:cloudp_90_z",
        "ann_maxdhw_z:secc3m_z",
        "ann_maxdhw_z:mcur_90_z",
        "ann_maxdhw_z:histmDHW6_z"
    )
    interactions <- c(
        interactions,
        'ann_maxdhw_z:dhw10_load4_z',
        'ann_maxdhw_z:dhw_novelty10_z',
        'ann_maxdhw_z:secc3m_p10_z'
    )
    paste(
        c(
            main_effects,
            interactions,
            "(1 | reef_effect)",
            "(1 | event_effect)",
            "(1 | region_effect)"
        ),
        collapse = " + "
    )
}

prediction_metrics <- function(data) {
    has_occurrence <- 'predicted_occurrence' %in% names(data)
    data |>
        group_by(programme_key, learner, scheme) |>
        summarise(
            n = n(),
            mae = mean(abs(observed_mortality - predicted_mortality)),
            rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
            predictive_r2 = 1 -
                sum((observed_mortality - predicted_mortality)^2) /
                sum((observed_mortality - mean(observed_mortality))^2),
            bias = mean(predicted_mortality - observed_mortality),
            observed_mean = mean(observed_mortality),
            predicted_mean = mean(predicted_mortality),
            zero_mae = mean(
                abs(observed_mortality - predicted_mortality)[
                    observed_mortality == 0
                ]
            ),
            positive_mae = mean(
                abs(observed_mortality - predicted_mortality)[
                    observed_mortality > 0
                ]
            ),
            balanced_mae = mean(c(zero_mae, positive_mae)),
            severe_mae = mean(
                abs(observed_mortality - predicted_mortality)[
                    observed_mortality >= 0.2
                ]
            ),
            occurrence_brier = if (has_occurrence) mean(
                ((observed_mortality > 0) - predicted_occurrence)^2
            ) else NA_real_,
            .groups = "drop"
        )
}

safe_fold_label <- function(fold) {
    gsub("[^A-Za-z0-9_-]+", "_", as.character(fold))
}
