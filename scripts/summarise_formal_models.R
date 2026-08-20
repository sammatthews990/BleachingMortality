# Export compact summaries from the cached full-data BRMS and BRT models.

suppressPackageStartupMessages({
    library(brms)
    library(loo)
    library(dplyr)
    library(gbm)
    library(readr)
    library(tibble)
})

model_dir <- "output/models/formal"
output_dir <- "output/formal_models"
programme_keys <- c("manta", "ltmp", "mmp")

brt_influence <- tibble()
brt_settings <- tibble()
brms_coefficients <- tibble()
brms_production_diagnostics <- tibble()
brms_fit_statistics <- tibble()

for (programme_key in programme_keys) {
    brt_cache <- readRDS(file.path(
        model_dir, paste0("brt_", programme_key, "_production.rds")
    ))
    for (component in c("occurrence", "magnitude")) {
        component_model <- brt_cache$model[[component]]
        if (!identical(component_model$type, "gbm")) next
        influence <- summary(component_model$model, plotit = FALSE) |>
            as_tibble() |>
            transmute(
                programme_key,
                component,
                predictor = var,
                relative_influence = rel.inf
            )
        brt_influence <- bind_rows(brt_influence, influence)
    }
    brt_settings <- bind_rows(
        brt_settings,
        brt_cache$tuning$best |>
            mutate(programme_key, .before = 1)
    )

    brms_cache <- readRDS(file.path(
        model_dir, paste0("brms_", programme_key, "_production.rds")
    ))
    marginal_r2 <- bayes_R2(
        brms_cache$fit, re_formula = NA,
        summary = TRUE, probs = c(0.05, 0.95)
    )
    conditional_r2 <- bayes_R2(
        brms_cache$fit, re_formula = NULL,
        summary = TRUE, probs = c(0.05, 0.95)
    )
    loo_result <- loo(brms_cache$fit)
    pareto_k <- pareto_k_values(loo_result)
    brms_fit_statistics <- bind_rows(
        brms_fit_statistics,
        tibble(
            programme_key,
            marginal_bayesian_r2 = marginal_r2[1, 'Estimate'],
            marginal_bayesian_r2_q05 = marginal_r2[1, 'Q5'],
            marginal_bayesian_r2_q95 = marginal_r2[1, 'Q95'],
            conditional_bayesian_r2 = conditional_r2[1, 'Estimate'],
            conditional_bayesian_r2_q05 = conditional_r2[1, 'Q5'],
            conditional_bayesian_r2_q95 = conditional_r2[1, 'Q95'],
            elpd_loo = loo_result$estimates['elpd_loo', 'Estimate'],
            elpd_loo_se = loo_result$estimates['elpd_loo', 'SE'],
            p_loo = loo_result$estimates['p_loo', 'Estimate'],
            looic = loo_result$estimates['looic', 'Estimate'],
            looic_se = loo_result$estimates['looic', 'SE'],
            pareto_k_over_0_7 = sum(pareto_k > 0.7, na.rm = TRUE),
            max_pareto_k = max(pareto_k, na.rm = TRUE)
        )
    )

    coefficient_matrix <- brms::fixef(
        brms_cache$fit,
        probs = c(0.05, 0.95),
        robust = FALSE
    )
    brms_coefficients <- bind_rows(
        brms_coefficients,
        as_tibble(coefficient_matrix, rownames = "term") |>
            transmute(
                programme_key,
                term,
                estimate = Estimate,
                std_error = Est.Error,
                q05 = Q5,
                q95 = Q95
            )
    )
    brms_production_diagnostics <- bind_rows(
        brms_production_diagnostics,
        brms_cache$diagnostics |>
            mutate(
                programme_key,
                n = brms_cache$n,
                event_years = paste(brms_cache$event_years, collapse = ","),
                .before = 1
            )
    )
}

write_csv(brt_influence, file.path(output_dir, "brt_production_influence.csv"))
write_csv(brt_settings, file.path(output_dir, "brt_production_settings.csv"))
write_csv(
    brms_coefficients,
    file.path(output_dir, "brms_production_coefficients.csv")
)
write_csv(
    brms_production_diagnostics,
    file.path(output_dir, "brms_production_diagnostics.csv")
)

write_csv(
    brms_fit_statistics,
    file.path(output_dir, 'brms_fit_statistics.csv')
)

print(brt_influence |>
    group_by(programme_key, component) |>
    slice_max(relative_influence, n = 5, with_ties = FALSE) |>
    arrange(programme_key, component, desc(relative_influence)))
print(brms_production_diagnostics)
