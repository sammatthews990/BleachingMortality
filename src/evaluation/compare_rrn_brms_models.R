# Compare the matched RRN BRMS models with and without modelled COTS pressure.

suppressPackageStartupMessages({
    library(dplyr)
    library(loo)
    library(readr)
})

model_dir <- 'output/models/rrn_pressure'
output_dir <- 'output/rrn_pressure_assessment'
models <- c('weather_rrn', 'weather_rrn_cots')

caches <- lapply(models, function(model) {
    readRDS(file.path(
        model_dir, paste0('rrn_pressure_brms_', model, '.rds')
    ))
})
names(caches) <- models

comparison <- as.data.frame(loo_compare(lapply(caches, `[[`, 'loo'))) |>
    tibble::rownames_to_column('model') |>
    as_tibble()

fit_metrics <- bind_rows(lapply(models, function(model) {
    read_csv(
        file.path(output_dir, paste0('rrn_brms_fit_', model, '.csv')),
        show_col_types = FALSE
    )
})) |>
    left_join(comparison, by = 'model')

effects <- bind_rows(lapply(models, function(model) {
    read_csv(
        file.path(output_dir, paste0('rrn_brms_effects_', model, '.csv')),
        show_col_types = FALSE
    )
}))
diagnostics <- bind_rows(lapply(models, function(model) {
    read_csv(
        file.path(output_dir, paste0('rrn_brms_diagnostics_', model, '.csv')),
        show_col_types = FALSE
    )
}))

write_csv(fit_metrics, file.path(output_dir, 'rrn_brms_model_comparison.csv'))
write_csv(effects, file.path(output_dir, 'rrn_brms_effects.csv'))
write_csv(diagnostics, file.path(output_dir, 'rrn_brms_diagnostics.csv'))

print(fit_metrics)
