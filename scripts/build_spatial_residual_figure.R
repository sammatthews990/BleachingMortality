# Six-panel spatial residual diagnostic: five event-specific fields from the
# event-spatial AR1 screening model plus the pooled persistent field from the
# selected operational model. Coordinates are reported in decimal degrees.

suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(readr)
    library(scales)
    library(sf)
    library(tibble)
})
source('scripts/model_registry.R')
source('scripts/model_diagnostics.R')

root <- project_root()
registry <- read_model_registry(root)
selected_id <- registry$best_model$id
model_dir <- file.path(root, 'output', 'inla_spatiotemporal')

event_path <- file.path(model_dir, 'full_v2_event_spatial_ar1.rds')
selected_path <- resolve_fit_artifact(selected_id, root)
mesh <- readRDS(file.path(model_dir, 'mesh.rds'))
event_fit <- readRDS(event_path)$fit
selected_fit <- readRDS(selected_path)$fit

event_years <- c(2016L, 2017L, 2020L, 2022L, 2024L)
event_field <- event_fit$summary.random$event_spatial_field |>
    as_tibble() |>
    mutate(
        node = rep(seq_len(mesh$n), times = length(event_years)),
        panel = rep(as.character(event_years), each = mesh$n),
        field_source = 'Event-spatial AR1 screening model'
    )
overall_field <- selected_fit$summary.random$spatial_field |>
    as_tibble() |>
    mutate(
        node = seq_len(mesh$n),
        panel = 'Overall',
        field_source = paste('Selected persistent field:', selected_id)
    )

node_sf <- st_as_sf(
    tibble(
        node = seq_len(mesh$n),
        x_m = mesh$loc[, 1] * 1000,
        y_m = mesh$loc[, 2] * 1000
    ),
    coords = c('x_m', 'y_m'), crs = 3577
) |>
    st_transform(4326)
coordinates <- st_coordinates(node_sf) |>
    as_tibble() |>
    transmute(node = seq_len(n()), longitude = X, latitude = Y)

fields <- bind_rows(event_field, overall_field) |>
    left_join(coordinates, by = 'node') |>
    mutate(
        panel = factor(
            panel,
            levels = c(as.character(event_years), 'Overall')
        ),
        lower95 = .data[['0.025quant']],
        median = .data[['0.5quant']],
        upper95 = .data[['0.975quant']]
    ) |>
    select(
        panel, node, longitude, latitude, mean, sd,
        lower95, median, upper95, field_source
    )

colour_limit <- max(abs(fields$mean), na.rm = TRUE)
plot <- ggplot(fields, aes(longitude, latitude, colour = mean)) +
    geom_point(size = 1.35, alpha = 0.9) +
    facet_wrap(~ panel, ncol = 3) +
    coord_equal() +
    scale_x_continuous(
        labels = label_number(accuracy = 0.1, suffix = '\u00b0E'),
        breaks = breaks_pretty(n = 4)
    ) +
    scale_y_continuous(
        labels = label_number(accuracy = 0.1, suffix = '\u00b0'),
        breaks = breaks_pretty(n = 4)
    ) +
    scale_colour_gradient2(
        low = '#2166AC', mid = 'white', high = '#B2182B', midpoint = 0,
        limits = c(-colour_limit, colour_limit)
    ) +
    labs(
        x = 'Longitude', y = 'Latitude',
        colour = 'Residual effect\n(link scale)',
        title = 'Event-specific and overall spatial residual fields',
        subtitle = 'Red: mortality higher than covariates explain; blue: lower than covariates explain'
    ) +
    theme_bw(base_size = 11) +
    theme(
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = 'bold'),
        legend.position = 'right'
    )

write_csv(
    fields,
    file.path(model_dir, 'spatial_residual_fields_decimal_degrees.csv')
)
save_figure_bundle(
    plot, 'Fig-INLA-07_spatial_residual_fields', fields,
    'Posterior mean event-specific spatial residual fields for 2016, 2017, 2020, 2022 and 2024, with the selected model persistent field shown as the overall panel.',
    'Maps residual geographic structure after measured covariates. Repeated red or blue areas suggest unresolved persistent mechanisms; event-only patterns suggest event-specific processes or observation gaps.',
    'The five event panels come from the event-spatial AR1 screening model, whereas Overall is the selected operational model persistent field. These are latent link-scale effects, not observed mortality or causal attribution.',
    selected_id, 'INLA', 'spatial_residual_field',
    'spatial_model_diagnostic', root, TRUE,
    code_source = 'scripts/build_spatial_residual_figure.R',
    width = 12, height = 8.5
)
write_figure_readme(root)
message('Wrote six-panel residual field figure')
