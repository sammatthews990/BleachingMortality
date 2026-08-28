# Build the small reef-event grid used to test new remotely sensed features.
# The full extractor can also use the all-GBR environmental grid directly, but
# this validation grid avoids downloading values that cannot yet inform model
# comparison.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})
source('scripts/formal_model_helpers.R')

output_path <- 'data/processed/environment_screen_grid_validation.csv'

grid <- bind_rows(lapply(names(validation_files), function(programme) {
    load_programme_rows(programme) |>
        transmute(
            programme_key = .env$programme,
            ReefID, ReefName, year = event_year,
            lon = as.numeric(lon), lat = as.numeric(lat)
        )
})) |>
    filter(is.finite(lon), is.finite(lat)) |>
    distinct(ReefID, year, lon, lat, .keep_all = TRUE) |>
    arrange(year, ReefID)

coordinate_conflicts <- grid |>
    count(ReefID, year, name = 'coordinate_rows') |>
    filter(coordinate_rows > 1)
if (nrow(coordinate_conflicts) > 0) {
    stop('A reef-year maps to multiple validation coordinates')
}

write_csv(grid, output_path, na = '')
cat(
    'Wrote', nrow(grid), 'validation reef-years across',
    n_distinct(grid$ReefID), 'reefs to', output_path, '\n'
)
