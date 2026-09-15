# Build the small site-event grid used to compare satellite cooling with
# subdaily AIMS water-quality logger summaries.
suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

input <- 'data/processed/aims_wq_logger_hourly_event_summary.csv'
output <- 'data/processed/aims_wq_logger_thermal_grid.csv'
if (!file.exists(input)) {
    stop('Missing ', input, '; run src/data/extract_aims_logger_hourly.R')
}

grid <- read_csv(input, show_col_types = FALSE) |>
    filter(event_year >= 2012L, event_year <= 2024L) |>
    transmute(
        ReefID = site,
        year = as.integer(event_year),
        lon = site_longitude,
        lat = site_latitude
    ) |>
    distinct() |>
    arrange(year, ReefID)

if (anyDuplicated(grid[c('ReefID', 'year')])) {
    stop('Logger thermal grid has duplicate site-event keys')
}
write_csv(grid, output)
cat('Wrote', nrow(grid), 'logger site-events to', output, '\n')
