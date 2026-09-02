# Audit the taxonomic/morphological resolution of the local AIMS community data.

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
})

workspace <- 'data/processed/01_exploratory_workspace.RData'
output <- 'data/processed/acropora_resolution_audit.csv'

load(workspace)
if (!exists('df.AIMS.full')) stop('df.AIMS.full is absent from ', workspace)

categories <- df.AIMS.full |>
    transmute(
        category = as.character(reefpage_category),
        source_variable = as.character(variable),
        reef_name = as.character(reef_name),
        report_year = as.integer(report_year)
    )

patterns <- tibble::tribble(
    ~composition_field, ~pattern,
    'Total Acropora', '^Acropora$',
    'Tabular Acropora', 'tabular|table',
    'Staghorn Acropora', 'staghorn',
    'Branching Acropora', 'branching.*Acropora|Acropora.*branching'
)

audit <- bind_rows(lapply(seq_len(nrow(patterns)), function(i) {
    matched <- categories |>
        filter(grepl(
            patterns$pattern[[i]], category, ignore.case = TRUE
        ) | grepl(
            patterns$pattern[[i]], source_variable, ignore.case = TRUE
        ))
    tibble(
        composition_field = patterns$composition_field[[i]],
        available = nrow(matched) > 0,
        rows = nrow(matched),
        reefs = n_distinct(matched$reef_name),
        first_year = if (nrow(matched)) min(matched$report_year, na.rm = TRUE) else NA_integer_,
        last_year = if (nrow(matched)) max(matched$report_year, na.rm = TRUE) else NA_integer_,
        source_categories = paste(sort(unique(matched$category)), collapse = '; '),
        source_variables = paste(sort(unique(matched$source_variable)), collapse = '; ')
    )
}))

write_csv(audit, output, na = '')
print(audit)
