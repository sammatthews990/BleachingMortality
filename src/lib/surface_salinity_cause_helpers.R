# Classify field disturbance descriptions for the surface-salinity screen.
#
# The source flags are reef/report annotations, not independently verified
# causal assignments. In particular, every 2024 flood annotation is also
# cyclone-labelled. The flood-inclusive class therefore permits a cyclone
# co-label only when flood/freshwater/low-salinity wording is explicit.

classify_surface_salinity_causes <- function(rows) {
    required <- c(
        'ReefID', 'event_year', 'disturbance_text',
        'disturbance_has_bleaching', 'disturbance_has_cyclone',
        'disturbance_has_flood', 'disturbance_has_cots'
    )
    missing <- setdiff(required, names(rows))
    if (length(missing) > 0L) {
        stop('Missing disturbance fields: ', paste(missing, collapse = ', '))
    }

    rows |>
        dplyr::mutate(
            disturbance_text_clean = trimws(dplyr::coalesce(
                as.character(disturbance_text), ''
            )),
            disturbance_has_bleaching = dplyr::coalesce(
                as.logical(disturbance_has_bleaching), FALSE
            ),
            disturbance_has_cyclone = dplyr::coalesce(
                as.logical(disturbance_has_cyclone), FALSE
            ),
            disturbance_has_flood = dplyr::coalesce(
                as.logical(disturbance_has_flood), FALSE
            ),
            disturbance_has_cots = dplyr::coalesce(
                as.logical(disturbance_has_cots), FALSE
            ),
            disturbance_has_disease = grepl(
                'disease', disturbance_text_clean, ignore.case = TRUE
            ),
            target_description_class = dplyr::case_when(
                disturbance_has_cots ~ 'excluded_cots',
                disturbance_has_disease ~ 'excluded_disease',
                disturbance_has_flood ~ 'flooding_recorded',
                disturbance_has_cyclone ~
                    'excluded_cyclone_without_flood',
                disturbance_has_bleaching ~ 'bleaching_only',
                disturbance_text_clean == '' ~
                    'no_recorded_disturbance',
                TRUE ~ 'excluded_unclassified_text'
            ),
            target_description_strict_row = target_description_class %in% c(
                'bleaching_only', 'no_recorded_disturbance'
            ),
            target_description_flood_row = target_description_class %in% c(
                'bleaching_only', 'flooding_recorded',
                'no_recorded_disturbance'
            )
        ) |>
        dplyr::group_by(ReefID, event_year) |>
        dplyr::mutate(
            target_description_strict_reef = all(
                target_description_strict_row
            ),
            target_description_flood_reef = all(
                target_description_flood_row
            ),
            reef_has_flood_description = any(
                target_description_class == 'flooding_recorded'
            ),
            reef_has_competing_description = any(
                target_description_class %in% c(
                    'excluded_cots', 'excluded_disease',
                    'excluded_cyclone_without_flood',
                    'excluded_unclassified_text'
                )
            )
        ) |>
        dplyr::ungroup()
}

validate_surface_salinity_cause_classifier <- function() {
    example <- data.frame(
        ReefID = LETTERS[1:6], event_year = 2024L,
        disturbance_text = c(
            '', 'coral bleaching', 'flood waters, Cyclone Jasper',
            'coral bleaching, Cyclone Jasper',
            'crown-of-thorns starfish', 'DHW 8.1'
        ),
        disturbance_has_bleaching = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE),
        disturbance_has_cyclone = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE),
        disturbance_has_flood = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE),
        disturbance_has_cots = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)
    )
    observed <- classify_surface_salinity_causes(example)$target_description_class
    expected <- c(
        'no_recorded_disturbance', 'bleaching_only', 'flooding_recorded',
        'excluded_cyclone_without_flood', 'excluded_cots',
        'excluded_unclassified_text'
    )
    if (!identical(observed, expected)) {
        stop('Surface-salinity cause classifier self-check failed')
    }
    invisible(TRUE)
}
