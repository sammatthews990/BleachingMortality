suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(brms)
  library(patchwork)
})

# Load the cached environments
e_m <- new.env()
lazyLoad('DHW_Mortality_Attribution_cache/html/manta-model-data_d3d5fd24b8f548ff1210864d2d490674', envir = e_m)
dat.mod.mant <- e_m[['dat.mod.mant']] %>% st_drop_geometry()

e_b <- new.env()
lazyLoad('DHW_Mortality_Attribution_cache/html/benthic-data_3436ddf577b9428741a2f0fe11c82b03', envir = e_b)
dat.mod.bent.hc <- e_b[['dat.mod.bent.hc']] %>% st_drop_geometry()

zarr_full_raw <- read.csv('data/processed/cheung_recreated_gbr_full.csv')

d_m_all <- dat.mod.mant %>%
    ungroup() %>%
    filter(!is.na(ReefID), !is.na(DHWYear)) %>%
    mutate(DHWYear = as.numeric(DHWYear),
           ReefID_clean = toupper(trimws(as.character(ReefID)))) %>%
    left_join(
        zarr_full_raw %>% mutate(LABEL_ID_clean = toupper(trimws(as.character(LABEL_ID)))),
        by = c('ReefID_clean' = 'LABEL_ID_clean', 'DHWYear' = 'year')
    ) %>%
    mutate(Dataset = 'Manta Tow', survey_depth = 7.5)

d_l_all <- dat.mod.bent.hc %>%
    ungroup() %>%
    filter(project_code == 'LTMP', !is.na(ReefID), !is.na(DHWYear)) %>%
    mutate(DHWYear = as.numeric(DHWYear),
           ReefID_clean = toupper(trimws(as.character(ReefID)))) %>%
    left_join(
        zarr_full_raw %>% mutate(LABEL_ID_clean = toupper(trimws(as.character(LABEL_ID)))),
        by = c('ReefID_clean' = 'LABEL_ID_clean', 'DHWYear' = 'year')
    ) %>%
    mutate(Dataset = 'LTMP Benthic')

d_p_all <- dat.mod.bent.hc %>%
    ungroup() %>%
    filter(project_code == 'MMP', !is.na(ReefID), !is.na(DHWYear)) %>%
    mutate(DHWYear = as.numeric(DHWYear),
           ReefID_clean = toupper(trimws(as.character(ReefID)))) %>%
    left_join(
        zarr_full_raw %>% mutate(LABEL_ID_clean = toupper(trimws(as.character(LABEL_ID)))),
        by = c('ReefID_clean' = 'LABEL_ID_clean', 'DHWYear' = 'year')
    ) %>%
    mutate(Dataset = 'MMP Inshore')

dat.ml.expanded <- bind_rows(d_m_all, d_l_all, d_p_all) %>%
    mutate(MaxDHW_val = coalesce(ann_maxdhw, MaxDHW.mean)) %>%
    filter(!is.na(Rel.Change), !is.na(MaxDHW_val))

# Define Era: 3 eras
dat.ml.expanded <- dat.ml.expanded %>%
    mutate(
        Era = case_when(
            DHWYear <= 2002 ~ "Historical (1998-2002)",
            DHWYear >= 2016 & DHWYear <= 2020 ~ "Modern Baseline (2016-2020)",
            DHWYear >= 2022 ~ "Contemporary (2022-2024)",
            TRUE ~ "Intermediate"
        ),
        Era = factor(Era, levels = c("Historical (1998-2002)", "Modern Baseline (2016-2020)", "Contemporary (2022-2024)"))
    ) %>%
    filter(!is.na(Era))

# Also define discrete Major Event Year (1998, 2002, 2016, 2017, 2020, 2022, 2024)
dat.ml.expanded <- dat.ml.expanded %>%
    mutate(
        EventYear = as.character(DHWYear),
        EventYear = factor(EventYear, levels = c("1998", "2002", "2016", "2017", "2020", "2022", "2024"))
    )

# Prepare scaled variables
N_all <- nrow(dat.ml.expanded)
dat.ml.expanded <- dat.ml.expanded %>%
    mutate(
        Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
        Mort.prop.nudge = (Mort.prop * (N_all - 1) + 0.5) / N_all,
        MaxDHW_s = (MaxDHW_val - mean(MaxDHW_val, na.rm=TRUE)) / (2 * sd(MaxDHW_val, na.rm=TRUE)),
        Year_s = (DHWYear - 2016) / 10,
        ReefName_clean = coalesce(as.character(Reef_Name), as.character(ReefName), as.character(ReefID))
    )

cat('Clean dataset ready with N =', nrow(dat.ml.expanded), '\n')

# Quick GLM test with MaxDHW_s * Era
m_era_glm <- glm(Mort.prop ~ MaxDHW_s * Era + Dataset, data = dat.ml.expanded, family = quasibinomial)
cat('\n=== GLM with DHW * Era ===\n')
print(summary(m_era_glm))

# Quick GLM test with MaxDHW_s * Year_s (Continuous time trend)
m_time_glm <- glm(Mort.prop ~ MaxDHW_s * Year_s + Dataset, data = dat.ml.expanded, family = quasibinomial)
cat('\n=== GLM with DHW * Continuous Year ===\n')
print(summary(m_time_glm))

# Quick GLM test with MaxDHW_s * EventYear (Event-specific sensitivity)
dat_events <- dat.ml.expanded %>% filter(!is.na(EventYear))
m_event_glm <- glm(Mort.prop ~ MaxDHW_s * EventYear + Dataset, data = dat_events, family = quasibinomial)
cat('\n=== GLM with DHW * Event Year ===\n')
print(summary(m_event_glm))
