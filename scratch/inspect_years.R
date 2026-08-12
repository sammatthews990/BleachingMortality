suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
})
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

df_counts <- dat.ml.expanded %>%
    group_by(DHWYear, Dataset) %>%
    summarise(
        n_surveys = n(),
        mean_dhw = round(mean(MaxDHW_val, na.rm=TRUE), 2),
        max_dhw = round(max(MaxDHW_val, na.rm=TRUE), 2),
        mean_mort = round(mean(pmax(-Rel.Change, 0), na.rm=TRUE), 3),
        .groups = 'drop'
    )

print(as.data.frame(df_counts))

cat('\nTotal by DHWYear:\n')
df_yr <- dat.ml.expanded %>%
    group_by(DHWYear) %>%
    summarise(
        n_surveys = n(),
        mean_dhw = round(mean(MaxDHW_val, na.rm=TRUE), 2),
        max_dhw = round(max(MaxDHW_val, na.rm=TRUE), 2),
        mean_mort = round(mean(pmax(-Rel.Change, 0), na.rm=TRUE), 3),
        .groups = 'drop'
    )
print(as.data.frame(df_yr))
