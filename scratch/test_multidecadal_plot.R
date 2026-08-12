suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(patchwork)
  library(brms)
})

e_m <- new.env()
lazyLoad('DHW_Mortality_Attribution_cache/html/manta-model-data_d3d5fd24b8f548ff1210864d2d490674', envir = e_m)
dat.mod.mant <- e_m[['dat.mod.mant']] %>% st_drop_geometry()

e_b <- new.env()
lazyLoad('DHW_Mortality_Attribution_cache/html/benthic-data_3436ddf577b9428741a2f0fe11c82b03', envir = e_b)
dat.mod.bent.hc <- e_b[['dat.mod.bent.hc']] %>% st_drop_geometry()

e_br <- new.env()
lazyLoad('DHW_Mortality_Attribution_cache/html/brms-combined-model_8b1a9744b367cd90c612986cdf406444', envir = e_br)
fit_brms_comb <- e_br[['fit_brms_comb']]

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
    filter(!is.na(Rel.Change))

# Scale params comb from cache or calculated
scale_params_comb <- list(
    MaxDHW.mean = c(mean = mean(dat.ml.expanded[['MaxDHW.mean']], na.rm=TRUE), sd2 = 2 * sd(dat.ml.expanded[['MaxDHW.mean']], na.rm=TRUE)),
    secc3m = c(mean = 17.7, sd2 = 11.0),
    cloudp_90 = c(mean = 0.72, sd2 = 0.2),
    prop_acropora = c(mean = 0.35, sd2 = 0.3)
)

dat.ml.expanded[['MaxDHW_val']] <- ifelse(!is.na(dat.ml.expanded[['ann_maxdhw']]), 
                                          dat.ml.expanded[['ann_maxdhw']], 
                                          dat.ml.expanded[['MaxDHW.mean']])

N_exp <- nrow(dat.ml.expanded)
dat.ml.expanded[['Mort.prop']] <- pmin(pmax(-dat.ml.expanded[['Rel.Change']], 0), 1)
dat.ml.expanded[['Dataset']] <- factor(dat.ml.expanded[['Dataset']], levels = c('Manta Tow', 'LTMP Benthic', 'MMP Inshore'))

d_exp <- dat.ml.expanded

dhw_ext_exp <- max(d_exp[['MaxDHW_val']], na.rm = TRUE)
es_exp <- annotate('rect', xmin = dhw_ext_exp, xmax = 20, ymin = -Inf, ymax = Inf, fill = 'grey80', alpha = 0.35)
el_exp <- annotate('text', x = (dhw_ext_exp + 20) / 2, y = 0.05, label = 'Extrapolation', fontface = 'italic', colour = 'grey40', size = 3)

grid_a_exp <- data.frame(
    MaxDHW_s = (seq(0, 20, length.out = 200) - scale_params_comb$MaxDHW.mean['mean']) / scale_params_comb$MaxDHW.mean['sd2'],
    secc3m_s = 0, cloudp_90_s = 0, prop_acropora_s = 0, winyear_sd_s = 0,
    histmDHW6_s = 0, mcur_90_s = 0, winyear_mean_s = 0, yrsince6_s = 0,
    Dataset = factor('Manta Tow', levels = levels(d_exp$Dataset))
)
grid_a_exp$MaxDHW_val <- seq(0, 20, length.out = 200)

pp_a_exp <- posterior_epred(fit_brms_comb, newdata = grid_a_exp, re_formula = NA)
grid_a_exp$pred <- apply(pp_a_exp, 2, median)
grid_a_exp$pred_lo <- apply(pp_a_exp, 2, quantile, 0.05)
grid_a_exp$pred_hi <- apply(pp_a_exp, 2, quantile, 0.95)

pA_exp <- ggplot() +
    es_exp + el_exp +
    geom_point(data = d_exp, aes(x = MaxDHW_val, y = Mort.prop, colour = Dataset), alpha = 0.35, size = 1.2) +
    geom_ribbon(data = grid_a_exp, aes(x = MaxDHW_val, ymin = pred_lo, ymax = pred_hi), fill = '#d95f02', alpha = 0.2) +
    geom_line(data = grid_a_exp, aes(x = MaxDHW_val, y = pred), colour = '#d95f02', linewidth = 1) +
    scale_colour_brewer(palette = 'Set2') +
    labs(x = 'Max DHW (°C-weeks)', y = 'Mortality Proportion', title = '(A) Multi-Decadal DHW Dose-Response (1998–2024)') +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 11)

ggsave('output/plots/test_multidecadal_pA.png', pA_exp, width = 6, height = 5, dpi = 150)
cat('Successfully generated output/plots/test_multidecadal_pA.png!\n')
