library(dplyr)
library(tidyr)
library(brms)

# Load data
load("Data/aims_ltmp/aims_ltmp.RData")
load("data/Cheungetal2025/01_sstvar_blchrf.RData")

# 1. Recreate cheung_predictors
predictor_vars <- c("mcur_90", "cloudp_90", "secc3m", "cbclus2", "histmDHW6", "yrsince6", "histmDHW4", "yrsince4", "winyear_sd", "winyear_mean")
cheung_predictors <- sstvar_blch2 %>%
  select(LABEL_id, year, lat = Y, lon = X, Sector, all_of(predictor_vars)) %>%
  group_by(LABEL_id, year) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(year = as.numeric(year))

# 2. Recreate dat.ml.ltmp and dat.ml.mmp
dat.mod.bent.hc <- df.BENT.HC %>%
  ungroup() %>%
  filter(!is.na(ReefID), !is.na(DHWYear)) %>%
  mutate(DHWYear = as.numeric(DHWYear))

dat.ml.ltmp <- dat.mod.bent.hc %>%
  filter(project_code == "LTMP") %>%
  left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>%
  filter(!is.na(mcur_90))

dat.ml.mmp <- dat.mod.bent.hc %>%
  filter(project_code == "MMP") %>%
  left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>%
  filter(!is.na(mcur_90))

# 3. Recreate acrop_comp
acrop_comp <- df.AIMS.full %>%
  filter(
    data_type == "photo-transect", domain_category == "reef",
    purpose == "COMPOSITION", variable == "HARD CORAL"
  ) %>%
  mutate(Reef_Name_Clean = gsub(" REEF(S)?| ISLAND| IS", "", toupper(domain_name))) %>%
  left_join(dat.Ref.Clean %>% select(AIMS_REEF_NAME_Clean, ReefName),
    by = c("Reef_Name_Clean" = "AIMS_REEF_NAME_Clean")
  ) %>%
  filter(!is.na(ReefName)) %>%
  group_by(ReefName, report_year, depth) %>%
  summarise(
    total_hc = sum(mean, na.rm = TRUE),
    acrop_cover = sum(mean[grepl("Acropora", reefpage_category, ignore.case = TRUE)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(prop_acropora = ifelse(total_hc > 0, acrop_cover / total_hc, 0))

# 4. Prepare LTMP clean data
dat.ml.ltmp.clean <- dat.ml.ltmp %>%
  filter(!is.na(mcur_90), !is.na(Rel.Change)) %>%
  left_join(acrop_comp %>% select(ReefName, report_year, depth, prop_acropora),
    by = c("ReefName", "report_year", "depth")
  ) %>%
  mutate(prop_acropora = replace_na(prop_acropora, median(prop_acropora, na.rm = TRUE)))

N_ltmp <- nrow(dat.ml.ltmp.clean)
dat.ml.ltmp.clean <- dat.ml.ltmp.clean %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
    Mort.prop.nudge = (Mort.prop * (N_ltmp - 1) + 0.5) / N_ltmp
  )

# Scale LTMP
scale_params_ltmp <- list(
  MaxDHW.mean = c(mean = mean(dat.ml.ltmp.clean$MaxDHW.mean, na.rm = T), sd2 = 2 * sd(dat.ml.ltmp.clean$MaxDHW.mean, na.rm = T)),
  histmDHW6 = c(mean = mean(dat.ml.ltmp.clean$histmDHW6, na.rm = T), sd2 = 2 * sd(dat.ml.ltmp.clean$histmDHW6, na.rm = T)),
  secc3m = c(mean = mean(dat.ml.ltmp.clean$secc3m, na.rm = T), sd2 = 2 * sd(dat.ml.ltmp.clean$secc3m, na.rm = T)),
  cloudp_90 = c(mean = mean(dat.ml.ltmp.clean$cloudp_90, na.rm = T), sd2 = 2 * sd(dat.ml.ltmp.clean$cloudp_90, na.rm = T)),
  winyear_sd = c(mean = mean(dat.ml.ltmp.clean$winyear_sd, na.rm = T), sd2 = 2 * sd(dat.ml.ltmp.clean$winyear_sd, na.rm = T)),
  mcur_90 = c(mean = mean(dat.ml.ltmp.clean$mcur_90, na.rm = T), sd2 = 2 * sd(dat.ml.ltmp.clean$mcur_90, na.rm = T)),
  prop_acropora = c(mean = mean(dat.ml.ltmp.clean$prop_acropora, na.rm = T), sd2 = 2 * sd(dat.ml.ltmp.clean$prop_acropora, na.rm = T))
)

s2 <- function(x) (x - mean(x, na.rm = T)) / (2 * sd(x, na.rm = T))
dat.ml.ltmp.clean <- dat.ml.ltmp.clean %>%
  mutate(
    MaxDHW_s = s2(MaxDHW.mean), histmDHW6_s = s2(histmDHW6), secc3m_s = s2(secc3m),
    cloudp_90_s = s2(cloudp_90), winyear_sd_s = s2(winyear_sd), mcur_90_s = s2(mcur_90),
    prop_acropora_s = s2(prop_acropora), obs_id = 1:n(),
    n_trials = 100L, y_binom = as.integer(round(Mort.prop.nudge * 100))
  )

# Fit and save LTMP Binomial-OLRE model (replaces corrupted file)
cat("Fitting LTMP Binomial-OLRE model...\n")
fit_brmsb_ltmp <- brm(
  y_binom | trials(n_trials) ~ MaxDHW_s * cloudp_90_s + MaxDHW_s * histmDHW6_s +
    MaxDHW_s * mcur_90_s + MaxDHW_s * prop_acropora_s +
    secc3m_s + winyear_sd_s + (1 | ReefName) + (1 | obs_id),
  data = dat.ml.ltmp.clean, family = binomial(link = "logit"),
  prior = c(
    prior(normal(-1, 1.5), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(1), class = "sd")
  ),
  chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
  control = list(adapt_delta = 0.97)
)
saveRDS(fit_brmsb_ltmp, "output/models/brms_ltmp_binomial_olre.rds")
cat("Saved LTMP Binomial-OLRE model successfully.\n")

# 5. Prepare MMP clean data
dat.ml.mmp.clean <- dat.ml.mmp %>%
  filter(!is.na(mcur_90), !is.na(Rel.Change)) %>%
  left_join(acrop_comp %>% select(ReefName, report_year, depth, prop_acropora),
    by = c("ReefName", "report_year", "depth")
  ) %>%
  mutate(prop_acropora = replace_na(prop_acropora, median(prop_acropora, na.rm = TRUE)))

N_mmp <- nrow(dat.ml.mmp.clean)
dat.ml.mmp.clean <- dat.ml.mmp.clean %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
    Mort.prop.nudge = (Mort.prop * (N_mmp - 1) + 0.5) / N_mmp
  )

dat.ml.mmp.clean <- dat.ml.mmp.clean %>%
  mutate(
    MaxDHW_s = s2(MaxDHW.mean), cloudp_90_s = s2(cloudp_90),
    secc3m_s = s2(secc3m), prop_acropora_s = s2(prop_acropora),
    obs_id = 1:n(),
    n_trials = 100L, y_binom = as.integer(round(Mort.prop.nudge * 100))
  )

# Fit and save MMP Beta model
cat("Fitting MMP Beta model...\n")
fit_brms_mmp <- brm(
  Mort.prop.nudge ~ MaxDHW_s * cloudp_90_s + MaxDHW_s * prop_acropora_s +
    secc3m_s + depth + (1 | ReefName),
  data = dat.ml.mmp.clean, family = Beta(link = "logit"),
  prior = c(
    prior(normal(-1, 1.5), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(2), class = "sd")
  ),
  chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
  control = list(adapt_delta = 0.98)
)
saveRDS(fit_brms_mmp, "output/models/brms_mmp_beta.rds")
cat("Saved MMP Beta model successfully.\n")

# Fit and save MMP Binomial-OLRE model
cat("Fitting MMP Binomial-OLRE model...\n")
fit_brmsb_mmp <- brm(
  y_binom | trials(n_trials) ~ MaxDHW_s * cloudp_90_s + MaxDHW_s * prop_acropora_s +
    secc3m_s + depth + (1 | ReefName) + (1 | obs_id),
  data = dat.ml.mmp.clean, family = binomial(link = "logit"),
  prior = c(
    prior(normal(-1, 1.5), class = "Intercept"),
    prior(normal(0, 1), class = "b"),
    prior(exponential(2), class = "sd")
  ),
  chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
  control = list(adapt_delta = 0.98)
)
saveRDS(fit_brmsb_mmp, "output/models/brms_mmp_binomial_olre.rds")
cat("Saved MMP Binomial-OLRE model successfully.\n")
