# Standalone script to fit and cache new brms models
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(brms)
  library(lubridate)
  library(sf)
})

cat("=== Initializing Data for Model Fitting ===\n")

if (!dir.exists("output/models")) dir.create("output/models", recursive = TRUE)

# Load AIMS data & DHW data
load("Data/aims_ltmp/aims_ltmp.RData")
dat.AIMSRef <- read.csv("Data/AIMS-Reef_Reference.csv") %>% rename(Reef_Name = AIMS_REEF_NAME)
Disturbances_manta <- read.csv("Data/Disturbances_manta.csv") %>%
  left_join(select(dat.AIMSRef, Reef_Name, ReefName), by = c("AIMS_REEF_NAME" = "Reef_Name")) %>%
  rename("report_year" = "REPORT_YEAR")
dat.DHW <- read.csv("Data/DHW_1985_2024_GBRReefs.csv") %>% rename("ReefName" = "LOC_NAME_S")

# Growth parameters
dat.grid <- read.csv("Data/COTSModParams.csv")
dat.gridSum <- dat.grid %>%
  mutate(
    REEF_ID = ifelse(grepl("14-116", REEF_ID, fixed = TRUE), "14-116", REEF_ID),
    REEF_NAME = ifelse(grepl("14-116", REEF_NAME, fixed = TRUE), "Lizard Island Reef (14-116)", REEF_NAME)
  ) %>%
  group_by(REEF_NAME, REEF_ID, SECTOR) %>%
  summarise(
    B0 = mean(pred.b0.mean),
    HC.asym = mean(pred.HCmax.mean),
    WQ1 = mean(Primary), WQ2 = mean(Secondary), WQ3 = mean(Tertiary),
    .groups = "drop"
  ) %>%
  mutate(WQ = WQ1 + WQ2 + WQ3)

doCoralGrowth.Simp <- function(CoralCover, B0, HC.asym, WQ) {
  WQ.mn.sd <- c(-0.68, 0.03)
  b0.wq <- B0 + WQ * WQ.mn.sd[1]
  b1.wq <- b0.wq / log(HC.asym)
  CoralCover <- log(CoralCover)
  b0.wq <- as.numeric(b0.wq)
  b1.wq <- as.numeric(b1.wq)
  CoralCover <- as.numeric(CoralCover)
  CoralCover.t1 <- (b0.wq + (1 - b1.wq) * CoralCover)
  CoralCover.t1 <- exp(CoralCover.t1)
  Change.t1 <- CoralCover.t1 - exp(CoralCover)
  return(c(exp(CoralCover), CoralCover.t1, Change.t1))
}

# Agincourt DHW fill
dat.agin <- dat.DHW %>%
  filter(grepl("Agincourt", ReefName)) %>%
  group_by(Year) %>%
  summarise(MaxDHW.mean = mean(MaxDHW.mean, na.rm = TRUE), .groups = "drop")

dat.DHW <- dat.DHW %>%
  left_join(dat.agin, by = "Year") %>%
  mutate(MaxDHW.mean = ifelse(grepl("Agincourt", ReefName) & is.na(MaxDHW.mean.x), MaxDHW.mean.y, MaxDHW.mean.x)) %>%
  select(-MaxDHW.mean.x, -MaxDHW.mean.y)

df.AIMS.Mant <- reef_manta_df %>%
  select(domain_name, project_code, report_year, date, depth, lower:median, mean, reef_zone, id, shelf, reefpage_category) %>%
  mutate(Reef_Name = toupper(domain_name))

dat.AIMSRef_key <- dat.AIMSRef %>%
  filter(!ReefName %in% c("Round-Russell Reef (17-013)", "Snake Reef (14-087)"))

df.MANT <- df.AIMS.Mant %>%
  left_join(select(dat.AIMSRef_key, AIMS_REEF_NAME_cap, ReefID, ReefName, SECT_NAME, SECTOR),
            by = join_by("Reef_Name" == "AIMS_REEF_NAME_cap")) %>%
  mutate(
    year = floor(as.numeric(date)),
    frac = as.numeric(date) - year,
    date_approx = ymd(paste0(year, "-01-01")) +
      round(frac * as.numeric(ymd(paste0(year + 1, "-01-01")) - ymd(paste0(year, "-01-01"))))
  )

dat.yrs.mant <- df.MANT %>%
  ungroup() %>%
  filter(!is.na(ReefName)) %>%
  mutate(
    Date = date_approx,
    Year = year(Date), Month = month(Date),
    DHWYear = ifelse(Month >= 7, Year, Year - 1)
  ) %>%
  left_join(select(dat.DHW, ReefName, Year, MaxDHW.mean), by = c("ReefName", "DHWYear" = "Year"))

df.MANT <- df.MANT %>%
  arrange(ReefName, Reef_Name, reefpage_category, depth, report_year) %>%
  group_by(ReefName, Reef_Name, depth, reefpage_category) %>%
  mutate(
    Change = mean - lag(mean, n = 1),
    Rel.Change = Change / lag(mean, n = 1),
    Lag = report_year - lag(report_year, n = 1)
  ) %>%
  left_join(dat.yrs.mant)

df.MANT <- df.MANT %>%
  mutate(Change.Cor = Change, Prev.HC = mean - Change.Cor) %>%
  filter(ReefName %in% dat.gridSum$REEF_NAME, !is.na(Lag), Reef_Name != "MILLN REEF") %>%
  rowwise() %>%
  mutate(
    Pred.Grth = doCoralGrowth.Simp(
      CoralCover = Prev.HC * 100,
      B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"],
      WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"],
      HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"]
    )[3] / 100
  ) %>%
  mutate(
    Pred.Grth2 = ifelse(Lag == 2,
                        Pred.Grth + doCoralGrowth.Simp(
                          CoralCover = (Prev.HC + Pred.Grth) * 100,
                          B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"],
                          WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"],
                          HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"]
                        )[3] / 100,
                        Pred.Grth),
    Abs.Change = Change.Cor - Pred.Grth2,
    Rel.Change = Abs.Change / (Prev.HC + Pred.Grth2)
  )

dat.mod.mant <- df.MANT %>%
  filter(Lag < 3, project_code != "MMP", abs(Change) > 0.01) %>%
  left_join(select(Disturbances_manta, ReefName, report_year, DISTURBANCE_TYPE),
            by = c("ReefName", "report_year")) %>%
  filter(
    is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n"),
    report_year %in% c(1999, 2003, 2017, 2018, 2021, 2023, 2024, 2025)
  ) %>%
  mutate(
    Era = case_when(
      DHWYear %in% c(1998, 2002) ~ "1998/2002",
      DHWYear %in% c(2016, 2017) ~ "2016/2017",
      DHWYear %in% c(2020, 2022) ~ "2020/2022",
      DHWYear == 2024 ~ "2024",
      TRUE ~ NA_character_
    )
  )

# Add Cheung predictors for N=110 subset
load("data/Cheungetal2025/01_sstvar_blchrf.RData")
predictor_vars <- c("mcur_90", "cloudp_90", "secc3m", "cbclus2", "histmDHW6", "yrsince6", "histmDHW4", "yrsince4", "winyear_sd", "winyear_mean")
cheung_predictors <- sstvar_blch2 %>%
  select(LABEL_id, year, lat = Y, lon = X, Sector, all_of(predictor_vars)) %>%
  group_by(LABEL_id, year) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(year = as.numeric(year))

dat.ml.mant <- dat.mod.mant %>%
  ungroup() %>%
  filter(!is.na(ReefID), !is.na(DHWYear)) %>%
  mutate(DHWYear = as.numeric(DHWYear)) %>%
  left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>%
  filter(!is.na(mcur_90))

dat.ml.mant.clean <- dat.ml.mant
N_obs <- nrow(dat.ml.mant.clean)

dat.ml.mant.clean <- dat.ml.mant.clean %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
    Mort.prop.nudge = (Mort.prop * (N_obs - 1) + 0.5) / N_obs
  )

scale_2sd <- function(x) {
  (x - mean(x, na.rm = TRUE)) / (2 * sd(x, na.rm = TRUE))
}

dat.ml.mant.clean <- dat.ml.mant.clean %>%
  mutate(
    MaxDHW_s     = scale_2sd(MaxDHW.mean),
    winyear_sd_s = scale_2sd(winyear_sd),
    histmDHW6_s  = scale_2sd(histmDHW6),
    mcur_90_s    = scale_2sd(mcur_90),
    cloudp_90_s  = scale_2sd(cloudp_90),
    secc3m_s     = scale_2sd(secc3m)
  )

dat.ml.mant.clean$n_trials <- 100L
dat.ml.mant.clean$y_binom <- as.integer(round(dat.ml.mant.clean$Mort.prop.nudge * 100))

cat("N=110 subset loaded with", nrow(dat.ml.mant.clean), "rows.\n")
cat("N=465 full dataset loaded with", nrow(dat.mod.mant), "rows.\n")

# ─────────────────────────────────────────────────────────────
# MODEL 1: Zero-Inflated Beta Model (fit3c_brms)
# ─────────────────────────────────────────────────────────────
cache_fit3c <- "output/models/brms_fit3c_zero_inflated_beta.rds"
if (file.exists(cache_fit3c)) {
  cat("fit3c (Zero-Inflated Beta) already exists.\n")
} else {
  cat("Fitting fit3c (Zero-Inflated Beta)... \n")
  fit3c_brms <- brm(
    formula = bf(
      Mort.prop ~ MaxDHW_s * secc3m_s + MaxDHW_s * cloudp_90_s +
        MaxDHW_s * histmDHW6_s + MaxDHW_s * mcur_90_s +
        winyear_sd_s + (1 | ReefName),
      zi ~ MaxDHW_s
    ),
    data = dat.ml.mant.clean,
    family = zero_inflated_beta(link = "logit", link_zi = "logit"),
    prior = c(
      prior(normal(-1, 1.5), class = "Intercept"),
      prior(normal(0.3, 0.3), class = "b", coef = "MaxDHW_s"),
      prior(normal(0, 1), class = "b"),
      prior(exponential(1), class = "sd")
    ),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.97)
  )
  saveRDS(fit3c_brms, cache_fit3c)
  cat("Saved fit3c to", cache_fit3c, "\n")
}

# ─────────────────────────────────────────────────────────────
# MODEL 2: Standard Binomial Model (No OLRE) (fit3d_brms)
# ─────────────────────────────────────────────────────────────
cache_fit3d <- "output/models/brms_fit3d_binomial_standard.rds"
if (file.exists(cache_fit3d)) {
  cat("fit3d (Standard Binomial) already exists.\n")
} else {
  cat("Fitting fit3d (Standard Binomial without OLRE)... \n")
  fit3d_brms <- brm(
    formula = y_binom | trials(n_trials) ~ MaxDHW_s * secc3m_s + MaxDHW_s * cloudp_90_s +
      MaxDHW_s * histmDHW6_s + MaxDHW_s * mcur_90_s +
      winyear_sd_s + (1 | ReefName),
    data = dat.ml.mant.clean,
    family = binomial(link = "logit"),
    prior = c(
      prior(normal(-1, 1.5), class = "Intercept"),
      prior(normal(0.3, 0.3), class = "b", coef = "MaxDHW_s"),
      prior(normal(0, 1), class = "b"),
      prior(exponential(1), class = "sd")
    ),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.97)
  )
  saveRDS(fit3d_brms, cache_fit3d)
  cat("Saved fit3d to", cache_fit3d, "\n")
}

# ─────────────────────────────────────────────────────────────
# MODEL 3: Full-Dataset Zero-Inflated Beta (N=465) (fit_full_brms)
# ─────────────────────────────────────────────────────────────
cache_fit_full <- "output/models/brms_fit_full_manta_n465.rds"
if (file.exists(cache_fit_full)) {
  cat("fit_full (N=465 Zero-Inflated Beta) already exists.\n")
} else {
  cat("Fitting fit_full (Full N=465 Dataset Zero-Inflated Beta)... \n")
  dat.mod.mant.prep <- dat.mod.mant %>%
    ungroup() %>%
    filter(is.finite(MaxDHW.mean), is.finite(Rel.Change)) %>%
    mutate(
      Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
      MaxDHW_s = scale_2sd(MaxDHW.mean)
    )

  fit_full_brms <- brm(
    formula = bf(
      Mort.prop ~ MaxDHW_s + Era + (1 | ReefName),
      zi ~ MaxDHW_s
    ),
    data = dat.mod.mant.prep,
    family = zero_inflated_beta(link = "logit", link_zi = "logit"),
    prior = c(
      prior(normal(-1, 1.5), class = "Intercept"),
      prior(normal(0.3, 0.3), class = "b", coef = "MaxDHW_s"),
      prior(normal(0, 1), class = "b"),
      prior(exponential(1), class = "sd")
    ),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.97)
  )
  saveRDS(fit_full_brms, cache_fit_full)
  cat("Saved fit_full to", cache_fit_full, "\n")
}

cat("=== Model Fitting Script Finished ===\n")
