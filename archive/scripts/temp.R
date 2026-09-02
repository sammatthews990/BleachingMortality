## ----setup--------------------------------------------------------------------
#| label: setup

# ── Core libraries ──
library(tidyr)
library(dplyr)
library(ggplot2)
library(DT)
library(arrow)
library(stringr)
library(janitor)
library(lubridate)
library(readr)
library(sf)
library(lwgeom)
library(scales)
library(patchwork)

# ── Modelling libraries ──
library(mgcv)
library(betareg)
library(MASS)
library(lme4)
library(sandwich)
library(lmtest)

# ── Shared colour palettes ──
dist_colours <- c(
  "b" = "darkorange", "B" = "darkorange",
  "c" = "maroon4",    "C" = "maroon4",
  "d" = "darkolivegreen3", "D" = "darkolivegreen3",
  "f" = "burlywood4", "F" = "burlywood4",
  "s" = "dodgerblue2", "S" = "dodgerblue2",
  "u" = "grey",       "U" = "grey"
)

cols_exposure <- c("First" = "#E69F00", "Repeat" = "#0072B2")
cols_outbreak <- c("COTS Outbreak" = "#D55E00", "No Outbreak" = "#009E73")

secorder <- c("CG", "PC", "CL", "CA", "IN", "TO", "CU", "WH", "PO", "SW", "CB")
secorderlabs <- c("Cape Grenville", "Princess Charlotte Bay", "Cooktown / Lizard Island",
                  "Cairns", "Innisfail", "Townsville", "Cape Upstart", "Whitsunday",
                  "Pompey", "Swain", "Capricorn Bunker")


## ----helper-functions---------------------------------------------------------
#| label: helper-functions

# ── Gompertz coral growth function (simplified, single time-step) ──
# Returns c(current_cover, next_cover, change)
doCoralGrowth.Simp <- function(CoralCover, B0, HC.asym, WQ) {
  WQ.mn.sd <- c(-0.68, 0.03)
  b0.wq <- B0 + WQ * WQ.mn.sd
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

# ── Multi-step growth function (n years from starting cover) ──
doCoralGrowth <- function(n, CoralCover, B0, HC.asym, WQ) {
  WQ.mn.sd <- c(-0.68, 0.03)
  b0.wq <- B0 + WQ * WQ.mn.sd
  b1.wq <- b0.wq / log(HC.asym)
  CoralCover <- log(CoralCover)
  b0.wq <- as.numeric(b0.wq)
  b1.wq <- as.numeric(b1.wq)
  CoralCover <- as.numeric(CoralCover)
  df <- matrix(nrow = n + 1, ncol = 3)
  df[1, 1] <- 0
  df[1, 2] <- exp(CoralCover)
  for (i in 1:n) {
    CoralCover.t1 <- (b0.wq + (1 - b1.wq) * CoralCover)
    df[i + 1, 1] <- i
    df[i + 1, 2] <- exp(CoralCover.t1)
    df[i + 1, 3] <- exp(CoralCover.t1) - exp(CoralCover)
    CoralCover <- CoralCover.t1
  }
  colnames(df) <- c("Year", "CoralCover", "CoralGrowth")
  return(df)
}

# ── Inverse logit ──
inv_logit <- function(x) 1 / (1 + exp(-x))


## ----growth-params------------------------------------------------------------
#| label: growth-params
#| cache: true

dat.grid <- read.csv("Data/COTSModParams.csv")

dat.gridSum <- dat.grid %>%
  mutate(
    REEF_ID = ifelse(grepl("14-116", REEF_ID, fixed = TRUE), "14-116", REEF_ID),
    REEF_NAME = ifelse(grepl("14-116", REEF_NAME, fixed = TRUE),
                       "Lizard Island Reef (14-116)", REEF_NAME)
  ) %>%
  group_by(REEF_NAME, REEF_ID, SECTOR) %>%
  summarise(
    B0 = mean(pred.b0.mean),
    HC.asym = mean(pred.HCmax.mean),
    WQ1 = mean(Primary),
    WQ2 = mean(Secondary),
    WQ3 = mean(Tertiary),
    .groups = "drop"
  ) %>%
  mutate(
    WQ = WQ1 + WQ2 + WQ3,
    HC90 = 0.8 * HC.asym
  )


## ----growth-curves------------------------------------------------------------
#| label: growth-curves
#| cache: true
#| fig-cap: "Modelled coral recovery trajectories from 0.1% cover, by sector. Ribbons show inter-reef variability."

# Simulate 30-year recovery from 0.1% cover for all reefs
CC.tmp <- data.frame()
for (i in seq_len(nrow(dat.gridSum))) {
  CCtmp2 <- as.data.frame(
    doCoralGrowth(30, 0.1, dat.gridSum$B0[i], dat.gridSum$HC.asym[i], dat.gridSum$WQ[i])
  )
  CCtmp2$ReefName <- dat.gridSum$REEF_NAME[i]
  CCtmp2$REEF_NAME <- dat.gridSum$REEF_NAME[i]
  CCtmp2$ReefID <- dat.gridSum$REEF_ID[i]
  CCtmp2$Sector <- dat.gridSum$SECTOR[i]
  CC.tmp <- bind_rows(CC.tmp, CCtmp2)
}

CC.tmp <- CC.tmp %>% left_join(select(dat.gridSum, REEF_NAME, HC90), by = "REEF_NAME")

library(ggdist)
ggplot(CC.tmp, aes(x = Year, y = CoralCover)) +
  stat_lineribbon(alpha = 0.3) +
  facet_wrap(~Sector) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "Years since disturbance", y = "Coral cover (%)")


## ----aims-timeseries----------------------------------------------------------
#| label: aims-timeseries
#| cache: false

dat.AIMSRef <- read.csv("Data/AIMS-Reef_Reference.csv") %>%
  rename(Reef_Name = AIMS_REEF_NAME)

dat.AIMS <- read.csv("Data/AIMS_ModelledCoralCover_Reef.csv") %>%
  left_join(dat.AIMSRef, by = "Reef_Name") %>%
  mutate(
    REEF_ID = ifelse(grepl("14-116", ReefID, fixed = TRUE), "14-116", ReefID),
    ReefName = ifelse(grepl("14-116", ReefName, fixed = TRUE),
                      "Lizard Island Reef (14-116)", ReefName),
    domain_name = factor(SECT_NAME),
    date = as.Date(date)
  ) %>%
  rename(
    median_HC = median, lower_HC_CI = lower,
    upper_HC_CI = upper, mean_HC = mean
  ) %>%
  mutate(across(c(median_HC, lower_HC_CI, upper_HC_CI, mean_HC), ~ . * 100)) %>%
  relocate(ReefName, .after = Reef_Name)

datatable(
  dat.AIMS %>% mutate(across(where(is.numeric), ~ round(., 4))),
  style = "default", caption = "AIMS modelled coral cover time-series"
)


## ----aims-plots---------------------------------------------------------------
#| label: aims-plots
#| fig-cap: "Modelled coral cover trajectories for selected reefs with credible intervals."

ggplot(dat.AIMS %>% filter(ReefID %in% c("14-116", "17-034", "18-075", "23-082a")),
       aes(x = report_year, y = median_HC, colour = ReefName, fill = ReefName)) +
  geom_line(colour = "black") +
  geom_pointrange(aes(ymin = lower_HC_CI, ymax = upper_HC_CI)) +
  geom_ribbon(aes(ymin = lower_HC_CI, ymax = upper_HC_CI), alpha = 0.3, colour = NA) +
  facet_wrap(~ReefName, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Year", y = "Hard coral cover (%)")


## ----disturbance-records------------------------------------------------------
#| label: disturbance-records

dat.DIST <- read.csv("Data/Disturbances_manta.csv")

df.DIST <- dat.DIST %>%
  left_join(dat.AIMSRef, by = "FULLREEF_ID") %>%
  mutate(DISTURBANCE_TYPE = replace(DISTURBANCE_TYPE, DISTURBANCE_TYPE == "m", "u")) %>%
  filter(!SECT_NAME %in% c("Outside AIMS field sectors", "PO", "CG", "PB")) %>%
  filter(DISTURBANCE_TYPE != "n") %>%
  mutate(
    DISTURBANCE_TYPE = gsub(".*\\b[b]\\b.*", "B", DISTURBANCE_TYPE),
    DISTURBANCE_TYPE = gsub(".*\\b[c]\\b.*", "C", DISTURBANCE_TYPE),
    DISTURBANCE_TYPE = gsub(".*\\b[d]\\b.*", "D", DISTURBANCE_TYPE),
    DISTURBANCE_TYPE = gsub(".*\\b[f]\\b.*", "F", DISTURBANCE_TYPE),
    DISTURBANCE_TYPE = gsub(".*\\b[s]\\b.*", "S", DISTURBANCE_TYPE),
    DISTURBANCE_TYPE = gsub(".*\\b[u]\\b.*", "U", DISTURBANCE_TYPE),
    DISTURBANCE_TYPE = as.factor(DISTURBANCE_TYPE)
  ) %>%
  group_by(SECT_NAME, REPORT_YEAR, DISTURBANCE_TYPE) %>%
  distinct(AIMS_REEF_NAME, REEF_NAME, REPORT_YEAR, DISTURBANCE_TYPE) %>%
  mutate(SECTOR = factor(SECT_NAME, levels = secorderlabs, labels = secorderlabs))


## ----disturbance-benthic------------------------------------------------------
#| label: disturbance-benthic

dat.mort.att <- read_csv("Data/240522 dist.lookup_AIMSLTMP.csv") %>%
  mutate(
    Decade = plyr::round_any(year, 10, floor),
    AllTime = "AllTime"
  )

dat.mort.att_filtered <- dat.mort.att %>%
  filter(change < 0) %>%
  group_by(Decade, DISTURBANCE_TYPE) %>%
  summarize(mean_change = mean(change, na.rm = TRUE), .groups = "drop")


## ----disturbance-plot---------------------------------------------------------
#| label: disturbance-attribution-plot
#| fig-cap: "Proportional attribution of coral cover decline by disturbance type across decades."

ggplot(dat.mort.att %>% filter(change < 0),
       aes(x = Decade, y = change, colour = DISTURBANCE_TYPE, fill = DISTURBANCE_TYPE)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_bw() +
  scale_fill_manual(values = dist_colours) +
  scale_colour_manual(values = dist_colours) +
  scale_y_continuous(labels = percent) +
  facet_wrap(~A_SECTOR) +
  labs(x = "Decade", y = "Proportion of coral decline", fill = "Disturbance", colour = "Disturbance")


## ----dhw-mortality-comparison-------------------------------------------------
#| label: dhw-mortality-comparison

dat.bleach <- read.csv("Data/DHW_vs_Coral_Cover_Dataset.csv") %>%
  mutate(Mort.bin = ifelse(Change >= 0, 0.001, -Change / 100))

mod.beta <- betareg(Mort.bin ~ DHW, data = dat.bleach)
mod.bin <- glm(Mort.bin ~ DHW, data = dat.bleach, family = binomial(link = "logit"))

# Predictions
newdata <- with(dat.bleach, data.frame(DHW = seq(0, 16, len = 100)))
fit_bin <- predict(mod.bin, newdata = newdata, type = "link", se = TRUE)
q <- qt(0.975, df = df.residual(mod.bin))
newdata <- cbind(newdata,
  fit = binomial()$linkinv(fit_bin$fit),
  lower = binomial()$linkinv(fit_bin$fit - q * fit_bin$se.fit),
  upper = binomial()$linkinv(fit_bin$fit + q * fit_bin$se.fit)
)

newdata2 <- with(dat.bleach, data.frame(DHW = seq(0, 16, len = 100)))
fit_beta <- predict(mod.beta, newdata = newdata2, type = "link")
newdata2 <- cbind(newdata2, fit = binomial()$linkinv(fit_beta))


## ----dhw-mortality-plot-------------------------------------------------------
#| label: dhw-mortality-plot
#| fig-cap: "Comparison of binomial (black) and beta regression (red) fits to the DHW–mortality relationship from Hughes et al. (2018) data."

ggplot(data = newdata, aes(y = fit, x = DHW)) +
  geom_point(data = dat.bleach, aes(y = Mort.bin)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) +
  geom_line() +
  geom_line(data = newdata2, aes(y = fit, x = DHW), colour = "red") +
  theme_classic() +
  labs(x = "DHW (°C-weeks)", y = "Mortality proportion",
       caption = "Black = binomial GLM, Red = beta regression")


## ----ltmp-mortality-----------------------------------------------------------
#| label: ltmp-mortality

Disturbances_manta <- read.csv("Data/Disturbances_manta.csv") %>%
  left_join(select(dat.AIMSRef, Reef_Name, ReefName), by = c("AIMS_REEF_NAME" = "Reef_Name")) %>%
  rename("report_year" = "REPORT_YEAR")

dat.mort <- read.csv("Data/DisturbanceMortality.csv") %>%
  filter(!is.na(Lag), !Lag > 2) %>%
  left_join(select(dat.AIMSRef, Reef_Name, ReefName), by = "Reef_Name") %>%
  select(-DISTURBANCE_TYPE) %>%
  distinct() |>
  left_join(select(Disturbances_manta, ReefName, report_year, DISTURBANCE_TYPE),
            by = c("ReefName", "report_year"))

dat.mort.bl <- dat.mort %>%
  filter(
    !DISTURBANCE_TYPE %in% c("c", "d", "s", "f"),
    report_year %in% c(1999, 2003, 2017, 2018, 2021, 2024)
  )

# DHW data
dat.DHW <- read.csv("Data/DHW_1985_2024_GBRReefs.csv") %>%
  rename("ReefName" = "LOC_NAME_S")

# Fill Agincourt reefs (missing individual DHW)
dat.agin <- dat.DHW |>
  filter(grepl("Agincourt", ReefName)) |>
  group_by(Year) |>
  summarise(MaxDHW.mean = mean(MaxDHW.mean, na.rm = TRUE), .groups = "drop")

dat.DHW <- dat.DHW |>
  left_join(dat.agin, by = "Year") |>
  mutate(MaxDHW.mean = ifelse(grepl("Agincourt", ReefName) & is.na(MaxDHW.mean.x),
                              MaxDHW.mean.y, MaxDHW.mean.x)) |>
  select(-MaxDHW.mean.x, -MaxDHW.mean.y)

# Max DHW across current and previous year
dat.yrs <- dat.mort.bl %>%
  select(ReefName, report_year) |>
  mutate(Prev = report_year - 1) |>
  left_join(select(dat.DHW, ReefName, Year, MaxDHW.mean),
            by = c("ReefName", "report_year" = "Year")) |>
  left_join(select(dat.DHW, ReefName, Year, MaxDHW.mean),
            by = c("ReefName", "Prev" = "Year")) |>
  rowwise() |>
  mutate(MaxDHW = max(MaxDHW.mean.x, MaxDHW.mean.y)) |>
  select(-MaxDHW.mean.x, -MaxDHW.mean.y)


## ----growth-adjusted----------------------------------------------------------
#| label: growth-adjusted

dat.grth <- dat.mort.bl |>
  mutate(Prev.HC = mean_HC - Abs.Change.Cor) |>
  filter(ReefName %in% dat.gridSum$REEF_NAME) |>
  rowwise() %>%
  mutate(
    Pred.Grth = doCoralGrowth.Simp(
      CoralCover = Prev.HC,
      B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"],
      WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"],
      HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"]
    )[3]
  ) |>
  mutate(
    Abs.Change.Cor = Abs.Change.Cor - Pred.Grth,
    Rel.Change.Cor = Abs.Change.Cor / (Prev.HC + Pred.Grth)
  )

dat.mod <- dat.grth |>
  left_join(dat.yrs, by = c("ReefName", "report_year")) |>
  mutate(
    Era = case_when(
      report_year %in% c(1999, 2003) ~ "1998/2002",
      report_year %in% c(2017, 2018) ~ "2016/2017",
      report_year %in% c(2021, 2023) ~ "2020/2022",
      TRUE ~ NA_character_
    )
  )

# Combine with Hughes data
dat.mod.Terry <- dat.mod |>
  bind_rows(
    rename(dat.bleach, "Rel.Change.Cor" = "Change", "MaxDHW" = "DHW") |>
      mutate(Era = "Hughes et al 2018", Rel.Change.Cor = Rel.Change.Cor / 100)
  ) |>
  mutate(
    Region = ifelse(SECTOR_NAME %in% c("CG", "PC", "CL"), "Northern",
             ifelse(SECTOR_NAME %in% c("CA", "IN", "TO", "WH"), "Central",
             ifelse(SECTOR_NAME %in% c("PO", "SW", "CB"), "Southern", "NA")))
  )


## ----growth-adjusted-plot-----------------------------------------------------
#| label: growth-adjusted-era-plot
#| fig-cap: "Growth-adjusted relative coral cover change vs DHW, stratified by bleaching era. Includes Hughes et al. (2018) data for comparison."

ggplot(dat.mod.Terry, aes(x = MaxDHW, y = Rel.Change.Cor, fill = Era, colour = Era)) +
  geom_point(shape = 17, alpha = 0.5) +
  scale_fill_viridis_d(end = 0.98, option = "D", begin = 0.1, direction = -1) +
  scale_color_viridis_d(end = 0.98, option = "D", begin = 0.1, direction = -1) +
  coord_cartesian(xlim = c(0, 12)) +
  geom_smooth(method = "gam") +
  scale_y_continuous(labels = percent) +
  theme_bw() +
  labs(x = "DHW (°C-weeks)", y = "Relative change (growth-adjusted)")


## ----manta-tow-data-----------------------------------------------------------
#| label: manta-tow-data
#| cache: true

load("Data/aims_ltmp/aims_ltmp.RData")

df.AIMS.Mant <- reef_manta_df |>
  select(domain_name, project_code, report_year, date, depth,
         lower:median, mean, reef_zone, id, shelf, reefpage_category) |>
  mutate(Reef_Name = toupper(domain_name))

library(fuzzyjoin)

dat.AIMSRef_key <- dat.AIMSRef %>%
  filter(!ReefName %in% c("Round-Russell Reef (17-013)", "Snake Reef (14-087)")) %>%
  mutate(AIMS_REEF_NAME_case = str_to_title(AIMS_REEF_NAME_cap))

df.MANT <- df.AIMS.Mant |>
  left_join(select(dat.AIMSRef_key, AIMS_REEF_NAME_cap, ReefID, ReefName, SECT_NAME, SECTOR),
            by = join_by("Reef_Name" == "AIMS_REEF_NAME_cap")) %>%
  mutate(
    year = floor(as.numeric(date)),
    frac = as.numeric(date) - year,
    date_approx = ymd(paste0(year, "-01-01")) +
      round(frac * as.numeric(ymd(paste0(year + 1, "-01-01")) - ymd(paste0(year, "-01-01"))))
  ) %>%
  select(-year, -frac)

# DHW year assignment
dat.yrs.mant <- df.MANT |> ungroup() |>
  filter(!is.na(ReefName)) |>
  mutate(
    Date = date_approx, Prev = report_year - 1,
    Year = year(Date), Month = month(Date),
    DHWYear = ifelse(Month >= 7, Year, Year - 1)
  ) |>
  left_join(dplyr::select(dat.DHW, ReefName, Year, MaxDHW.mean),
            by = c("ReefName", "DHWYear" = "Year"))

# Calculate change per reef
df.MANT <- df.MANT %>%
  arrange(ReefName, Reef_Name, reefpage_category, depth, report_year) |>
  group_by(ReefName, Reef_Name, depth, reefpage_category) |>
  mutate(
    Change = mean - lag(mean, n = 1),
    Rel.Change = Change / lag(mean, n = 1),
    Lag = report_year - lag(report_year, n = 1)
  ) |>
  relocate(c("Change", "Lag", "Rel.Change"), .before = median) |>
  left_join(dat.yrs.mant)

# Growth correction
df.MANT <- df.MANT |>
  mutate(Change.Cor = Change, Prev.HC = mean - Change.Cor) |>
  filter(ReefName %in% dat.gridSum$REEF_NAME, !is.na(Lag),
         Reef_Name != "MILLN REEF") |>
  rowwise() |>
  mutate(
    Pred.Grth = doCoralGrowth.Simp(
      CoralCover = Prev.HC * 100,
      B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"],
      WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"],
      HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"]
    )[3] / 100
  ) |>
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


## ----manta-model-data---------------------------------------------------------
#| label: manta-model-data

# Load spatial regions
gbr <- st_read("Data/GBR_AIMS") |> st_make_valid()
data_bucket <- s3_bucket("s3://gbr-dms-data-public/gbrmpa-complete-gbr-features/data.parquet")
reefs <- open_dataset(data_bucket) %>%
  filter(FEAT_NAME == "Reef", DATASET == "GBR Features") |>
  distinct(UNIQUE_ID, GBR_NAME, LOC_NAME_S, geometry, FEAT_NAME) %>%
  collect() %>% st_as_sf(crs = 4283) |> st_make_valid()
reefs_centroids <- st_centroid(reefs)
reefs_join <- st_join(reefs_centroids, gbr)

dat.mod.mant <- df.MANT %>%
  filter(Lag < 3, project_code != "MMP", abs(Change) > 0.01) %>%
  left_join(select(Disturbances_manta, ReefName, report_year, DISTURBANCE_TYPE),
            by = c("ReefName", "report_year")) %>%
  filter(
    !DISTURBANCE_TYPE %in% c("c", "d", "s", "f"),
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
  ) |>
  left_join(select(reefs_join, LOC_NAME_S, Region), by = c("ReefName" = "LOC_NAME_S"))


## ----manta-plots--------------------------------------------------------------
#| label: manta-plots
#| fig-cap: "Manta tow DHW vs relative change by era and region."

ggplot(dat.mod.mant,
       aes(x = MaxDHW.mean, y = Rel.Change, fill = Era, colour = Era)) +
  geom_point(alpha = 0.5) +
  scale_fill_viridis_d(end = 0.8) +
  scale_color_viridis_d(end = 0.8) +
  coord_cartesian(xlim = c(0, 12), ylim = c(-1, 1)) +
  scale_y_continuous(labels = scales::percent) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4)) +
  theme_bw() +
  labs(x = "DHW (°C-weeks)", y = "Relative change (growth-adjusted)")


## ----lm-vs-gam----------------------------------------------------------------
#| label: lm-vs-gam
#| fig-cap: "LM vs GAM fits stratified by bleaching era."

dat.mod.mantTerry <- dat.mod.mant |>
  bind_rows(rename(dat.bleach, "Rel.Change" = "Change", "MaxDHW.mean" = "DHW") |>
              mutate(Era = "Hughes et al 2018", Rel.Change = Rel.Change / 100))

dat_plot <- dat.mod.mantTerry %>%
  filter(is.finite(MaxDHW.mean), is.finite(Rel.Change)) %>%
  mutate(Region = case_when(
    Era == "Hughes et al 2018" ~ "Hughes et al 2018",
    TRUE ~ Region
  )) %>%
  mutate(Era = factor(Era), Region = factor(Region)) %>%
  filter(!is.na(Era), !is.na(Region))

# Prediction function for LM + GAM
fit_predict <- function(dat_sub, grid_sub, k = 4) {
  m_lm <- lm(Rel.Change ~ MaxDHW.mean, data = dat_sub)
  pr_lm <- predict(m_lm, newdata = grid_sub, se.fit = TRUE)
  out_lm <- grid_sub %>%
    mutate(model = "LM", fit = pr_lm$fit,
           lwr = pr_lm$fit - 1.96 * pr_lm$se.fit,
           upr = pr_lm$fit + 1.96 * pr_lm$se.fit)

  m_gam <- mgcv::gam(Rel.Change ~ s(MaxDHW.mean, k = k), data = dat_sub, method = "REML")
  pr_gam <- predict(m_gam, newdata = grid_sub, se.fit = TRUE)
  out_gam <- grid_sub %>%
    mutate(model = paste0("GAM (k=", k, ")"), fit = pr_gam$fit,
           lwr = pr_gam$fit - 1.96 * pr_gam$se.fit,
           upr = pr_gam$fit + 1.96 * pr_gam$se.fit)

  bind_rows(out_lm, out_gam)
}

grid <- dat_plot %>%
  distinct(Era) %>%
  tidyr::expand_grid(MaxDHW.mean = seq(0, 20, by = 0.05))

pred.era <- dat_plot %>%
  group_by(Era) %>%
  group_modify(~{
    grid_sub <- grid %>% filter(Era == .y$Era) %>% dplyr::select(MaxDHW.mean)
    fit_predict(.x, grid_sub, k = 3)
  }) %>% ungroup()

x_max <- max(dat_plot$MaxDHW.mean, na.rm = TRUE)

Fig3.Eras <- ggplot(dat_plot, aes(MaxDHW.mean, Rel.Change, colour = Era, fill = Era)) +
  annotate("rect", xmin = x_max, xmax = 20, ymin = -Inf, ymax = Inf,
           fill = "grey80", alpha = 0.5) +
  geom_vline(xintercept = x_max, linetype = "dotted", linewidth = 0.8) +
  geom_point(shape = 17, alpha = 0.5) +
  geom_ribbon(data = pred.era,
              aes(x = MaxDHW.mean, y = fit, ymin = lwr, ymax = upr, group = Era),
              alpha = 0.12, colour = NA, show.legend = FALSE) +
  geom_line(data = pred.era, aes(x = MaxDHW.mean, y = fit), linewidth = 1) +
  scale_color_viridis_d(option = "A", begin = 0.2, end = 0.9) +
  scale_fill_viridis_d(option = "A", begin = 0.2, end = 0.9) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(xlim = c(0, 20), ylim = c(-1, 1)) +
  theme_classic() +
  facet_grid(~model) +
  labs(x = "DHW (°C-weeks)", y = "Relative change (growth-adjusted)",
       colour = NULL, fill = NULL)

Fig3.Eras
ggsave("output/plots/Fig3_Eras.png", Fig3.Eras, width = 8, height = 5, dpi = 300)


## ----benthic-data-------------------------------------------------------------
#| label: benthic-data
#| cache: true

df.AIMS.Bent <- df.AIMS.full |> 
  mutate(across(where(is.character), as.factor)) |> 
  filter(data_type == "photo-transect", domain_category == "reef",
         project_code %in% c("LTMP", "MMP"), purpose == "GROUP_LEVEL", variable == "HARD CORAL") |> 
  select(domain_name, project_code, date, depth, lower:median, reef_zone:fid)

df.BENT.HC <- df.AIMS.Bent |> 
  rename("Reef_Name" = "domain_name") |>
  left_join(select(dat.AIMSRef, AIMS_REEF_NAME_cap, ReefID, ReefName, SECT_NAME, SECTOR), 
            by = c("Reef_Name" = "AIMS_REEF_NAME_cap"))

# Find max DHW for each reef year
dat.yrs.bent.hc <- df.BENT.HC |> ungroup() |> 
  distinct(ReefName, date, report_year) |>
  filter(!is.na(ReefName)) |> 
  mutate(
    Date = as.Date(as.character(date)), 
    Prev = report_year - 1,
    Year = year(date), Month = month(date),
    DHWYear = ifelse(Month >= 7, Year, Year - 1)
  ) |> 
  left_join(dplyr::select(dat.DHW, ReefName, Year, MaxDHW.mean), 
            by = c("ReefName", "DHWYear" = "Year"))

# Calculate % decline for each reef
df.BENT.HC <- df.BENT.HC %>% 
  arrange(ReefName, Reef_Name, reefpage_category, depth, report_year) |>
  group_by(ReefName, Reef_Name, depth, reefpage_category) |> 
  mutate(
    Change = mean - lag(mean, n = 1),
    Rel.Change = Change / lag(mean, n = 1),
    Lag = report_year - lag(report_year, n = 1)
  ) |> 
  relocate(c("Change", "Lag", "Rel.Change"), .before = median) |> 
  left_join(dat.yrs.bent.hc, by = c("ReefName", "date", "report_year"))

# Apply growth correction
df.BENT.HC <- df.BENT.HC |> 
  mutate(Change.Cor = Change, Prev.HC = mean - Change.Cor) |> 
  filter(ReefName %in% dat.gridSum$REEF_NAME, !is.na(Lag),
         Reef_Name != "MILLN REEF") |> 
  rowwise() |> 
  mutate(
    Pred.Grth = doCoralGrowth.Simp(CoralCover = Prev.HC * 100,
      B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"],
      WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"], 
      HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"])[3] / 100
  ) |> 
  mutate(
    Pred.Grth2 = ifelse(Lag == 2, 
      Pred.Grth + doCoralGrowth.Simp(CoralCover = (Prev.HC + Pred.Grth) * 100,
        B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"],
        WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"], 
        HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"])[3] / 100,
      Pred.Grth),
    Abs.Change = Change.Cor - Pred.Grth2,
    Rel.Change = Abs.Change / (Prev.HC + Pred.Grth2)
  )

dat.mod.bent.hc <- df.BENT.HC %>%
  filter(Lag < 3) %>% 
  left_join(select(Disturbances_manta, ReefName, report_year, DISTURBANCE_TYPE), 
            by = c("ReefName", "report_year")) %>%
  filter(
    !DISTURBANCE_TYPE %in% c("c", "d", "s", "f"),
    report_year %in% c(1999, 2003, 2017, 2018, 2021, 2023, 2024)
  ) %>%
  mutate(
    Era = case_when(
      DHWYear %in% c(1998, 2002) ~ "1998/2002",
      DHWYear %in% c(2016, 2017) ~ "2016/2017",
      DHWYear %in% c(2020, 2022) ~ "2020/2022",
      DHWYear == 2024 ~ "2024",
      TRUE ~ NA_character_
    )
  ) |> 
  left_join(select(reefs_join, LOC_NAME_S, Region), by = c("ReefName" = "LOC_NAME_S"))


## ----benthic-plots------------------------------------------------------------
#| label: benthic-plots
#| fig-cap: "Benthic transect DHW vs relative change by era and region."

ggplot(dat.mod.bent.hc, aes(x = MaxDHW.mean, y = Rel.Change, fill = Era, colour = Era)) +
  geom_point() +
  scale_fill_viridis_d(end = 0.8) + 
  scale_color_viridis_d(end = 0.8) +
  coord_cartesian(xlim = c(0, 12), ylim = c(-1, 1)) +
  scale_y_continuous(labels = scales::percent) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) +
  theme_bw()


## ----repeat-exposure----------------------------------------------------------
#| label: repeat-exposure

dat.mod.repeat <- dat.mod.mant %>%
  filter(report_year > 2015) %>%
  group_by(Reef_Name) %>%
  arrange(report_year, .by_group = TRUE) %>%
  mutate(
    exposed4 = MaxDHW.mean > 4,
    ExposureNo = cumsum(exposed4),                 
    ExposureCat = case_when(
      exposed4 & ExposureNo == 1 ~ "First",
      exposed4 & ExposureNo  > 1 ~ "Repeat",
      TRUE ~ "NotExposed"
    ),
    StressCat = case_when(
      MaxDHW.mean < 4 ~ "0-4DHW",
      MaxDHW.mean < 8 ~ "4-8DHW",
      TRUE ~ "8-14DHW"
    )
  ) %>%
  ungroup()

dat_exp <- dat.mod.repeat %>%
  filter(ExposureCat %in% c("First", "Repeat")) %>%
  filter(is.finite(Abs.Change), is.finite(MaxDHW.mean))

# Mixed effects model to account for repeated measures on reefs
m1_re <- lmer(Abs.Change ~ ExposureCat + MaxDHW.mean + Region + (1|Reef_Name), data = dat_exp)
summary(m1_re)


## ----repeat-exposure-plots----------------------------------------------------
#| label: repeat-exposure-plots
#| fig-cap: "Relative coral cover change by first vs repeat DHW exposure >4."

Fig4Repeat <- ggplot(dat_exp |> filter(!is.na(ExposureCat)), 
       aes(x = MaxDHW.mean, y = Rel.Change, fill = ExposureCat, shape = ExposureCat, colour = ExposureCat)) +
  geom_point(alpha = 0.6) +
  scale_fill_manual(values = cols_exposure) +
  scale_color_manual(values = cols_exposure) +
  coord_cartesian(xlim = c(4, 12)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 3)) + 
  scale_y_continuous(labels = scales::percent) +
  labs(x = "DHW (°C-weeks)", y = "Relative coral cover change",
       fill = "Exposure", shape = "Exposure", colour = "Exposure") +
  theme_bw() +
  theme(legend.position = "bottom")

Fig4Repeat


## ----frequency-analysis-------------------------------------------------------
#| label: frequency-analysis

# Total cover change vs frequency of >8 DHW since 2015
reef_keys <- dat.mod.repeat %>%
  distinct(ReefName) %>% filter(!is.na(ReefName))

dhw_2015 <- dat.DHW %>%
  filter(ReefName %in% reef_keys$ReefName, Year >= 2015) %>%
  arrange(ReefName, Year) %>%
  group_by(ReefName) %>%
  mutate(ev8 = MaxDHW.mean > 8) %>%
  ungroup()

exposure_summary <- dhw_2015 %>%
  group_by(ReefName) %>%
  summarise(
    dhw_years = n_distinct(Year),
    n_ev8 = sum(ev8, na.rm = TRUE),
    freq_ev8 = n_ev8 / dhw_years,
    .groups = "drop"
  )

reef_summary <- dat.mod.repeat %>%
  filter(report_year > 2015) %>%
  filter(!is.na(ReefName), is.finite(mean)) %>%
  group_by(ReefName) %>%
  summarise(
    Region = first(Region),
    coral_years = n_distinct(report_year),
    year_start = min(report_year, na.rm = TRUE),
    year_end   = max(report_year, na.rm = TRUE),
    cover_start = mean[which.min(report_year)],
    cover_end   = mean[which.max(report_year)],
    rel_total_change   = (cover_end - cover_start) / cover_start,
    total_cover_change = cover_end - cover_start,
    .groups = "drop"
  ) %>%
  filter(coral_years >= 3, is.finite(rel_total_change), cover_start > 0, year_end > 2023)

reef_disturbance_df <- read_csv("Data/aims_ltmp/reef_disturbance.csv")
dist_summary <- reef_disturbance_df %>%
  filter(year > 2015, sample_type == "MANTA", disturbance == "c") %>%
  group_by(aims_reef_name, disturbance) %>%
  summarise(n_disturbances = n(), .groups = "drop") %>%
  left_join(dat.AIMSRef_key |> dplyr::select(aims_reef_name = AIMS_REEF_NAME_case, ReefName), 
            by = "aims_reef_name") 

reef_summary2 <- reef_summary %>%
  left_join(exposure_summary, by = "ReefName") %>%
  left_join(dist_summary, by = "ReefName") %>%
  mutate(
    COTSOutbreak = case_when(
      n_disturbances > 0 ~ "COTS Outbreak", 
      is.na(n_disturbances) ~ "No Outbreak"
    ),
    BleachingFreq = case_when(
      freq_ev8 == 0 ~ "0",
      freq_ev8 < 0.2 ~ "0-0.2/year",
      TRUE ~ "0.2/year +"
    )
  )

Fig4Freq <- ggplot(reef_summary2, aes(BleachingFreq, total_cover_change, fill = COTSOutbreak)) +
  geom_boxplot() +
  scale_fill_manual(values = cols_outbreak) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "Frequency of >8 DHW years since 2015",
       y = "Total coral cover change since 2015", fill = NULL) +
  theme(legend.position = "bottom")

Fig4Freq


## ----spatial-mapping----------------------------------------------------------
#| label: spatial-mapping
#| fig-cap: "Spatial distribution of thermal stress frequency ( >8 DHW) across the GBR."

shp <- list.files("data/SectorShapefile", pattern = "\\.shp$", full.names = TRUE)[1]
sectors <- st_read(shp, quiet = TRUE) |> st_make_valid()
sect_id <- intersect(names(sectors), c("SECT_NAME","SECTOR","Sector","NAME","Name"))[1]

target_crs <- st_crs(sectors)
coast_sp <- st_transform(coast, target_crs)
reefs_sp <- st_transform(reefs, target_crs)
reef_id <- intersect(names(reefs_sp), c("ReefName","REEF_NAME","LOC_NAME_S","domain_name","name"))[1]

reef_pts <- st_point_on_surface(reefs_sp) |>
  mutate(ReefName = .data[[reef_id]]) |>
  left_join(reef_summary2, by = "ReefName") |>
  filter(is.finite(freq_ev8), is.finite(rel_total_change))

reef_pts2 <- st_join(reef_pts, sectors |> dplyr::select(!!sect_id), left = TRUE)

sector_stats <- reef_pts2 |>
  st_drop_geometry() |>
  group_by(sector = .data[[sect_id]]) |>
  summarise(
    n_reefs = n(),
    freq_mean = mean(freq_ev8, na.rm = TRUE),
    freq_max  = max(freq_ev8, na.rm = TRUE),
    chg_mean  = mean(total_cover_change, na.rm = TRUE),
    chg_se    = sd(total_cover_change, na.rm = TRUE) / sqrt(sum(is.finite(total_cover_change))),
    .groups = "drop"
  )

sectors2 <- sectors |>
  mutate(sector = .data[[sect_id]]) |>
  left_join(sector_stats, by = "sector")

sector_labs <- sectors2 |>
  filter(is.finite(freq_mean)) |>
  mutate(
    lab_pt = st_point_on_surface(geometry),
    label = sprintf("freq>8: %.2f\nΔcover: %+.0f ± %.0f%%", 
                    freq_mean, 100*chg_mean, 100*chg_se)
  )

pal <- c("#2A9D8F", "#E9C46A", "#E76F51")

FigMapFreq <- ggplot() +
  geom_sf(data = coast_sp, fill = "grey92", colour = "grey50", linewidth = 0.2) +
  geom_sf(data = sectors2, aes(fill = freq_mean), colour = "grey25", linewidth = 0.25, alpha = 0.95) +
  geom_sf(data = reef_pts2, aes(colour = freq_ev8, size = freq_ev8*10), alpha = 0.9) +
  geom_sf_text(data = sector_labs, aes(geometry = lab_pt, label = label), size = 3.0, lineheight = 0.95) +
  scale_fill_gradientn(colours = pal, limits = c(0, 0.5), name = "Freq >8 DHW") +
  scale_colour_gradientn(colours = pal, limits = c(0, 0.5), name = "Freq >8 DHW") +
  scale_size(range = c(0.8, 4.2), limits = c(0, 1), guide = "none") +
  coord_sf(expand = FALSE) +
  theme_classic(base_size = 12) +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank(), legend.position = "right")

FigMapFreq


## ----ml-load-cheung-----------------------------------------------------------
#| label: ml-load-cheung
#| cache: true

load("data/Cheungetal2025/01_sstvar_blchrf.RData")

# Target predictors from Cheung et al.
predictor_vars <- c(
  "DHW",                                              # DHW
 "mcur_90",                                           # current speed (90-day)
 "cloudp_90",                                         # cloud fraction (90-day)
 "secc3m",                                            # Secchi depth (3-month)
 "cbclus2", "cbclus3", "cbclus4", "cbclus5",          # SST trajectory clusters
 "cbclus6", "cbclus7", "cbclus8", "cbclus9", "cbclus10",
 "histmDHW6", "yrsince6",                             # bleaching history (>6 DHW)
 "histmDHW4", "yrsince4",                             # bleaching history (>4 DHW)
 "ann_maxdhw_mx", "ann_maxdhw",                       # annual max DHW stats
 "winyear_sd", "winyear_mean"                         # winter SST variability
)

# Aggregate to reef-level means per LABEL_id × year
cheung_predictors <- sstvar_blch2 %>%
  select(LABEL_id, year, lat = Y, lon = X, Sector, all_of(predictor_vars)) %>%
  group_by(LABEL_id, year) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(year = as.numeric(year))

cat("Cheung predictors:", nrow(cheung_predictors), "reef-years,",
    length(unique(cheung_predictors$LABEL_id)), "unique reefs\n")
cat("Years:", paste(sort(unique(cheung_predictors$year)), collapse = ", "), "\n")
cat("Missing data summary:\n")
print(colSums(is.na(cheung_predictors)))


## ----ml-data-join-------------------------------------------------------------
#| label: ml-data-join

# 1. Manta tow data
dat.ml.mant <- dat.mod.mant %>%
  ungroup() %>%
  filter(!is.na(ReefID), !is.na(DHWYear)) %>%
  mutate(DHWYear = as.numeric(DHWYear)) %>%
  left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>%
  filter(!is.na(DHW)) %>%
  mutate(survey_depth = 7.5, survey_type = "manta") # Assume manta tow is 6-9m

# 2. Benthic + MMP data
dat.ml.bent <- dat.mod.bent.hc %>%
  ungroup() %>%
  filter(!is.na(ReefID), !is.na(DHWYear)) %>%
  mutate(DHWYear = as.numeric(DHWYear)) %>%
  left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>%
  filter(!is.na(DHW)) %>%
  mutate(survey_depth = as.numeric(depth), survey_type = "benthic")

# 3. Combined dataset
dat.ml.all <- bind_rows(dat.ml.mant, dat.ml.bent) %>%
  mutate(survey_type = as.factor(survey_type))

cat("Manta ML dataset:", nrow(dat.ml.mant), "observations\n")
cat("Benthic ML dataset:", nrow(dat.ml.bent), "observations\n")
cat("Combined ML dataset:", nrow(dat.ml.all), "observations\n")


## ----ml-join-summary----------------------------------------------------------
#| label: ml-join-summary

datatable(
  dat.ml.all %>%
    select(ReefName, DHWYear, survey_type, survey_depth, Abs.Change, Rel.Change,
           DHW, mcur_90, cloudp_90, histmDHW6, ann_maxdhw) %>%
    slice_head(n = 50) %>%
    mutate(across(where(is.numeric), ~ round(., 3))),
  caption = "Sample of Combined ML dataset with Cheung predictors"
)


## ----ml-correlation-----------------------------------------------------------
#| label: ml-correlation
#| fig-cap: "Correlation matrix of environmental predictors and mortality response (Combined dataset)."
#| fig-height: 8

library(corrplot)

cor_vars <- c("Abs.Change", "Rel.Change", "MaxDHW.mean", "survey_depth",
              predictor_vars[predictor_vars %in% names(dat.ml.all)])

cor_data <- dat.ml.all %>%
  ungroup() %>%
  select(all_of(cor_vars)) %>%
  select(where(~ sum(!is.na(.x)) > 10)) %>%
  drop_na()

if (nrow(cor_data) > 10) {
  cor_mat <- cor(cor_data, use = "pairwise.complete.obs")
  corrplot(cor_mat, method = "color", type = "lower", tl.cex = 0.7,
           order = "hclust", addrect = 4,
           title = "Predictor Correlations (Combined Data)", mar = c(0, 0, 2, 0))
} else {
  cat("Insufficient matched data for correlation analysis.\n")
}


## ----ml-modelling-------------------------------------------------------------
#| label: ml-modelling
#| cache: true

library(ranger)
library(gbm)
library(tidymodels)

# Features to include in the models
features <- c("MaxDHW.mean", "survey_depth", predictor_vars[predictor_vars %in% names(dat.ml.all)])

# Function to fit models and extract metrics
fit_eval_models <- function(dataset, response_var, dataset_name) {
  df <- dataset %>% 
    select(all_of(c(response_var, features, "Region"))) %>%
    drop_na()
  
  if(nrow(df) < 30) {
    return(data.frame(Dataset = dataset_name, Response = response_var, Model = "Insufficient Data", R2 = NA))
  }
  
  # 1. Random Forest (via ranger)
  rf_formula <- as.formula(paste(response_var, "~", paste(features, collapse = " + ")))
  rf_fit <- ranger(rf_formula, data = df, importance = "permutation", num.trees = 500, seed = 42)
  
  # 2. Baseline DHW-only RF
  rf_baseline_formula <- as.formula(paste(response_var, "~ MaxDHW.mean"))
  rf_base <- ranger(rf_baseline_formula, data = df, num.trees = 500, seed = 42)
  
  # 3. GBM
  gbm_fit <- gbm(rf_formula, data = df, distribution = "gaussian", 
                 n.trees = 500, interaction.depth = 3, shrinkage = 0.01, cv.folds = 5, n.cores = 1)
  best_iter <- gbm.perf(gbm_fit, method = "cv", plot.it = FALSE)
  gbm_pred <- predict(gbm_fit, df, n.trees = best_iter)
  gbm_r2 <- cor(df[[response_var]], gbm_pred)^2
  
  # Return metrics
  metrics <- data.frame(
    Dataset = dataset_name,
    Response = response_var,
    RF_Multi_R2 = rf_fit$r.squared,
    RF_Baseline_R2 = rf_base$r.squared,
    GBM_R2 = gbm_r2
  )
  
  # Store the best RF model in global env for later plotting
  assign(paste0("rf_", tolower(dataset_name), "_", tolower(sub("\\.", "", response_var))), rf_fit, envir = .GlobalEnv)
  
  return(metrics)
}

results <- bind_rows(
  fit_eval_models(dat.ml.mant, "Abs.Change", "Manta"),
  fit_eval_models(dat.ml.mant, "Rel.Change", "Manta"),
  fit_eval_models(dat.ml.bent, "Abs.Change", "Benthic"),
  fit_eval_models(dat.ml.bent, "Rel.Change", "Benthic"),
  fit_eval_models(dat.ml.all, "Abs.Change", "Combined"),
  fit_eval_models(dat.ml.all, "Rel.Change", "Combined")
)

datatable(results %>% mutate(across(where(is.numeric), ~ round(., 3))),
          caption = "Machine Learning Model Performance (OOB/CV R²)")


## ----ml-variable-importance---------------------------------------------------
#| label: ml-variable-importance
#| fig-cap: "Variable importance from Random Forest models (Combined Dataset) for absolute and relative change."

if (exists("rf_combined_abschange") && exists("rf_combined_relchange")) {
  vi_abs <- data.frame(
    Variable = names(rf_combined_abschange$variable.importance),
    Importance = rf_combined_abschange$variable.importance,
    Response = "Absolute Change"
  )
  vi_rel <- data.frame(
    Variable = names(rf_combined_relchange$variable.importance),
    Importance = rf_combined_relchange$variable.importance,
    Response = "Relative Change"
  )

  vi_all <- bind_rows(vi_abs, vi_rel)

  ggplot(vi_all, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    facet_wrap(~Response, scales = "free_x") +
    theme_bw() +
    labs(x = NULL, y = "Permutation Importance",
         title = "Variable Importance: Combined RF Models")
}


## ----ml-partial-dependence----------------------------------------------------
#| label: ml-partial-dependence
#| fig-cap: "Partial dependence of absolute coral cover change on key predictors."
#| fig-height: 8
#| fig-width: 10

library(pdp)

if (exists("rf_combined_abschange")) {
  # Get top 4 predictors for Absolute Change
  vi <- rf_combined_abschange$variable.importance
  top_vars <- names(sort(vi, decreasing = TRUE))[1:4]
  
  df_pd <- dat.ml.all %>% select(Abs.Change, all_of(features)) %>% drop_na()
  
  pd_plots <- list()
  for (var in top_vars) {
    pd <- partial(rf_combined_abschange, pred.var = var, train = df_pd)
    p <- ggplot(pd, aes_string(x = var, y = "yhat")) +
      geom_line(size = 1, color = "darkred") +
      theme_bw() +
      labs(y = "Partial Dependence (Abs. Change)")
    pd_plots[[var]] <- p
  }
  
  library(patchwork)
  wrap_plots(pd_plots, ncol = 2)
}

