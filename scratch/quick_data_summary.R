# Mock reefs_join and run data prep using dat.AIMSRef and COTSModParams
library(tidyr)
library(dplyr)
library(janitor)
library(lubridate)
library(readr)
library(stringr)

load("Data/aims_ltmp/aims_ltmp.RData")
load("data/Cheungetal2025/01_sstvar_blchrf.RData")

dat.AIMSRef <- read_csv("Data/AIMS-Reef_Reference.csv") %>%
  rename(Reef_Name = AIMS_REEF_NAME)

Disturbances_manta <- read_csv("Data/Disturbances_manta.csv") %>%
  left_join(select(dat.AIMSRef, Reef_Name, ReefName), by = c("AIMS_REEF_NAME" = "Reef_Name")) %>%
  rename("report_year" = "REPORT_YEAR")

# DHW data
dat.DHW <- read_csv("Data/DHW_1985_2024_GBRReefs.csv") %>%
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

df.AIMS.Mant <- reef_manta_df |>
  select(domain_name, project_code, report_year, date, depth,
         lower:median, mean, reef_zone, id, shelf, reefpage_category) |>
  mutate(Reef_Name = toupper(domain_name))

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

# Growth params
dat.grid <- read_csv("Data/COTSModParams.csv")
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

doCoralGrowth.Simp <- function(CoralCover, B0, HC.asym, WQ) {
  WQ.mn.sd <- c(-0.68, 0.03)
  b0.wq <- B0 + WQ * WQ.mn.sd[1]
  b1.wq <- b0.wq / log(HC.asym)
  CoralCover <- log(CoralCover)
  CoralCover[which(!is.finite(CoralCover))] <- NA
  cc_next <- exp(log(HC.asym) * (1 - exp(b1.wq)) + CoralCover * exp(b1.wq))
  cc_next[which(is.na(cc_next))] <- 0
  return(c(exp(CoralCover), cc_next, cc_next - exp(CoralCover)))
}

print(paste("Rows in df.MANT before filter:", nrow(df.MANT)))
print(paste("Unique reefs in df.MANT:", length(unique(df.MANT$ReefName))))
print(paste("Unique reefs in dat.gridSum:", length(unique(dat.gridSum$REEF_NAME))))

df.MANT <- df.MANT |>
  mutate(Change.Cor = Change, Prev.HC = mean - Change.Cor) |>
  filter(ReefName %in% dat.gridSum$REEF_NAME, !is.na(Lag),
         Reef_Name != "MILLN REEF")

print(paste("Rows in df.MANT after filter:", nrow(df.MANT)))

df.MANT <- df.MANT |>
  rowwise() |>
  mutate(
    idx = which(dat.gridSum$REEF_NAME == ReefName)[1],
    B0_val = dat.gridSum$B0[idx],
    WQ_val = dat.gridSum$WQ[idx],
    HC_val = dat.gridSum$HC.asym[idx],
    Pred.Grth = doCoralGrowth.Simp(
      CoralCover = Prev.HC * 100,
      B0 = B0_val,
      WQ = WQ_val,
      HC.asym = HC_val
    )[3] / 100
  ) |>
  mutate(
    Pred.Grth2 = ifelse(Lag == 2,
      Pred.Grth + doCoralGrowth.Simp(
        CoralCover = (Prev.HC + Pred.Grth) * 100,
        B0 = B0_val,
        WQ = WQ_val,
        HC.asym = HC_val
      )[3] / 100,
      Pred.Grth),
    Abs.Change = Change.Cor - Pred.Grth2,
    Rel.Change = Abs.Change / (Prev.HC + Pred.Grth2)
  )

# MOCK reefs_join!
reefs_join <- data.frame(
  LOC_NAME_S = unique(df.MANT$ReefName),
  Region = "Central"
)

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

# Prepare cheung predictors
predictor_vars <- c(
  "mcur_90", "cloudp_90", "secc3m", "cbclus2", "histmDHW6", "yrsince6", "histmDHW4", "yrsince4", "winyear_sd", "winyear_mean"
)
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
  filter(!is.na(mcur_90)) %>%
  mutate(survey_depth = 7.5, survey_type = "manta")

dat.ml.mant.clean <- dat.ml.mant %>%
  filter(!is.na(mcur_90))

print("=== Manta Tow Data Cleaned ===")
print(nrow(dat.ml.mant.clean))
print(summary(dat.ml.mant.clean$MaxDHW.mean))

# Scale params
scale_params <- list(
  MaxDHW.mean = c(mean = mean(dat.ml.mant.clean$MaxDHW.mean, na.rm = TRUE),
                  sd2  = 2 * sd(dat.ml.mant.clean$MaxDHW.mean, na.rm = TRUE)),
  histmDHW6   = c(mean = mean(dat.ml.mant.clean$histmDHW6, na.rm = TRUE),
                  sd2  = 2 * sd(dat.ml.mant.clean$histmDHW6, na.rm = TRUE))
)
print("=== scale_params ===")
print(scale_params)
