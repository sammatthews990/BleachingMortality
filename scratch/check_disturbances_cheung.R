library(dplyr)

# Load data
load("Data/aims_ltmp/aims_ltmp.RData")
dat.AIMSRef <- read.csv("Data/AIMS-Reef_Reference.csv") %>% rename(Reef_Name = AIMS_REEF_NAME)
Disturbances_manta <- read.csv("Data/Disturbances_manta.csv") %>%
  left_join(select(dat.AIMSRef, Reef_Name, ReefName), by = c("AIMS_REEF_NAME" = "Reef_Name")) %>%
  rename("report_year" = "REPORT_YEAR")
dat.DHW <- read.csv("Data/DHW_1985_2024_GBRReefs.csv") %>% rename("ReefName" = "LOC_NAME_S")

# Growth params
dat.grid <- read.csv("Data/COTSModParams.csv")
dat.gridSum <- dat.grid %>%
  mutate(REEF_ID = ifelse(grepl("14-116", REEF_ID, fixed = TRUE), "14-116", REEF_ID),
         REEF_NAME = ifelse(grepl("14-116", REEF_NAME, fixed = TRUE), "Lizard Island Reef (14-116)", REEF_NAME)) %>%
  group_by(REEF_NAME, REEF_ID, SECTOR) %>%
  summarise(B0 = mean(pred.b0.mean), HC.asym = mean(pred.HCmax.mean),
            WQ1 = mean(Primary), WQ2 = mean(Secondary), WQ3 = mean(Tertiary), .groups = "drop") %>%
  mutate(WQ = WQ1 + WQ2 + WQ3)

doCoralGrowth.Simp <- function(CoralCover, B0, HC.asym, WQ) {
  WQ.mn.sd <- c(-0.68, 0.03)
  b0.wq <- B0 + WQ * WQ.mn.sd[1]
  b1.wq <- b0.wq / log(HC.asym)
  CoralCover <- log(CoralCover)
  b0.wq <- as.numeric(b0.wq); b1.wq <- as.numeric(b1.wq); CoralCover <- as.numeric(CoralCover)
  CoralCover.t1 <- (b0.wq + (1 - b1.wq) * CoralCover)
  CoralCover.t1 <- exp(CoralCover.t1)
  Change.t1 <- CoralCover.t1 - exp(CoralCover)
  return(c(exp(CoralCover), CoralCover.t1, Change.t1))
}

# Agincourt DHW fill
dat.agin <- dat.DHW %>% filter(grepl("Agincourt", ReefName)) %>% group_by(Year) %>% summarise(MaxDHW.mean = mean(MaxDHW.mean, na.rm = TRUE), .groups = "drop")
dat.DHW <- dat.DHW %>% left_join(dat.agin, by = "Year") %>% mutate(MaxDHW.mean = ifelse(grepl("Agincourt", ReefName) & is.na(MaxDHW.mean.x), MaxDHW.mean.y, MaxDHW.mean.x)) %>% select(-MaxDHW.mean.x, -MaxDHW.mean.y)

df.AIMS.Mant <- reef_manta_df %>% select(domain_name, project_code, report_year, date, depth, lower:median, mean, reef_zone, id, shelf, reefpage_category) %>% mutate(Reef_Name = toupper(domain_name))
dat.AIMSRef_key <- dat.AIMSRef %>% filter(!ReefName %in% c("Round-Russell Reef (17-013)", "Snake Reef (14-087)"))
df.MANT <- df.AIMS.Mant %>% left_join(select(dat.AIMSRef_key, AIMS_REEF_NAME_cap, ReefID, ReefName, SECT_NAME, SECTOR), by = join_by("Reef_Name" == "AIMS_REEF_NAME_cap")) %>% mutate(year = floor(as.numeric(date)), frac = as.numeric(date) - year, date_approx = lubridate::ymd(paste0(year, "-01-01")) + round(frac * as.numeric(lubridate::ymd(paste0(year + 1, "-01-01")) - lubridate::ymd(paste0(year, "-01-01")))))
dat.yrs.mant <- df.MANT %>% ungroup() %>% filter(!is.na(ReefName)) %>% mutate(Date = date_approx, Year = lubridate::year(Date), Month = lubridate::month(Date), DHWYear = ifelse(Month >= 7, Year, Year - 1)) %>% left_join(select(dat.DHW, ReefName, Year, MaxDHW.mean), by = c("ReefName", "Year" = "Year"))
df.MANT <- df.MANT %>% arrange(ReefName, Reef_Name, reefpage_category, depth, report_year) %>% group_by(ReefName, Reef_Name, depth, reefpage_category) %>% mutate(Change = mean - lag(mean, n = 1), Rel.Change = Change / lag(mean, n = 1), Lag = report_year - lag(report_year, n = 1)) %>% left_join(dat.yrs.mant)
df.MANT <- df.MANT %>% mutate(Change.Cor = Change, Prev.HC = mean - Change.Cor) %>% filter(ReefName %in% dat.gridSum$REEF_NAME, !is.na(Lag), Reef_Name != "MILLN REEF") %>% rowwise() %>% mutate(Pred.Grth = doCoralGrowth.Simp(CoralCover = Prev.HC * 100, B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"], WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"], HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"])[3] / 100) %>% mutate(Pred.Grth2 = ifelse(Lag == 2, Pred.Grth + doCoralGrowth.Simp(CoralCover = (Prev.HC + Pred.Grth) * 100, B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"], WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"], HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"])[3] / 100, Pred.Grth), Abs.Change = Change.Cor - Pred.Grth2, Rel.Change = Abs.Change / (Prev.HC + Pred.Grth2))

dat.mod.mant.raw <- df.MANT %>% filter(Lag < 3, project_code != "MMP", abs(Change) > 0.01) %>% left_join(select(Disturbances_manta, ReefName, report_year, DISTURBANCE_TYPE), by = c("ReefName", "report_year")) %>% filter(report_year %in% c(1999, 2003, 2017, 2018, 2021, 2023, 2024, 2025))

load("data/Cheungetal2025/01_sstvar_blchrf.RData")
predictor_vars <- c("mcur_90", "cloudp_90", "secc3m", "cbclus2", "histmDHW6", "yrsince6", "histmDHW4", "yrsince4", "winyear_sd", "winyear_mean")
cheung_predictors <- sstvar_blch2 %>% select(LABEL_id, year, lat = Y, lon = X, Sector, all_of(predictor_vars)) %>% group_by(LABEL_id, year) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>% mutate(year = as.numeric(year))

dat.ml.mant.raw <- dat.mod.mant.raw %>% ungroup() %>% filter(!is.na(ReefID), !is.na(DHWYear)) %>% mutate(DHWYear = as.numeric(DHWYear)) %>% left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>% filter(!is.na(mcur_90))

cat("=== Disturbance Types in Full Dataset (N=516 raw) ===\n")
print(table(dat.mod.mant.raw$DISTURBANCE_TYPE, useNA = "always"))

cat("\n=== Disturbance Types in Cheung Predictor Subset (N=128 raw) ===\n")
print(table(dat.ml.mant.raw$DISTURBANCE_TYPE, useNA = "always"))

cat("\n=== Comparing Filtering Strategies on Cheung Subset ===\n")
cat("1. Old Filter (!c, d, s, f): N =", nrow(dat.ml.mant.raw %>% filter(!DISTURBANCE_TYPE %in% c("c", "d", "s", "f"))), "\n")
cat("2. New Explicit Filter (is.na | b, m, n): N =", nrow(dat.ml.mant.raw %>% filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n"))), "\n")
cat("3. Bleaching Only Filter (b, m): N =", nrow(dat.ml.mant.raw %>% filter(DISTURBANCE_TYPE %in% c("b", "m"))), "\n")
cat("4. Bleaching + No Disturbance (is.na | b, m, n): N =", nrow(dat.ml.mant.raw %>% filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n"))), "\n")
