library(dplyr)

Disturbances_manta <- read.csv("Data/Disturbances_manta.csv")
cat("=== Disturbances_manta Summary ===\n")
print(table(Disturbances_manta$DISTURBANCE_TYPE, useNA = "always"))

load("Data/aims_ltmp/aims_ltmp.RData")
dat.AIMSRef <- read.csv("Data/AIMS-Reef_Reference.csv") %>% rename(Reef_Name = AIMS_REEF_NAME)
Disturbances_manta <- Disturbances_manta %>%
  left_join(select(dat.AIMSRef, Reef_Name, ReefName), by = c("AIMS_REEF_NAME" = "Reef_Name")) %>%
  rename("report_year" = "REPORT_YEAR")

# Load full Manta Tow dataset prepared in QMD
dat.DHW <- read.csv("Data/DHW_1985_2024_GBRReefs.csv") %>% rename("ReefName" = "LOC_NAME_S")

# Check what disturbance codes exist and how current filter affects observations
cat("\n=== Checking Filter: !DISTURBANCE_TYPE %in% c('c','d','s','f') ===\n")
# vs Strict filter: is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m")

df.AIMS.Mant <- reef_manta_df %>%
  select(domain_name, project_code, report_year, date, depth, lower:median, mean, reef_zone, id, shelf, reefpage_category) %>%
  mutate(Reef_Name = toupper(domain_name))

dat.AIMSRef_key <- dat.AIMSRef %>% filter(!ReefName %in% c("Round-Russell Reef (17-013)", "Snake Reef (14-087)"))

df.MANT <- df.AIMS.Mant %>%
  left_join(select(dat.AIMSRef_key, AIMS_REEF_NAME_cap, ReefID, ReefName, SECT_NAME, SECTOR),
            by = join_by("Reef_Name" == "AIMS_REEF_NAME_cap")) %>%
  mutate(year = floor(as.numeric(date))) %>%
  arrange(ReefName, Reef_Name, reefpage_category, depth, report_year) %>%
  group_by(ReefName, Reef_Name, depth, reefpage_category) %>%
  mutate(
    Change = mean - lag(mean, n = 1),
    Rel.Change = Change / lag(mean, n = 1),
    Lag = report_year - lag(report_year, n = 1)
  )

df.MANT_dist <- df.MANT %>%
  filter(Lag < 3, project_code != "MMP", abs(Change) > 0.01) %>%
  left_join(select(Disturbances_manta, ReefName, report_year, DISTURBANCE_TYPE),
            by = c("ReefName", "report_year")) %>%
  filter(report_year %in% c(1999, 2003, 2017, 2018, 2021, 2023, 2024, 2025))

cat("Total observations in target report years:", nrow(df.MANT_dist), "\n")
cat("\nBreakdown of DISTURBANCE_TYPE in df.MANT_dist:\n")
print(table(df.MANT_dist$DISTURBANCE_TYPE, useNA = "always"))

df_curr <- df.MANT_dist %>% filter(!DISTURBANCE_TYPE %in% c("c", "d", "s", "f"))
cat("\nCurrent filter (!c, d, s, f) N =", nrow(df_curr), "\n")

df_strict <- df.MANT_dist %>% filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "B", "M"))
cat("Strict bleaching filter (NA or 'b' or 'm') N =", nrow(df_strict), "\n")

# Check if there are other codes like "u", "o", etc. that were slipping through
dropped_by_strict <- setdiff(rownames(df_curr), rownames(df_strict))
cat("\nObservations in current filter but NOT in strict bleaching filter:\n")
cat("Count:", nrow(df_curr) - nrow(df_strict), "\n")
print(table(df_curr$DISTURBANCE_TYPE, useNA = "always"))
