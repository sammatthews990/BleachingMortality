# Diagnostic stocktake script
suppressPackageStartupMessages({
  library(dplyr)
  library(brms)
})

# Load datasets if pre-saved, or run quick inspection of QMD environment
qmd_path <- "DHW_Mortality_Attribution.qmd"

cat("=== Checking Dataset Sizes and Filtering ===\n")
# Read QMD lines to inspect dat.mod.mant and dat.ml.mant
# We can source or run R code directly from data files
load("Data/aims_ltmp/aims_ltmp.RData")

dat.DHW <- read.csv("Data/DHW_1985_2024_GBRReefs.csv") %>% rename("ReefName" = "LOC_NAME_S")
dat.AIMSRef <- read.csv("Data/AIMS-Reef_Reference.csv") %>% rename(Reef_Name = AIMS_REEF_NAME)
Disturbances_manta <- read.csv("Data/Disturbances_manta.csv") %>%
  left_join(dplyr::select(dat.AIMSRef, Reef_Name, ReefName), by = c("AIMS_REEF_NAME" = "Reef_Name")) %>%
  rename("report_year" = "REPORT_YEAR")

# Count rows in dat.mod.mant
df.AIMS.Mant <- reef_manta_df %>%
  dplyr::select(domain_name, project_code, report_year, date, depth, lower:median, mean, reef_zone, id, shelf, reefpage_category) %>%
  mutate(Reef_Name = toupper(domain_name))

dat.AIMSRef_key <- dat.AIMSRef %>%
  filter(!ReefName %in% c("Round-Russell Reef (17-013)", "Snake Reef (14-087)"))

df.MANT <- df.AIMS.Mant %>%
  left_join(dplyr::select(dat.AIMSRef_key, AIMS_REEF_NAME_cap, ReefID, ReefName, SECT_NAME, SECTOR),
            by = join_by("Reef_Name" == "AIMS_REEF_NAME_cap")) %>%
  mutate(
    year = floor(as.numeric(date)),
    frac = as.numeric(date) - year,
    date_approx = as.Date(paste0(year, "-01-01")) +
      round(frac * as.numeric(as.Date(paste0(year + 1, "-01-01")) - as.Date(paste0(year, "-01-01"))))
  )

dat.yrs.mant <- df.MANT %>%
  ungroup() %>%
  filter(!is.na(ReefName)) %>%
  mutate(
    Date = date_approx,
    Year = lubridate::year(Date), Month = lubridate::month(Date),
    DHWYear = ifelse(Month >= 7, Year, Year - 1)
  ) %>%
  left_join(dplyr::select(dat.DHW, ReefName, Year, MaxDHW.mean), by = c("ReefName", "DHWYear" = "Year"))

df.MANT <- df.MANT %>%
  arrange(ReefName, Reef_Name, reefpage_category, depth, report_year) %>%
  group_by(ReefName, Reef_Name, depth, reefpage_category) %>%
  mutate(
    Change = mean - lag(mean, n = 1),
    Rel.Change = Change / lag(mean, n = 1),
    Lag = report_year - lag(report_year, n = 1)
  ) %>%
  left_join(dat.yrs.mant)

dat.mod.mant <- df.MANT %>%
  filter(Lag < 3, project_code != "MMP", abs(Change) > 0.01) %>%
  left_join(dplyr::select(Disturbances_manta, ReefName, report_year, DISTURBANCE_TYPE),
            by = c("ReefName", "report_year")) %>%
  filter(
    !DISTURBANCE_TYPE %in% c("c", "d", "s", "f"),
    report_year %in% c(1999, 2003, 2017, 2018, 2021, 2023, 2024, 2025)
  )

cat("dat.mod.mant total rows:", nrow(dat.mod.mant), "\n")
cat("dat.mod.mant years present:", paste(sort(unique(dat.mod.mant$DHWYear)), collapse = ", "), "\n")

# Now check Cheung predictors join
cheung_files <- list.files("Data", pattern = "cheung|Cheung", full.names = TRUE)
cat("Cheung predictor files found:", paste(cheung_files, collapse = ", "), "\n")

# Load Cheung predictors if arrow/parquet available
if (file.exists("Data/cheung_predictors.csv") || file.exists("Data/cheung_predictors.parquet")) {
  if (file.exists("Data/cheung_predictors.parquet")) {
    cheung <- arrow::read_parquet("Data/cheung_predictors.parquet")
  } else {
    cheung <- read.csv("Data/cheung_predictors.csv")
  }
  cat("Cheung years:", paste(sort(unique(cheung$year)), collapse = ", "), "\n")
} else {
  # Search inside Data for parquet/csv files
  data_files <- list.files("Data", recursive = TRUE, full.names = TRUE)
  cheung_match <- grep("cheung", data_files, ignore.case = TRUE, value = TRUE)
  cat("Matched Cheung files:", paste(cheung_match, collapse = ", "), "\n")
}
