library(dplyr)
library(gbm)

source("scratch/check_disturbances_cheung.R")

# Clean Manta Tow Training Dataset (N=110)
dat.train <- dat.ml.mant.raw %>%
  filter(is.na(DISTURBANCE_TYPE) | DISTURBANCE_TYPE %in% c("b", "m", "n")) %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1)
  ) %>%
  filter(is.finite(Mort.prop), is.finite(MaxDHW.mean))

predictor_vars <- c("MaxDHW.mean", "secc3m", "winyear_mean", "histmDHW6", "mcur_90", "winyear_sd", "cloudp_90", "yrsince6")
f_str <- as.formula(paste("Mort.prop ~", paste(predictor_vars, collapse = " + ")))

# 1. Unregularized BRT (Depth = 3)
set.seed(42)
brt_orig <- gbm(f_str, data = dat.train, distribution = "gaussian", n.trees = 1000, interaction.depth = 3, shrinkage = 0.01, n.minobsinnode = 5, bag.fraction = 0.5, cv.folds = 5)
best_orig <- gbm.perf(brt_orig, method = "cv", plot.it = FALSE)

# 2. Regularized Pairwise BRT (Depth = 2)
set.seed(42)
brt_reg2 <- gbm(f_str, data = dat.train, distribution = "gaussian", n.trees = 1500, interaction.depth = 2, shrinkage = 0.005, n.minobsinnode = 10, bag.fraction = 0.70, cv.folds = 5)
best_reg2 <- gbm.perf(brt_reg2, method = "cv", plot.it = FALSE)

# 3. Regularized Additive BRT (Depth = 1)
set.seed(42)
brt_reg1 <- gbm(f_str, data = dat.train, distribution = "gaussian", n.trees = 1500, interaction.depth = 1, shrinkage = 0.005, n.minobsinnode = 10, bag.fraction = 0.70, cv.folds = 5)
best_reg1 <- gbm.perf(brt_reg1, method = "cv", plot.it = FALSE)

load("Data/aims_ltmp/aims_ltmp.RData")
df.AIMS.full <- reef_photo_df

df.AIMS.Bent <- df.AIMS.full %>%
  filter(data_type == "photo-transect", domain_category == "reef", project_code %in% c("LTMP", "MMP"), purpose == "GROUP_LEVEL", variable == "HARD CORAL") %>%
  mutate(Reef_Name = toupper(domain_name), Reef_Name_Clean = gsub(" REEF(S)?| ISLAND| IS", "", Reef_Name)) %>%
  select(Reef_Name, Reef_Name_Clean, project_code, report_year, date, depth, lower:median, mean, reef_zone, id, shelf, reefpage_category)

dat.AIMSRef_key <- dat.AIMSRef %>% filter(!ReefName %in% c("Round-Russell Reef (17-013)", "Snake Reef (14-087)"))

df.BENT <- df.AIMS.Bent %>%
  left_join(select(dat.AIMSRef_key, AIMS_REEF_NAME_cap, ReefID, ReefName, SECT_NAME, SECTOR), by = join_by("Reef_Name" == "AIMS_REEF_NAME_cap")) %>%
  mutate(year = floor(as.numeric(date)), frac = as.numeric(date) - year, date_approx = lubridate::ymd(paste0(year, "-01-01")) + round(frac * as.numeric(lubridate::ymd(paste0(year + 1, "-01-01")) - lubridate::ymd(paste0(year, "-01-01")))))

dat.yrs.bent <- df.BENT %>% ungroup() %>% filter(!is.na(ReefName)) %>% mutate(Date = date_approx, Year = lubridate::year(Date), Month = lubridate::month(Date), DHWYear = ifelse(Month >= 7, Year, Year - 1)) %>% left_join(select(dat.DHW, ReefName, Year, MaxDHW.mean), by = c("ReefName", "Year" = "Year"))

df.BENT <- df.BENT %>% arrange(ReefName, Reef_Name, depth, reefpage_category, report_year) %>% group_by(ReefName, Reef_Name, depth, reefpage_category) %>% mutate(Change = mean - lag(mean, n = 1), Rel.Change = Change / lag(mean, n = 1), Lag = report_year - lag(report_year, n = 1)) %>% left_join(dat.yrs.bent)

df.BENT <- df.BENT %>% mutate(Change.Cor = Change, Prev.HC = mean - Change.Cor) %>% filter(ReefName %in% dat.gridSum$REEF_NAME, !is.na(Lag), Reef_Name != "MILLN REEF") %>% rowwise() %>% mutate(Pred.Grth = doCoralGrowth.Simp(CoralCover = Prev.HC * 100, B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"], WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"], HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"])[3] / 100) %>% mutate(Pred.Grth2 = ifelse(Lag == 2, Pred.Grth + doCoralGrowth.Simp(CoralCover = (Prev.HC + Pred.Grth) * 100, B0 = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "B0"], WQ = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "WQ"], HC.asym = dat.gridSum[which(dat.gridSum$REEF_NAME == ReefName), "HC.asym"])[3] / 100, Pred.Grth), Abs.Change = Change.Cor - Pred.Grth2, Rel.Change = Abs.Change / (Prev.HC + Pred.Grth2))

dat.ml.ltmp <- df.BENT %>% ungroup() %>% filter(project_code == "LTMP", !is.na(ReefID), !is.na(DHWYear)) %>% mutate(DHWYear = as.numeric(DHWYear)) %>% left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>% filter(!is.na(mcur_90))

dat.ml.ltmp.clean <- dat.ml.ltmp %>% filter(!is.na(mcur_90), !is.na(Rel.Change), !is.na(MaxDHW.mean)) %>% mutate(Mort.prop = pmin(pmax(-Rel.Change, 0), 1)) %>% filter(is.finite(Mort.prop))

Y_test <- dat.ml.ltmp.clean$Mort.prop
SS_tot_test <- sum((Y_test - mean(Y_test))^2)
test_df <- dat.ml.ltmp.clean %>% select(all_of(predictor_vars))

eval_test <- function(fit, b_iter, name) {
  p <- predict(fit, test_df, n.trees = b_iter, type = "response")
  p_b <- pmin(pmax(p, 0), 1)
  dev <- (1 - sum((Y_test - p_b)^2) / SS_tot_test) * 100
  r2 <- max(0, cor(Y_test, p_b)^2)
  rmse <- sqrt(mean((Y_test - p_b)^2))
  data.frame(Model = name, Test_Deviance = sprintf("%.1f%%", dev), Test_R2 = sprintf("%.3f", r2), Test_RMSE = sprintf("%.4f", rmse))
}

res_test <- bind_rows(
  eval_test(brt_orig, best_orig, "Original Unregularized BRT (Depth=3)"),
  eval_test(brt_reg2, best_reg2, "Regularized Pairwise BRT (Depth=2)"),
  eval_test(brt_reg1, best_reg1, "Regularized Additive BRT (Depth=1)")
)

cat("\n=======================================================================\n")
cat("      INDEPENDENT LTMP TEST EVALUATION OF REGULARIZED BRTs (N=107)\n")
cat("=======================================================================\n")
print(res_test)

saveRDS(res_test, "output/models/brt_regularized_ltmp_test.rds")
