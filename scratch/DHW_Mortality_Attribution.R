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

# Force dplyr functions over MASS/stats
select <- dplyr::select
filter <- dplyr::filter


## ----helper-functions---------------------------------------------------------
#| label: helper-functions

# ── Gompertz coral growth function (simplified, single time-step) ──
# Returns c(current_cover, next_cover, change)
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

# ── Multi-step growth function (n years from starting cover) ──
doCoralGrowth <- function(n, CoralCover, B0, HC.asym, WQ) {
  WQ.mn.sd <- c(-0.68, 0.03)
  b0.wq <- B0 + WQ * WQ.mn.sd[1]
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


## ----manta-tow-data-----------------------------------------------------------
#| label: manta-tow-data
#| cache: true

load("Data/aims_ltmp/aims_ltmp.RData")

dat.AIMSRef <- read.csv("Data/AIMS-Reef_Reference.csv") %>%
  rename(Reef_Name = AIMS_REEF_NAME)

Disturbances_manta <- read.csv("Data/Disturbances_manta.csv") %>%
  left_join(select(dat.AIMSRef, Reef_Name, ReefName), by = c("AIMS_REEF_NAME" = "Reef_Name")) %>%
  rename("report_year" = "REPORT_YEAR")

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
#| cache: true

# Load spatial regions
gbr <- st_read("Data/GBR_AIMS") |> st_make_valid()
data_bucket <- s3_bucket("s3://gbr-dms-data-public/gbrmpa-complete-gbr-features/data.parquet")
reefs <- open_dataset(data_bucket) %>%
  filter(FEAT_NAME == "Reef", DATASET == "GBR Features") |>
  distinct(UNIQUE_ID, GBR_NAME, LOC_NAME_S, geometry, FEAT_NAME) %>%
  collect() %>% st_as_sf(crs = 4283) |> st_make_valid()

coast <- open_dataset(data_bucket) %>%
  filter(FEAT_NAME %in% c("Mainland", "Island"), DATASET == "GBR Features") |>
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


## ----model-extrapolation-fits-------------------------------------------------
#| label: model-extrapolation-fits

# ── ReefMod mortality function (Bozec et al.) ──
reefmod_M <- function(DHW, d = 7, s = 1, threshold = 3) {
  w <- 1 / (0.420 + 0.272 * d)
  m <- w * s * (exp(0.168 + 0.347 * DHW) - 1) / 100
  m <- pmin(1, pmax(0, m))
  m[DHW < threshold] <- 0
  M <- 1 - (1 - m)^6
  pmin(1, pmax(0, M))
}

# ── Helper: fit all four models and return prediction grid ──
fit_four_models <- function(dat, dhw_col, mort_col, max_dhw = 20) {
  d <- dat %>%
    transmute(DHW = .data[[dhw_col]], Mort = .data[[mort_col]]) %>%
    filter(is.finite(DHW), is.finite(Mort), Mort >= 0, Mort <= 1)
  
  # Nudge exact 0/1 for beta regression
  n <- nrow(d)
  d <- d %>% mutate(
    Mort_beta = pmin(pmax(Mort, 0), 1),
    Mort_beta = (Mort_beta * (n - 1) + 0.5) / n
  )
  
  grid <- data.frame(DHW = seq(0, max_dhw, by = 0.1))
  obs_max <- max(d$DHW, na.rm = TRUE)
  
  # 1. Unbounded GAM
  m1 <- gam(Mort ~ s(DHW, k = 4), data = d, method = "REML")
  pr1 <- predict(m1, newdata = grid, se.fit = TRUE)
  grid$unb_fit <- pr1$fit
  grid$unb_lwr <- pr1$fit - 1.96 * pr1$se.fit
  grid$unb_upr <- pr1$fit + 1.96 * pr1$se.fit
  
  # 2. Bounded GAM (quasibinomial logit)
  m2 <- gam(Mort ~ s(DHW, k = 4), data = d, family = quasibinomial(link = "logit"), method = "REML")
  pr2 <- predict(m2, newdata = grid, type = "link", se.fit = TRUE)
  grid$bnd_fit <- inv_logit(pr2$fit)
  grid$bnd_lwr <- inv_logit(pr2$fit - 1.96 * pr2$se.fit)
  grid$bnd_upr <- inv_logit(pr2$fit + 1.96 * pr2$se.fit)
  
  # 3. Beta regression
  m3 <- tryCatch(betareg(Mort_beta ~ DHW, data = d, link = "logit"), error = function(e) NULL)
  if (!is.null(m3)) {
    grid$beta_fit <- predict(m3, newdata = grid, type = "response")
    # Bootstrap CIs via coefficient simulation
    Xg <- model.matrix(~ DHW, data = grid)
    b_mean <- coef(m3)[c("(Intercept)", "DHW")]
    V_mean <- vcov(m3)[c("(Intercept)", "DHW"), c("(Intercept)", "DHW")]
    set.seed(42)
    b_sim <- MASS::mvrnorm(2000, mu = b_mean, Sigma = V_mean)
    mu_sim <- inv_logit(b_sim %*% t(Xg))
    grid$beta_lwr <- apply(mu_sim, 2, quantile, 0.025)
    grid$beta_upr <- apply(mu_sim, 2, quantile, 0.975)
  } else {
    grid$beta_fit <- grid$beta_lwr <- grid$beta_upr <- NA
  }
  
  # 4. Quasi-binomial GLM
  m4 <- glm(Mort ~ DHW, data = d, family = quasibinomial(link = "logit"))
  pr4 <- predict(m4, newdata = grid, type = "link", se.fit = TRUE)
  grid$glm_fit <- inv_logit(pr4$fit)
  grid$glm_lwr <- inv_logit(pr4$fit - 1.96 * pr4$se.fit)
  grid$glm_upr <- inv_logit(pr4$fit + 1.96 * pr4$se.fit)
  
  # ReefMod reference curves
  sens_s <- c(1.4, 1.5, 1.6, 1.7)
  grid$reefmod_sens     <- reefmod_M(grid$DHW, d = 7, s = mean(sens_s))
  grid$reefmod_sens_lwr <- reefmod_M(grid$DHW, d = 7, s = min(sens_s))
  grid$reefmod_sens_upr <- reefmod_M(grid$DHW, d = 7, s = max(sens_s))
  grid$reefmod_tol      <- reefmod_M(grid$DHW, d = 7, s = 0.25)
  
  # Calculate R2 values
  r2_vals <- list(
    unb_gam = summary(m1)$r.sq,
    bnd_gam = summary(m2)$r.sq,
    beta    = if (!is.null(m3)) m3$pseudo.r.squared else NA,
    glm     = 1 - (m4$deviance / m4$null.deviance)
  )
  
  list(data = d, grid = grid, obs_max = obs_max, r2 = r2_vals,
       models = list(unbounded_gam = m1, bounded_gam = m2, beta = m3, glm = m4))
}

# ── Fit to Hughes et al. (2018) data ──
hughes_fits <- fit_four_models(dat.bleach, "DHW", "Mort.bin", max_dhw = 20)

# ── Fit to AIMS LTMP manta tow data ──
# Transform growth-adjusted Rel.Change to mortality proportion
dat_aims_mort <- dat.mod.mant %>%
  ungroup() %>%
  filter(is.finite(MaxDHW.mean), is.finite(Rel.Change)) %>%
  mutate(Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
         Mort.prop = ifelse(Mort.prop == 0, 0.001, Mort.prop))

aims_fits <- fit_four_models(dat_aims_mort, "MaxDHW.mean", "Mort.prop", max_dhw = 25)


## ----model-extrapolation-plot-------------------------------------------------
#| label: model-extrapolation-plot
#| fig-cap: "Model choice & extrapolation sensitivity. Left: Hughes et al. (2018) single-event data. Right: AIMS LTMP manta tow data (growth-corrected, multi-event). Grey shading indicates extrapolation beyond the observed DHW range. Black lines show the ReefMod process-based mortality curve (solid = heat-sensitive, dashed = heat-tolerant)."
#| fig-height: 7
#| fig-width: 12

# ── Helper: build a plot from fits ──
build_extrap_plot <- function(fits, title_label, caption_text = NULL, x_limit = 20) {
  g <- fits$grid %>% filter(DHW <= x_limit)
  d <- fits$data
  x_max <- fits$obs_max
  r2 <- fits$r2
  
  # Reshape to long format for plotting
  pred_long <- bind_rows(
    g %>% transmute(DHW, model = "Bounded GAM (logit)", fit = bnd_fit, lwr = bnd_lwr, upr = bnd_upr),
    g %>% transmute(DHW, model = "Beta regression", fit = beta_fit, lwr = beta_lwr, upr = beta_upr),
    g %>% transmute(DHW, model = "Quasi-binomial GLM", fit = glm_fit, lwr = glm_lwr, upr = glm_upr)
  ) %>%
    mutate(model = factor(model, levels = c("Bounded GAM (logit)", "Beta regression", "Quasi-binomial GLM")))
  
  # Create R2 annotation text
  r2_text <- paste0(
    "R\u00B2 values:\n",
    "Bounded GAM: ", round(r2$bnd_gam, 3), "\n",
    "Beta Reg: ", round(r2$beta, 3), "\n",
    "GLM: ", round(r2$glm, 3)
  )
  
  p <- ggplot(d, aes(DHW, Mort)) +
    # Extrapolation shading
    annotate("rect", xmin = x_max, xmax = x_limit, ymin = -Inf, ymax = Inf,
             fill = "grey80", alpha = 0.5) +
    geom_vline(xintercept = x_max, linetype = "dotted", linewidth = 0.8) +
    annotate("text", x = x_max + (x_limit - x_max) * 0.5, y = 0.05,
             label = "Extrapolation", colour = "grey30", size = 3) +
    
    # Data points
    geom_point(shape = 17, colour = "grey30", size = 1.2, alpha = 0.3) +
    
    # Model ribbons and lines
    geom_ribbon(data = pred_long,
                aes(x = DHW, ymin = lwr, ymax = upr, fill = model),
                alpha = 0.12, inherit.aes = FALSE) +
    geom_line(data = pred_long,
              aes(x = DHW, y = fit, colour = model),
              linewidth = 0.9, inherit.aes = FALSE) +
    
    # ReefMod reference (heat-sensitive)
    geom_ribbon(data = g, aes(x = DHW, ymin = reefmod_sens_lwr, ymax = reefmod_sens_upr),
                fill = "black", alpha = 0.06, inherit.aes = FALSE) +
    geom_line(data = g, aes(x = DHW, y = reefmod_sens),
              colour = "black", linewidth = 1, inherit.aes = FALSE) +
    # ReefMod (heat-tolerant)
    geom_line(data = g, aes(x = DHW, y = reefmod_tol),
              colour = "black", linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE) +
    
    # R2 Annotation
    annotate("label", x = x_limit, y = 0, label = r2_text, 
             hjust = 1.05, vjust = -0.1, size = 3, label.size = NA, fill = "white", alpha = 0.7) +
    
    scale_color_viridis_d(option = "A", begin = 0.15, end = 0.85) +
    scale_fill_viridis_d(option = "A", begin = 0.15, end = 0.85) +
    coord_cartesian(xlim = c(0, x_limit), ylim = c(0, 1.05)) +
    theme_classic(base_size = 11) +
    labs(x = "DHW (°C-weeks)", y = "Mortality proportion",
         colour = NULL, fill = NULL, title = title_label,
         caption = caption_text)
  
  return(p)
}

p_hughes <- build_extrap_plot(
  hughes_fits,
  "A) Hughes et al. (2018) — Single event (2016)"
)

p_aims <- build_extrap_plot(
  aims_fits,
  "B) AIMS LTMP Manta Tow \u2014 Multi-event (1998\u20132024)",
  caption_text = "Note: AIMS mortality proportions are growth-corrected.\nReefs that gained cover are set to mortality \u2248 0.",
  x_limit = 25
)

Fig_Extrapolation <- p_hughes + p_aims +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

Fig_Extrapolation
ggsave("output/plots/Fig_ModelChoice_Extrapolation.png", Fig_Extrapolation,
       width = 12, height = 7, dpi = 300)


## ----model-extrapolation-table------------------------------------------------
#| label: model-extrapolation-table

# Extract predictions at key DHW values
extract_at_dhw <- function(fits, dataset_label, dhw_vals = c(8, 10, 12, 15, 20)) {
  g <- fits$grid %>% filter(round(DHW, 1) %in% dhw_vals)
  
  bind_rows(
    g %>% transmute(Dataset = dataset_label, DHW, Model = "Bounded GAM", Predicted_Mortality = bnd_fit),
    g %>% transmute(Dataset = dataset_label, DHW, Model = "Beta regression", Predicted_Mortality = beta_fit),
    g %>% transmute(Dataset = dataset_label, DHW, Model = "Quasi-binomial GLM", Predicted_Mortality = glm_fit),
    g %>% transmute(Dataset = dataset_label, DHW, Model = "ReefMod (sensitive)", Predicted_Mortality = reefmod_sens)
  )
}

extrap_table <- bind_rows(
  extract_at_dhw(hughes_fits, "Hughes et al. 2018"),
  extract_at_dhw(aims_fits, "AIMS LTMP Manta Tow")
) %>%
  tidyr::pivot_wider(names_from = DHW, values_from = Predicted_Mortality,
                     names_prefix = "DHW_") %>%
  mutate(across(starts_with("DHW_"), ~ round(., 3)))

datatable(extrap_table,
          caption = "Predicted mortality proportion at key DHW thresholds.",
          options = list(pageLength = 12, dom = "t"))


## ----benthic-data-------------------------------------------------------------
#| label: benthic-data
#| cache: true

load("Data/aims_ltmp/aims_ltmp.RData")
df.AIMS.full <- reef_photo_df

df.AIMS.Bent <- df.AIMS.full |> 
  filter(data_type == "photo-transect", domain_category == "reef",
         project_code %in% c("LTMP", "MMP"), purpose == "GROUP_LEVEL", variable == "HARD CORAL") |> 
  mutate(
    Reef_Name = toupper(domain_name),
    Reef_Name_Clean = gsub(" REEF(S)?| ISLAND| IS", "", Reef_Name)
  ) |>
  select(Reef_Name, Reef_Name_Clean, project_code, report_year, date, depth, lower:median, mean, reef_zone, id, shelf, reefpage_category)

# Create a cleaned version of the reference for matching, ensuring uniqueness
dat.Ref.Clean <- dat.AIMSRef %>%
  mutate(AIMS_REEF_NAME_Clean = gsub(" REEF(S)?| ISLAND| IS", "", toupper(Reef_Name))) %>%
  group_by(AIMS_REEF_NAME_Clean) %>%
  slice(1) %>% # Take the first match to avoid many-to-many explosion
  ungroup()

df.BENT.HC <- df.AIMS.Bent |> 
  left_join(select(dat.Ref.Clean, AIMS_REEF_NAME_Clean, ReefID, ReefName, SECT_NAME, SECTOR), 
            by = c("Reef_Name_Clean" = "AIMS_REEF_NAME_Clean"))

# Find max DHW for each reef year
dat.yrs.bent.hc <- df.BENT.HC |> ungroup() |> 
  distinct(ReefName, date, report_year) |>
  filter(!is.na(ReefName)) |> 
  mutate(
    year_int = floor(as.numeric(date)),
    frac = as.numeric(date) - year_int,
    Date = ymd(paste0(year_int, "-01-01")) +
      round(frac * as.numeric(ymd(paste0(year_int + 1, "-01-01")) - ymd(paste0(year_int, "-01-01")))),
    Prev = report_year - 1,
    Year = year(Date), Month = month(Date),
    DHWYear = ifelse(Month >= 7, Year, Year - 1)
  ) |> 
  select(-year_int, -frac) |> 
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
    idx = which(dat.gridSum$REEF_NAME == ReefName)[1],
    Pred.Grth = doCoralGrowth.Simp(CoralCover = Prev.HC * 100,
      B0 = dat.gridSum$B0[idx],
      WQ = dat.gridSum$WQ[idx], 
      HC.asym = dat.gridSum$HC.asym[idx])[3] / 100
  ) |> 
  mutate(
    Pred.Grth2 = ifelse(Lag == 2, 
      Pred.Grth + doCoralGrowth.Simp(CoralCover = (Prev.HC + Pred.Grth) * 100,
        B0 = dat.gridSum$B0[idx],
        WQ = dat.gridSum$WQ[idx], 
        HC.asym = dat.gridSum$HC.asym[idx])[3] / 100,
      Pred.Grth),
    Abs.Change = Change.Cor - Pred.Grth2,
    Rel.Change = Abs.Change / (Prev.HC + Pred.Grth2 + 0.001) # Small epsilon to avoid division by zero
  )

dat.mod.bent.hc <- df.BENT.HC %>%
  filter(Lag < 3) %>% 
  left_join(select(Disturbances_manta, ReefName, report_year, DISTURBANCE_TYPE), 
            by = c("ReefName", "report_year")) %>%
  filter(
    is.na(DISTURBANCE_TYPE) | !DISTURBANCE_TYPE %in% c("c", "d", "s", "f"),
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
  facet_wrap(~project_code) +
  theme_bw() +
  labs(x = "DHW (°C-weeks)", y = "Relative change (growth-adjusted)",
       title = "Benthic Transect: LTMP (Offshore) vs MMP (Inshore)")


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

load("data/Cheungetal2025/01_sstvar_blchrf.RData")

# Target predictors from Cheung et al.
predictor_vars <- c(
  # "DHW",                                              # Removed as it's redundant with MaxDHW.mean
  "mcur_90",                                           # current speed (90-day)
  "cloudp_90",                                         # cloud fraction (90-day)
  "secc3m",                                            # Secchi depth (3-month)
  "cbclus2",                                           # SST trajectory cluster (representative)
  "histmDHW6", "yrsince6",                             # bleaching history (>6 DHW)
  "histmDHW4", "yrsince4",                             # bleaching history (>4 DHW)
  # "ann_maxdhw_mx", "ann_maxdhw",                     # Removed: redundant with MaxDHW.mean
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
if (max(cheung_predictors$year) < 2024) {
  cat("WARNING: Cheung predictors only go up to", max(cheung_predictors$year), 
      ". Observations from 2024/2025 will be excluded from ML models.\n")
}
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
  filter(!is.na(mcur_90)) %>%
  mutate(survey_depth = 7.5, survey_type = "manta") # Assume manta tow is 6-9m

# 2. Benthic LTMP data
dat.ml.ltmp <- dat.mod.bent.hc %>%
  ungroup() %>%
  filter(project_code == "LTMP") %>%
  filter(!is.na(ReefID), !is.na(DHWYear)) %>%
  mutate(DHWYear = as.numeric(DHWYear)) %>%
  left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>%
  filter(!is.na(mcur_90)) %>%
  mutate(survey_type = "ltmp_benthic")

# 3. MMP (Inshore) data
dat.ml.mmp <- dat.mod.bent.hc %>%
  ungroup() %>%
  filter(project_code == "MMP") %>%
  filter(!is.na(ReefID), !is.na(DHWYear)) %>%
  mutate(DHWYear = as.numeric(DHWYear)) %>%
  left_join(cheung_predictors, by = c("ReefID" = "LABEL_id", "DHWYear" = "year")) %>%
  filter(!is.na(mcur_90)) %>%
  mutate(survey_type = "mmp_inshore")

# 4. Combined dataset
dat.ml.all <- bind_rows(dat.ml.mant, dat.ml.ltmp, dat.ml.mmp) %>%
  mutate(survey_type = as.factor(survey_type))

cat("Manta ML dataset:", nrow(dat.ml.mant), "observations\n")
cat("LTMP Benthic ML dataset:", nrow(dat.ml.ltmp), "observations\n")
cat("MMP Inshore ML dataset:", nrow(dat.ml.mmp), "observations\n")
cat("Combined ML dataset:", nrow(dat.ml.all), "observations\n")


## ----ml-join-summary----------------------------------------------------------
#| label: ml-join-summary

datatable(
  dat.ml.all %>%
    select(ReefName, DHWYear, survey_type, Abs.Change, Rel.Change,
           mcur_90, cloudp_90, histmDHW6) %>%
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
  drop_na() %>%
  select(where(~ var(.x, na.rm = TRUE) > 0))

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

library(ranger)
library(gbm)
library(tidymodels)

# Features to include in the models
features <- c("MaxDHW.mean", predictor_vars[predictor_vars %in% names(dat.ml.all)])

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
  
  # 3. GBM (with adaptive CV folds based on N)
  cv_folds <- if(nrow(df) > 200) 5 else if(nrow(df) > 50) 3 else 0
  
  gbm_r2 <- NA
  gbm_fit <- tryCatch({
    fit <- gbm(rf_formula, data = df, distribution = "gaussian", 
               n.trees = 500, interaction.depth = 3, shrinkage = 0.01, 
               cv.folds = cv_folds, n.cores = 1, bag.fraction = 0.5)
    
    # Extract R2 (using CV if available, otherwise training R2)
    best_iter <- if(cv_folds > 0) gbm.perf(fit, method = "cv", plot.it = FALSE) else 500
    gbm_pred <- predict(fit, df, n.trees = best_iter)
    
    # Calculate R2 safely using traditional formula to handle zero variance
    ss_res <- sum((df[[response_var]] - gbm_pred)^2, na.rm = TRUE)
    ss_tot <- sum((df[[response_var]] - mean(df[[response_var]], na.rm = TRUE))^2, na.rm = TRUE)
    gbm_r2 <- if(ss_tot > 0) 1 - (ss_res / ss_tot) else 0
    fit
  }, error = function(e) {
    message("GBM failed for ", dataset_name, ": ", e$message)
    NULL
  })
  
  message("Resulting GBM R2 for ", dataset_name, ": ", gbm_r2)
  
  # Return metrics
  metrics <- data.frame(
    Dataset = dataset_name,
    Response = response_var,
    Null_R2 = 0,
    RF_DHW_Only_R2 = rf_base$r.squared,
    RF_Multi_R2 = rf_fit$r.squared,
    GBM_R2 = gbm_r2
  )
  
  # Sanitize name for global assignment and filename
  safe_name <- tolower(dataset_name)
  safe_name <- gsub("[^a-z0-9]", "_", safe_name)
  safe_name <- gsub("_+", "_", safe_name)
  safe_name <- gsub("^_|_$", "", safe_name)
  
  # Store the best RF model in global env for later plotting
  mod_obj_name <- paste0("rf_", safe_name, "_", tolower(sub("\\.", "", response_var)))
  message("Fitting model: ", mod_obj_name, " (N=", nrow(df), ")")
  assign(mod_obj_name, rf_fit, envir = .GlobalEnv)
  
  # Save model objects
  if (!dir.exists("output/models")) dir.create("output/models", recursive = TRUE)
  saveRDS(rf_fit, paste0("output/models/", mod_obj_name, ".rds"))
  if (!is.null(gbm_fit)) {
    saveRDS(gbm_fit, paste0("output/models/", gsub("rf_", "gbm_", mod_obj_name), ".rds"))
  }
  
  return(metrics)
}

results <- bind_rows(
  fit_eval_models(dat.ml.mant, "Abs.Change", "Manta (All)"),
  fit_eval_models(dat.ml.ltmp, "Abs.Change", "LTMP Benthic"),
  fit_eval_models(dat.ml.mmp, "Abs.Change", "MMP Inshore"),
  fit_eval_models(dat.ml.all, "Abs.Change", "Combined (All)"),
  
  # Era-specific models for Manta Tow Data
  fit_eval_models(dat.ml.mant %>% filter(Era == "1998/2002"), "Abs.Change", "Manta (98/02)"),
  fit_eval_models(dat.ml.mant %>% filter(Era == "2016/2017"), "Abs.Change", "Manta (16/17)"),
  fit_eval_models(dat.ml.mant %>% filter(Era == "2020/2022"), "Abs.Change", "Manta (20/22)")
)

datatable(results %>% mutate(across(where(is.numeric), ~ round(., 3))),
          caption = "Machine Learning Model Performance (OOB/CV R²). 
          'Null' = mean only, 'DHW Only' = MaxDHW.mean predictor only, 'Multi' = all predictors.")


## ----ml-variable-importance---------------------------------------------------
#| label: ml-variable-importance
#| fig-cap: "Variable importance from Random Forest models (Manta vs Combined) for absolute change."
#| fig-width: 10
#| fig-height: 7

get_vi <- function(mod_name, label) {
  if (exists(mod_name)) {
    mod <- get(mod_name)
    data.frame(
      Variable = names(mod$variable.importance),
      Importance = mod$variable.importance,
      Source = label
    )
  } else {
    NULL
  }
}

vi_all <- bind_rows(
  get_vi("rf_manta_all_abschange", "Manta (All)"),
  get_vi("rf_ltmp_benthic_abschange", "LTMP Benthic"),
  get_vi("rf_mmp_inshore_abschange", "MMP Inshore"),
  get_vi("rf_manta_16_17_abschange", "Manta (16/17)"),
  get_vi("rf_manta_20_22_abschange", "Manta (20/22)")
)

if (nrow(vi_all) > 0) {
  # Order sources logically
  vi_all$Source <- factor(vi_all$Source, levels = c("Manta (All)", "LTMP Benthic", "MMP Inshore", 
                                                   "Manta (16/17)", "Manta (20/22)"))
  
  ggplot(vi_all, aes(x = reorder(Variable, Importance), y = Importance)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    facet_wrap(~Source, scales = "free_x", ncol = 2) +
    theme_bw() +
    labs(x = NULL, y = "Permutation Importance",
         title = "Variable Importance Comparison across Models and Eras")
}


## ----ml-partial-dependence----------------------------------------------------
#| label: ml-partial-dependence
#| fig-cap: "Partial dependence of absolute coral cover change on key predictors."
#| fig-height: 8
#| fig-width: 10

library(pdp)

if (exists("rf_combined_all_abschange")) {
  # Get top 4 predictors for Absolute Change
  vi <- rf_combined_all_abschange$variable.importance
  top_vars <- names(sort(vi, decreasing = TRUE))[1:4]
  
  df_pd <- dat.ml.all %>% select(Abs.Change, all_of(features)) %>% drop_na()
  
  pd_plots <- list()
  for (var in top_vars) {
    pd <- partial(rf_combined_all_abschange, pred.var = var, train = df_pd)
    p <- ggplot(pd, aes_string(x = var, y = "yhat")) +
      geom_line(size = 1, color = "darkred") +
      theme_bw() +
      labs(y = "Partial Dependence (Abs. Change)")
    pd_plots[[var]] <- p
  }
  
  library(patchwork)
  wrap_plots(pd_plots, ncol = 2)
}


## ----brms-pca-----------------------------------------------------------------
#| label: brms-pca
#| fig-cap: "PCA Biplot of environmental predictors from Cheung et al. (2025). The selected parsimonious variables are highlighted in red."
#| fig-width: 8
#| fig-height: 6

predictor_vars <- c("mcur_90", "cloudp_90", "secc3m", "cbclus2", "histmDHW6", 
                    "yrsince6", "histmDHW4", "yrsince4", "winyear_sd", "winyear_mean")

# Subset clean data for PCA
pca_data <- sstvar_blch2 %>%
  select(all_of(predictor_vars)) %>%
  drop_na()

pca_fit <- prcomp(pca_data, scale. = TRUE)

# Custom PCA biplot using ggplot
library(ggfortify)
autoplot(pca_fit, data = pca_data, loadings = TRUE, loadings.label = TRUE,
         loadings.colour = "grey50", loadings.label.colour = "black", alpha = 0.1) +
  theme_bw(base_size = 12) +
  labs(title = "PCA of Environmental Predictors",
       x = paste0("PC1 (", round(summary(pca_fit)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_fit)$importance[2, 2] * 100, 1), "%)"))


## ----brms-data-prep-----------------------------------------------------------
#| label: brms-data-prep

# 1. Select manta tow matched dataset with Cheung predictors
dat.ml.mant.clean <- dat.ml.mant %>%
  filter(!is.na(mcur_90))

N_obs <- nrow(dat.ml.mant.clean)

# Bounded mortality proportion with nudging
dat.ml.mant.clean <- dat.ml.mant.clean %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
    Mort.prop.nudge = (Mort.prop * (N_obs - 1) + 0.5) / N_obs
  )

# Gelman scaling (center and scale by 2 SD) — save params for back-transform
scale_2sd <- function(x) {
  (x - mean(x, na.rm = TRUE)) / (2 * sd(x, na.rm = TRUE))
}

# Store scaling parameters for back-transformation on plots
scale_params <- list(
  MaxDHW.mean = c(mean = mean(dat.ml.mant.clean$MaxDHW.mean, na.rm = TRUE),
                  sd2  = 2 * sd(dat.ml.mant.clean$MaxDHW.mean, na.rm = TRUE)),
  winyear_sd  = c(mean = mean(dat.ml.mant.clean$winyear_sd, na.rm = TRUE),
                  sd2  = 2 * sd(dat.ml.mant.clean$winyear_sd, na.rm = TRUE)),
  histmDHW6   = c(mean = mean(dat.ml.mant.clean$histmDHW6, na.rm = TRUE),
                  sd2  = 2 * sd(dat.ml.mant.clean$histmDHW6, na.rm = TRUE)),
  mcur_90     = c(mean = mean(dat.ml.mant.clean$mcur_90, na.rm = TRUE),
                  sd2  = 2 * sd(dat.ml.mant.clean$mcur_90, na.rm = TRUE)),
  cloudp_90   = c(mean = mean(dat.ml.mant.clean$cloudp_90, na.rm = TRUE),
                  sd2  = 2 * sd(dat.ml.mant.clean$cloudp_90, na.rm = TRUE)),
  secc3m      = c(mean = mean(dat.ml.mant.clean$secc3m, na.rm = TRUE),
                  sd2  = 2 * sd(dat.ml.mant.clean$secc3m, na.rm = TRUE))
)

dat.ml.mant.clean <- dat.ml.mant.clean %>%
  mutate(
    MaxDHW_s     = scale_2sd(MaxDHW.mean),
    winyear_sd_s = scale_2sd(winyear_sd),
    histmDHW6_s  = scale_2sd(histmDHW6),
    mcur_90_s    = scale_2sd(mcur_90),
    cloudp_90_s  = scale_2sd(cloudp_90),
    secc3m_s     = scale_2sd(secc3m)
  )

cat("Final scaled Manta Tow model dataset:", N_obs, "observations across", 
    length(unique(dat.ml.mant.clean$ReefName)), "unique reefs.\n")


## ----brms-modelling-----------------------------------------------------------
#| label: brms-modelling
#| cache: true

library(brms)
library(loo)
library(bayesplot)
library(tidybayes)

# Create models directory if not exists
if (!dir.exists("output/models")) dir.create("output/models", recursive = TRUE)

# 1. Fit Additive Baseline Model
cache_fit1 <- "output/models/brms_fit1_additive.rds"
if (file.exists(cache_fit1)) {
  fit1_brms <- readRDS(cache_fit1)
  cat("Loaded cached Additive Baseline Model.\n")
} else {
  cat("Fitting Additive Baseline Model...\n")
  fit1_brms <- brm(
    formula = Mort.prop.nudge ~ MaxDHW_s + winyear_sd_s + histmDHW6_s + mcur_90_s + cloudp_90_s + secc3m_s + (1 | Region / ReefName),
    data = dat.ml.mant.clean,
    family = Beta(link = "logit"),
    prior = c(
      prior(normal(-1, 1.5), class = "Intercept"),
      prior(normal(0.3, 0.3), class = "b", coef = "MaxDHW_s"),
      prior(normal(0, 1), class = "b")
    ),
    chains = 4,
    iter = 3000,
    warmup = 1500,
    cores = 4,
    seed = 42,
    control = list(adapt_delta = 0.95)
  )
  saveRDS(fit1_brms, cache_fit1)
}

# 2. Fit Interactive Model
cache_fit2 <- "output/models/brms_fit2_interactive.rds"
if (file.exists(cache_fit2)) {
  fit2_brms <- readRDS(cache_fit2)
  cat("Loaded cached Interactive Modulator Model.\n")
} else {
  cat("Fitting Interactive Modulator Model...\n")
  fit2_brms <- brm(
    formula = Mort.prop.nudge ~ MaxDHW_s * secc3m_s + MaxDHW_s * cloudp_90_s + winyear_sd_s + histmDHW6_s + mcur_90_s + (1 | Region / ReefName),
    data = dat.ml.mant.clean,
    family = Beta(link = "logit"),
    prior = c(
      prior(normal(-1, 1.5), class = "Intercept"),
      prior(normal(0.3, 0.3), class = "b", coef = "MaxDHW_s"),
      prior(normal(0, 0.5), class = "b", coef = "MaxDHW_s:secc3m_s"),
      prior(normal(0, 0.5), class = "b", coef = "MaxDHW_s:cloudp_90_s"),
      prior(normal(0, 1), class = "b")
    ),
    chains = 4,
    iter = 3000,
    warmup = 1500,
    cores = 4,
    seed = 42,
    control = list(adapt_delta = 0.95)
  )
  saveRDS(fit2_brms, cache_fit2)
}

# 3. Full Interaction Model — Beta
cache_fit3 <- "output/models/brms_fit3_full_interactive.rds"
if (file.exists(cache_fit3)) {
  fit3_brms <- readRDS(cache_fit3)
  cat("Loaded cached Full Interactive Beta Model.\n")
} else {
  cat("Fitting Full Interactive Beta Model...\n")
  fit3_brms <- brm(
    formula = Mort.prop.nudge ~ MaxDHW_s * secc3m_s + MaxDHW_s * cloudp_90_s +
              MaxDHW_s * histmDHW6_s + MaxDHW_s * mcur_90_s +
              winyear_sd_s + (1 | ReefName),
    data = dat.ml.mant.clean,
    family = Beta(link = "logit"),
    prior = c(
      prior(normal(-1, 1.5), class = "Intercept"),
      prior(normal(0.3, 0.3), class = "b", coef = "MaxDHW_s"),
      prior(normal(0, 1), class = "b"),
      prior(exponential(1), class = "sd")
    ),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.97)
  )
  saveRDS(fit3_brms, cache_fit3)
}

# 4. Full Interaction Model — Binomial-OLRE (Bayesian quasibinomial)
dat.ml.mant.clean$obs_id <- 1:nrow(dat.ml.mant.clean)
# Binomial requires integer responses: convert proportion to pseudo-counts
# Using 100 pseudo-trials gives good resolution for proportions
dat.ml.mant.clean$n_trials <- 100L
dat.ml.mant.clean$y_binom <- as.integer(round(dat.ml.mant.clean$Mort.prop.nudge * 100))

cache_fit3b <- "output/models/brms_fit3b_binomial_olre.rds"
if (file.exists(cache_fit3b)) {
  fit3b_brms <- readRDS(cache_fit3b)
  cat("Loaded cached Binomial-OLRE Model.\n")
} else {
  cat("Fitting Binomial-OLRE Model (Bayesian quasibinomial)...\n")
  fit3b_brms <- brm(
    formula = y_binom | trials(n_trials) ~ MaxDHW_s * secc3m_s + MaxDHW_s * cloudp_90_s +
              MaxDHW_s * histmDHW6_s + MaxDHW_s * mcur_90_s +
              winyear_sd_s + (1 | ReefName) + (1 | obs_id),
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
  saveRDS(fit3b_brms, cache_fit3b)
}

# Print summaries
cat("\n=== Model 3: Full Interactive Beta ===\n")
print(summary(fit3_brms))
cat("\n=== Model 3b: Binomial-OLRE (Bayesian Quasibinomial) ===\n")
print(summary(fit3b_brms))


## ----brms-loo-compare---------------------------------------------------------
#| label: brms-loo-compare

library(brms)
library(dplyr)
library(loo)
loo1 <- loo(fit1_brms, moment_match = FALSE)
loo2 <- loo(fit2_brms, moment_match = FALSE)
loo3 <- loo(fit3_brms, moment_match = FALSE)
loo_comp_beta <- loo_compare(loo1, loo2, loo3)

cat("=== Beta Models LOO-CV Comparison ===\n")
print(loo_comp_beta)

# Beta vs Binomial-OLRE (same fixed effects)
# Note: direct LOO comparison across families is informative but not strictly apples-to-apples
cat("\n=== Full Interactive Model: Beta vs Binomial-OLRE ===\n")
cat("Beta R²:\n")
print(bayes_R2(fit3_brms))
cat("\nBinomial-OLRE R²:\n")
print(bayes_R2(fit3b_brms))

# Present table
datatable(
  as.data.frame(loo_comp_beta) %>%
    mutate(across(where(is.numeric), ~ round(., 2))),
  caption = "LOO-CV Model Comparison — Beta Models (Additive → Interactive → Full Interactive)."
)


## ----brms-diagnostics---------------------------------------------------------
#| label: brms-diagnostics
#| fig-cap: "Posterior predictive check (y vs yrep) demonstrating how closely simulated datasets from the posterior match the observed mortality distribution."
#| fig-width: 8
#| fig-height: 5

library(brms)
library(bayesplot)
pp_check(fit2_brms, ndraws = 50) +
  theme_classic(base_size = 12) +
  labs(title = "Posterior Predictive Check (Interactive Model)",
       x = "Mortality Proportion", y = "Density")


## ----brms-marginal-effects----------------------------------------------------
#| label: brms-marginal-effects
#| fig-cap: "Bayesian conditional effects (Interactive Model) on original predictor scales."
#| fig-width: 10
#| fig-height: 8

# Helper: back-transform standardized axis to original scale
backtransform_ce <- function(ce_data, x_var, orig_var) {
  sp <- scale_params[[orig_var]]
  ce_data[[x_var]] <- ce_data[[x_var]] * sp["sd2"] + sp["mean"]
  ce_data
}

# Extract conditional effects from fit2_brms
ce_dhw    <- conditional_effects(fit2_brms, effects = "MaxDHW_s")[[1]] %>% backtransform_ce("MaxDHW_s", "MaxDHW.mean")
ce_hist   <- conditional_effects(fit2_brms, effects = "histmDHW6_s")[[1]] %>% backtransform_ce("histmDHW6_s", "histmDHW6")
ce_secchi <- conditional_effects(fit2_brms, effects = "secc3m_s")[[1]] %>% backtransform_ce("secc3m_s", "secc3m")
ce_cloud  <- conditional_effects(fit2_brms, effects = "cloudp_90_s")[[1]] %>% backtransform_ce("cloudp_90_s", "cloudp_90")

# Also back-transform the interaction surface
ce_int <- conditional_effects(fit2_brms, effects = "MaxDHW_s:secc3m_s")[[1]]
ce_int$MaxDHW_s <- ce_int$MaxDHW_s * scale_params$MaxDHW.mean["sd2"] + scale_params$MaxDHW.mean["mean"]
# secc3m_s is a conditioning factor — back-transform its label
if ("secc3m_s" %in% names(ce_int)) {
  ce_int$secc3m_s <- round(ce_int$secc3m_s * scale_params$secc3m["sd2"] + scale_params$secc3m["mean"], 1)
  ce_int$secc3m_s <- factor(paste0(ce_int$secc3m_s, " m"))
}

p_dhw <- ggplot(ce_dhw, aes(x = MaxDHW_s, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "steelblue", alpha = 0.25) +
  geom_line(color = "steelblue", linewidth = 1) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion") +
  theme_bw(base_size = 11)

p_hist <- ggplot(ce_hist, aes(x = histmDHW6_s, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "coral", alpha = 0.25) +
  geom_line(color = "coral", linewidth = 1) +
  labs(x = "Prior Severe Bleaching (count >6 DHW)", y = "Mortality Proportion") +
  theme_bw(base_size = 11)

p_secchi <- ggplot(ce_secchi, aes(x = secc3m_s, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "darkgreen", alpha = 0.25) +
  geom_line(color = "darkgreen", linewidth = 1) +
  labs(x = "Secchi Depth (m)", y = "Mortality Proportion") +
  theme_bw(base_size = 11)

p_interact <- ggplot(ce_int, aes(x = MaxDHW_s, y = estimate__, colour = secc3m_s, fill = secc3m_s)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.15) +
  geom_line(linewidth = 1) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion",
       colour = "Secchi Depth", fill = "Secchi Depth") +
  theme_bw(base_size = 11)

library(patchwork)
(p_dhw + p_hist) / (p_secchi + p_interact)


## ----brms-halfeye-manta-------------------------------------------------------
#| label: brms-halfeye-manta
#| fig-cap: "Manta Tow — Posterior coefficient distributions (Full Interactive Beta Model)."
#| fig-width: 10
#| fig-height: 6

library(brms)
library(dplyr)
library(tidybayes)

draws_manta <- as_draws_df(fit3_brms) %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "term", values_to = "estimate") %>%
  filter(term != "b_Intercept") %>%
  mutate(term = gsub("b_", "", term),
         term = gsub("MaxDHW_s", "DHW", term),
         term = gsub("histmDHW6_s", "Prior Exposure", term),
         term = gsub("secc3m_s", "Secchi Depth", term),
         term = gsub("cloudp_90_s", "Cloud Cover", term),
         term = gsub("winyear_sd_s", "Winter SST SD", term),
         term = gsub("mcur_90_s", "Current Speed", term),
         term = gsub(":", " × ", term))

ggplot(draws_manta, aes(x = estimate, y = reorder(term, estimate))) +
  stat_halfeye(.width = c(0.66, 0.95), point_interval = median_qi,
               fill = "steelblue", alpha = 0.7, normalize = "xy") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Posterior Estimate (logit scale)", y = NULL,
       title = "Manta Tow — Posterior Effect Distributions",
       subtitle = "Full Interactive Beta Model with (1|ReefName) random intercept") +
  theme_bw(base_size = 12)


## ----brms-dose-manta----------------------------------------------------------
#| label: brms-dose-manta
#| fig-cap: "Manta Tow — 4-panel Bayesian dose-response summary."
#| fig-width: 12
#| fig-height: 10

library(brms)
library(dplyr)
d <- dat.ml.mant.clean; m <- fit3_brms
dhw_ext <- max(d$MaxDHW.mean, na.rm = TRUE)
dhw_ext_s <- (dhw_ext - scale_params$MaxDHW.mean["mean"]) / scale_params$MaxDHW.mean["sd2"]
dhw_max_s <- (20 - scale_params$MaxDHW.mean["mean"]) / scale_params$MaxDHW.mean["sd2"]

# Helper: build prediction grid in scaled space, predict, back-transform DHW axis
brms_pred <- function(mod, nd, sp = scale_params) {
  pp <- posterior_epred(mod, newdata = nd, re_formula = NA)
  nd$pred    <- apply(pp, 2, median)
  nd$pred_lo <- apply(pp, 2, quantile, 0.05)
  nd$pred_hi <- apply(pp, 2, quantile, 0.95)
  nd$MaxDHW.mean <- nd$MaxDHW_s * sp$MaxDHW.mean["sd2"] + sp$MaxDHW.mean["mean"]
  nd
}

es <- annotate("rect", xmin = dhw_ext, xmax = 20, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.35)
el <- annotate("text", x = (dhw_ext + 20)/2, y = 0.05, label = "Extrapolation", fontface = "italic", colour = "grey40", size = 3)

# Panel A: DHW dose-response
ga <- data.frame(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 200),
  secc3m_s = 0, cloudp_90_s = 0, histmDHW6_s = 0, winyear_sd_s = 0, mcur_90_s = 0)
ga <- brms_pred(m, ga)
pA <- ggplot() + es + el +
  geom_point(data = d, aes(x = MaxDHW.mean, y = Mort.prop.nudge), alpha = 0.35, size = 1.2) +
  geom_ribbon(data = ga, aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = 0.2) +
  geom_line(data = ga, aes(x = MaxDHW.mean, y = pred), colour = "steelblue", linewidth = 1) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", title = "(A) Bayesian DHW Dose-Response") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11)

# Panel B: DHW × Prior Exposure
hl_raw <- quantile(d$histmDHW6, c(0.1, 0.5, 0.9), na.rm = TRUE)
hl_s <- (hl_raw - scale_params$histmDHW6["mean"]) / scale_params$histmDHW6["sd2"]
gb <- expand.grid(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 200), histmDHW6_s = hl_s) %>%
  mutate(secc3m_s = 0, cloudp_90_s = 0, winyear_sd_s = 0, mcur_90_s = 0)
gb <- brms_pred(m, gb) %>%
  mutate(histmDHW6 = histmDHW6_s * scale_params$histmDHW6["sd2"] + scale_params$histmDHW6["mean"],
         lbl = factor(paste0(round(histmDHW6, 1), " prior"), levels = paste0(round(hl_raw, 1), " prior")))
pB <- ggplot(gb, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Prior Exposure", fill = "Prior Exposure",
       title = "(B) DHW × Prior Exposure") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel C: DHW × Cloud Cover
cl_raw <- quantile(d$cloudp_90, c(0.1, 0.5, 0.9), na.rm = TRUE)
cl_s <- (cl_raw - scale_params$cloudp_90["mean"]) / scale_params$cloudp_90["sd2"]
gc <- expand.grid(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 200), cloudp_90_s = cl_s) %>%
  mutate(secc3m_s = 0, histmDHW6_s = 0, winyear_sd_s = 0, mcur_90_s = 0)
gc <- brms_pred(m, gc) %>%
  mutate(cloudp_90 = cloudp_90_s * scale_params$cloudp_90["sd2"] + scale_params$cloudp_90["mean"],
         lbl = factor(paste0(round(cloudp_90 * 100, 0), "%"), levels = paste0(round(cl_raw * 100, 0), "%")))
pC <- ggplot(gc, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Cloud Cover", fill = "Cloud Cover",
       title = "(C) DHW × Cloud Cover") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel D: Worst vs Best
cl_05 <- (quantile(d$cloudp_90, 0.05, na.rm=T) - scale_params$cloudp_90["mean"]) / scale_params$cloudp_90["sd2"]
cl_95 <- (quantile(d$cloudp_90, 0.95, na.rm=T) - scale_params$cloudp_90["mean"]) / scale_params$cloudp_90["sd2"]
h_95 <- (quantile(d$histmDHW6, 0.95, na.rm=T) - scale_params$histmDHW6["mean"]) / scale_params$histmDHW6["sd2"]
h_00 <- (0 - scale_params$histmDHW6["mean"]) / scale_params$histmDHW6["sd2"]
gd <- data.frame(MaxDHW_s = rep(seq(-1.5, dhw_max_s, length.out = 200), 2),
  scenario = rep(c("Worst case", "Best case"), each = 200),
  cloudp_90_s = rep(c(cl_05, cl_95), each = 200),
  histmDHW6_s = rep(c(h_00, h_95), each = 200),
  secc3m_s = 0, winyear_sd_s = 0, mcur_90_s = 0)
gd <- brms_pred(m, gd); gd$scenario <- factor(gd$scenario, levels = c("Worst case", "Best case"))
ge <- gd %>% select(MaxDHW.mean, scenario, pred, pred_lo, pred_hi) %>%
  pivot_wider(id_cols = MaxDHW.mean, names_from = scenario, values_from = c(pred, pred_lo, pred_hi))
pD <- ggplot() + es + el +
  geom_ribbon(data = ge, aes(x = MaxDHW.mean, ymin = `pred_Best case`, ymax = `pred_Worst case`), fill = "grey70", alpha = 0.25) +
  geom_ribbon(data = gd %>% filter(scenario == "Worst case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "firebrick", alpha = 0.12) +
  geom_ribbon(data = gd %>% filter(scenario == "Best case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "forestgreen", alpha = 0.12) +
  geom_line(data = gd, aes(x = MaxDHW.mean, y = pred, colour = scenario), linewidth = 1.2) +
  scale_colour_manual(values = c("Worst case" = "firebrick", "Best case" = "forestgreen")) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Scenario",
       title = "(D) Worst vs Best Case", subtitle = "Worst: low cloud, no prior exposure | Best: high cloud, high prior exposure") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85), plot.subtitle = element_text(size = 8))

(pA + pB) / (pC + pD)


## ----qbin-full-model----------------------------------------------------------
#| label: qbin-full-model

# ── 1. Fit the full interaction quasibinomial model ──
fit_qbin <- glm(
  Mort.prop.nudge ~ MaxDHW.mean * secc3m + MaxDHW.mean * cloudp_90 +
                    MaxDHW.mean * histmDHW6 + MaxDHW.mean * mcur_90 + winyear_sd,
  data = dat.ml.mant.clean,
  family = quasibinomial(link = "logit")
)

cat("Quasibinomial GLM summary:\n")
print(summary(fit_qbin))


## ----qbin-effect-sizes--------------------------------------------------------
#| label: qbin-effect-sizes
#| fig-cap: "Coefficient estimates (logit scale) from the quasibinomial GLM with 90% confidence intervals. The dashed line marks zero effect."
#| fig-width: 8
#| fig-height: 5

# Extract coefficients and SEs (drop intercept)
coef_df <- as.data.frame(coef(summary(fit_qbin)))
names(coef_df) <- c("Estimate", "SE", "t", "p")
coef_df$Term <- rownames(coef_df)
coef_df <- coef_df[coef_df$Term != "(Intercept)", ]

# 90% CI
z90 <- qnorm(0.95)
coef_df$lo <- coef_df$Estimate - z90 * coef_df$SE
coef_df$hi <- coef_df$Estimate + z90 * coef_df$SE

# Clean labels
coef_df$Term <- gsub("MaxDHW.mean", "DHW", coef_df$Term)
coef_df$Term <- gsub("histmDHW6", "Prior Exposure (>6 DHW)", coef_df$Term)
coef_df$Term <- gsub("secc3m", "Secchi Depth", coef_df$Term)
coef_df$Term <- gsub("cloudp_90", "Cloud Cover", coef_df$Term)
coef_df$Term <- gsub("winyear_sd", "Winter SST SD", coef_df$Term)
coef_df$Term <- gsub("mcur_90", "Current Speed", coef_df$Term)
coef_df$Term <- gsub(":", " × ", coef_df$Term)

# Significance flag
coef_df$sig <- ifelse(coef_df$lo > 0 | coef_df$hi < 0, "Significant", "Not significant")

# Order by effect size
coef_df$Term <- reorder(coef_df$Term, coef_df$Estimate)

ggplot(coef_df, aes(x = Estimate, y = Term, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.7, fatten = 3) +
  scale_colour_manual(values = c("Significant" = "black", "Not significant" = "grey60")) +
  labs(x = "Coefficient Estimate (logit scale)", y = NULL, colour = NULL,
       title = "Quasibinomial GLM — Effect Sizes (90% CI)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")


## ----qbin-dose-response-------------------------------------------------------
#| label: qbin-dose-response
#| fig-cap: "Quasibinomial GLM dose-response: (A) Beta vs Quasibinomial comparison, (B) DHW × Prior Exposure interaction, (C) DHW × Water Clarity interaction, (D) Worst-case vs Best-case mortality envelopes."
#| fig-width: 12
#| fig-height: 10

library(patchwork)

# Helper: predict with 90% CI on response scale via link-scale SE
predict_qbin_ci <- function(model, newdata, level = 0.90) {
  z <- qnorm(1 - (1 - level) / 2)
  p <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
  newdata$pred    <- plogis(p$fit)
  newdata$pred_lo <- plogis(p$fit - z * p$se.fit)
  newdata$pred_hi <- plogis(p$fit + z * p$se.fit)
  newdata
}

# Extrapolation threshold: max observed DHW in the data
dhw_extrap <- max(dat.ml.mant.clean$MaxDHW.mean, na.rm = TRUE)

# Reusable extrapolation shading layer
extrap_shade <- annotate("rect",
  xmin = dhw_extrap, xmax = 20, ymin = -Inf, ymax = Inf,
  fill = "grey80", alpha = 0.35
)
extrap_label <- annotate("text",
  x = (dhw_extrap + 20) / 2, y = 0.05,
  label = "Extrapolation", fontface = "italic", colour = "grey40", size = 3
)

# ── 2. Median-level DHW grid for Beta vs Qbin comparison (Panel A) ──
dhw_grid <- data.frame(
  MaxDHW.mean = seq(0, 20, length.out = 200),
  secc3m      = median(dat.ml.mant.clean$secc3m, na.rm = TRUE),
  cloudp_90   = median(dat.ml.mant.clean$cloudp_90, na.rm = TRUE),
  winyear_sd  = median(dat.ml.mant.clean$winyear_sd, na.rm = TRUE),
  histmDHW6   = median(dat.ml.mant.clean$histmDHW6, na.rm = TRUE),
  mcur_90     = median(dat.ml.mant.clean$mcur_90, na.rm = TRUE)
)
dhw_grid <- predict_qbin_ci(fit_qbin, dhw_grid)

# Beta predictions on same grid
dhw_grid_s <- dhw_grid %>%
  mutate(
    MaxDHW_s     = (MaxDHW.mean - scale_params$MaxDHW.mean["mean"]) / scale_params$MaxDHW.mean["sd2"],
    secc3m_s     = (secc3m      - scale_params$secc3m["mean"])      / scale_params$secc3m["sd2"],
    cloudp_90_s  = (cloudp_90   - scale_params$cloudp_90["mean"])   / scale_params$cloudp_90["sd2"],
    winyear_sd_s = (winyear_sd  - scale_params$winyear_sd["mean"])  / scale_params$winyear_sd["sd2"],
    histmDHW6_s  = (histmDHW6   - scale_params$histmDHW6["mean"])   / scale_params$histmDHW6["sd2"],
    mcur_90_s    = (mcur_90     - scale_params$mcur_90["mean"])     / scale_params$mcur_90["sd2"],
    Region   = dat.ml.mant.clean$Region[1],
    ReefName = dat.ml.mant.clean$ReefName[1]
  )
pred_beta <- fitted(fit2_brms, newdata = dhw_grid_s, re_formula = NA, summary = TRUE)
dhw_grid$pred_beta    <- pred_beta[, "Estimate"]
dhw_grid$pred_beta_lo <- pred_beta[, "Q2.5"]
dhw_grid$pred_beta_hi <- pred_beta[, "Q97.5"]

pA <- ggplot() +
  extrap_shade + extrap_label +
  geom_point(data = dat.ml.mant.clean, aes(x = MaxDHW.mean, y = Mort.prop.nudge),
             alpha = 0.35, size = 1.2) +
  geom_ribbon(data = dhw_grid, aes(x = MaxDHW.mean, ymin = pred_beta_lo, ymax = pred_beta_hi),
              fill = "steelblue", alpha = 0.15) +
  geom_ribbon(data = dhw_grid, aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi),
              fill = "firebrick", alpha = 0.12) +
  geom_line(data = dhw_grid, aes(x = MaxDHW.mean, y = pred_beta, colour = "Beta (brms)"), linewidth = 1) +
  geom_line(data = dhw_grid, aes(x = MaxDHW.mean, y = pred, colour = "Quasibinomial GLM"),
            linewidth = 1, linetype = "dashed") +
  scale_colour_manual(values = c("Beta (brms)" = "steelblue", "Quasibinomial GLM" = "firebrick")) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Model",
       title = "(A) Beta vs Quasibinomial") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.3, 0.85))

# ── 3. DHW × histmDHW6 interaction (Panel B) ──
hist_levels <- quantile(dat.ml.mant.clean$histmDHW6, c(0.1, 0.5, 0.9), na.rm = TRUE)
grid_hist <- expand.grid(
  MaxDHW.mean = seq(0, 20, length.out = 200),
  histmDHW6   = hist_levels
) %>% mutate(
  secc3m     = median(dat.ml.mant.clean$secc3m, na.rm = TRUE),
  cloudp_90  = median(dat.ml.mant.clean$cloudp_90, na.rm = TRUE),
  winyear_sd = median(dat.ml.mant.clean$winyear_sd, na.rm = TRUE),
  mcur_90    = median(dat.ml.mant.clean$mcur_90, na.rm = TRUE)
)
grid_hist <- predict_qbin_ci(fit_qbin, grid_hist) %>%
  mutate(hist_label = factor(paste0(round(histmDHW6, 1), " prior events"),
                             levels = paste0(round(hist_levels, 1), " prior events")))

pB <- ggplot(grid_hist, aes(x = MaxDHW.mean, y = pred, colour = hist_label, fill = hist_label)) +
  extrap_shade + extrap_label +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion",
       colour = "Prior Exposure\n(>6 DHW events)", fill = "Prior Exposure\n(>6 DHW events)",
       title = "(B) DHW × Prior Exposure") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.35, 0.85))

# ── 4. DHW × cloudp_90 interaction (Panel C) ──
cloud_levels <- quantile(dat.ml.mant.clean$cloudp_90, c(0.1, 0.5, 0.9), na.rm = TRUE)
grid_cloud <- expand.grid(
  MaxDHW.mean = seq(0, 20, length.out = 200),
  cloudp_90   = cloud_levels
) %>% mutate(
  secc3m     = median(dat.ml.mant.clean$secc3m, na.rm = TRUE),
  histmDHW6  = median(dat.ml.mant.clean$histmDHW6, na.rm = TRUE),
  winyear_sd = median(dat.ml.mant.clean$winyear_sd, na.rm = TRUE),
  mcur_90    = median(dat.ml.mant.clean$mcur_90, na.rm = TRUE)
)
grid_cloud <- predict_qbin_ci(fit_qbin, grid_cloud) %>%
  mutate(cloud_label = factor(paste0(round(cloudp_90 * 100, 0), "%"),
                              levels = paste0(round(cloud_levels * 100, 0), "%")))

pC <- ggplot(grid_cloud, aes(x = MaxDHW.mean, y = pred, colour = cloud_label, fill = cloud_label)) +
  extrap_shade + extrap_label +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion",
       colour = "Cloud Cover", fill = "Cloud Cover",
       title = "(C) DHW × Cloud Cover") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.35, 0.85))

# ── 5. Worst-case vs Best-case (Panel D) ──
# Worst: low cloud (5th %ile), zero prior exposure
# Best:  high cloud (95th %ile), high prior exposure (95th %ile)
grid_scenario <- data.frame(
  MaxDHW.mean = rep(seq(0, 20, length.out = 200), 2),
  scenario = rep(c("Worst case", "Best case"), each = 200),
  cloudp_90 = rep(c(quantile(dat.ml.mant.clean$cloudp_90, 0.05, na.rm = TRUE),
                    quantile(dat.ml.mant.clean$cloudp_90, 0.95, na.rm = TRUE)), each = 200),
  histmDHW6 = rep(c(0,
                    quantile(dat.ml.mant.clean$histmDHW6, 0.95, na.rm = TRUE)), each = 200),
  secc3m     = median(dat.ml.mant.clean$secc3m, na.rm = TRUE),
  winyear_sd = median(dat.ml.mant.clean$winyear_sd, na.rm = TRUE),
  mcur_90    = median(dat.ml.mant.clean$mcur_90, na.rm = TRUE)
)
grid_scenario <- predict_qbin_ci(fit_qbin, grid_scenario)
grid_scenario$scenario <- factor(grid_scenario$scenario, levels = c("Worst case", "Best case"))

grid_envelope <- grid_scenario %>%
  select(MaxDHW.mean, scenario, pred, pred_lo, pred_hi) %>%
  pivot_wider(id_cols = MaxDHW.mean, names_from = scenario, values_from = c(pred, pred_lo, pred_hi))

pD <- ggplot() +
  extrap_shade + extrap_label +
  geom_ribbon(data = grid_envelope,
              aes(x = MaxDHW.mean, ymin = `pred_Best case`, ymax = `pred_Worst case`),
              fill = "grey70", alpha = 0.25) +
  geom_ribbon(data = grid_scenario %>% filter(scenario == "Worst case"),
              aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "firebrick", alpha = 0.12) +
  geom_ribbon(data = grid_scenario %>% filter(scenario == "Best case"),
              aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "forestgreen", alpha = 0.12) +
  geom_line(data = grid_scenario, aes(x = MaxDHW.mean, y = pred, colour = scenario), linewidth = 1.2) +
  scale_colour_manual(values = c("Worst case" = "firebrick", "Best case" = "forestgreen")) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion",
       colour = "Scenario",
       title = "(D) Worst vs Best Case",
       subtitle = "Worst: low cloud, no prior exposure | Best: high cloud, high prior exposure") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 11) +
  theme(legend.position = c(0.35, 0.85),
        plot.subtitle = element_text(size = 8))

(pA + pB) / (pC + pD)


## ----qbin-ltmp-model----------------------------------------------------------
#| label: qbin-ltmp-model

# Compute proportion Acropora from composition-level benthic data
# df.AIMS.full is loaded in the benthic-data chunk
acrop_comp <- df.AIMS.full %>%
  filter(data_type == "photo-transect", domain_category == "reef",
         purpose == "COMPOSITION", variable == "HARD CORAL") %>%
  mutate(Reef_Name_Clean = gsub(" REEF(S)?| ISLAND| IS", "", toupper(domain_name))) %>%
  left_join(dat.Ref.Clean %>% select(AIMS_REEF_NAME_Clean, ReefName),
            by = c("Reef_Name_Clean" = "AIMS_REEF_NAME_Clean")) %>%
  filter(!is.na(ReefName)) %>%
  group_by(ReefName, report_year, depth) %>%
  summarise(
    total_hc = sum(mean, na.rm = TRUE),
    acrop_cover = sum(mean[grepl("Acropora", reefpage_category, ignore.case = TRUE)], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(prop_acropora = ifelse(total_hc > 0, acrop_cover / total_hc, 0))

cat("Acropora composition summary:\n")
cat("  Rows:", nrow(acrop_comp), "\n")
cat("  Mean prop_acropora:", round(mean(acrop_comp$prop_acropora, na.rm = TRUE), 3), "\n")

# Prepare LTMP benthic data with bounded mortality proportion
dat.ml.ltmp.clean <- dat.ml.ltmp %>%
  filter(!is.na(mcur_90), !is.na(Rel.Change)) %>%
  left_join(acrop_comp %>% select(ReefName, report_year, depth, prop_acropora),
            by = c("ReefName", "report_year", "depth")) %>%
  mutate(prop_acropora = replace_na(prop_acropora, median(prop_acropora, na.rm = TRUE)))

N_ltmp <- nrow(dat.ml.ltmp.clean)
dat.ml.ltmp.clean <- dat.ml.ltmp.clean %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
    Mort.prop.nudge = (Mort.prop * (N_ltmp - 1) + 0.5) / N_ltmp
  )

cat("LTMP Benthic model dataset:", N_ltmp, "observations\n")
cat("  prop_acropora range:", round(range(dat.ml.ltmp.clean$prop_acropora, na.rm = TRUE), 3), "\n")

fit_qbin_ltmp <- glm(
  Mort.prop.nudge ~ MaxDHW.mean * secc3m + MaxDHW.mean * cloudp_90 +
                    MaxDHW.mean * histmDHW6 + MaxDHW.mean * mcur_90 +
                    MaxDHW.mean * prop_acropora + winyear_sd,
  data = dat.ml.ltmp.clean,
  family = quasibinomial(link = "logit")
)
cat("LTMP Quasibinomial GLM summary (with Acropora proportion):\n")
print(summary(fit_qbin_ltmp))


## ----qbin-ltmp-effects--------------------------------------------------------
#| label: qbin-ltmp-effects
#| fig-cap: "LTMP Benthic — Coefficient estimates (logit scale) with 90% CI."
#| fig-width: 8
#| fig-height: 5

coef_ltmp <- as.data.frame(coef(summary(fit_qbin_ltmp)))
names(coef_ltmp) <- c("Estimate", "SE", "t", "p")
coef_ltmp$Term <- rownames(coef_ltmp)
coef_ltmp <- coef_ltmp[coef_ltmp$Term != "(Intercept)", ]
z90 <- qnorm(0.95)
coef_ltmp$lo <- coef_ltmp$Estimate - z90 * coef_ltmp$SE
coef_ltmp$hi <- coef_ltmp$Estimate + z90 * coef_ltmp$SE
coef_ltmp$Term <- gsub("MaxDHW.mean", "DHW", coef_ltmp$Term)
coef_ltmp$Term <- gsub("histmDHW6", "Prior Exposure (>6 DHW)", coef_ltmp$Term)
coef_ltmp$Term <- gsub("secc3m", "Secchi Depth", coef_ltmp$Term)
coef_ltmp$Term <- gsub("cloudp_90", "Cloud Cover", coef_ltmp$Term)
coef_ltmp$Term <- gsub("winyear_sd", "Winter SST SD", coef_ltmp$Term)
coef_ltmp$Term <- gsub("mcur_90", "Current Speed", coef_ltmp$Term)
coef_ltmp$Term <- gsub("prop_acropora", "Prop. Acropora", coef_ltmp$Term)
coef_ltmp$Term <- gsub(":", " × ", coef_ltmp$Term)
coef_ltmp$sig <- ifelse(coef_ltmp$lo > 0 | coef_ltmp$hi < 0, "Significant", "Not significant")
coef_ltmp$Term <- reorder(coef_ltmp$Term, coef_ltmp$Estimate)

ggplot(coef_ltmp, aes(x = Estimate, y = Term, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.7, fatten = 3) +
  scale_colour_manual(values = c("Significant" = "black", "Not significant" = "grey60")) +
  labs(x = "Coefficient Estimate (logit scale)", y = NULL, colour = NULL,
       title = "LTMP Benthic — Effect Sizes (90% CI)") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")


## ----qbin-ltmp-dose-----------------------------------------------------------
#| label: qbin-ltmp-dose
#| fig-cap: "LTMP Benthic — 6-panel dose-response summary."
#| fig-width: 12
#| fig-height: 14

library(patchwork)
pci <- function(mod, nd, level = 0.90) {
  z <- qnorm(1 - (1 - level)/2); p <- predict(mod, newdata = nd, type = "link", se.fit = TRUE)
  nd$pred <- plogis(p$fit); nd$pred_lo <- plogis(p$fit - z*p$se.fit); nd$pred_hi <- plogis(p$fit + z*p$se.fit); nd }
dhw_ext <- max(dat.ml.ltmp.clean$MaxDHW.mean, na.rm = TRUE)
es <- annotate("rect", xmin = dhw_ext, xmax = 20, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.35)
el <- annotate("text", x = (dhw_ext + 20)/2, y = 0.05, label = "Extrapolation", fontface = "italic", colour = "grey40", size = 3)
d <- dat.ml.ltmp.clean; m <- fit_qbin_ltmp

# Panel A: DHW dose-response
ga <- data.frame(MaxDHW.mean = seq(0, 20, length.out = 200), secc3m = median(d$secc3m, na.rm=T),
  cloudp_90 = median(d$cloudp_90, na.rm=T), histmDHW6 = median(d$histmDHW6, na.rm=T),
  winyear_sd = median(d$winyear_sd, na.rm=T), mcur_90 = median(d$mcur_90, na.rm=T),
  prop_acropora = median(d$prop_acropora, na.rm=T))
ga <- pci(m, ga)
pA2 <- ggplot() + es + el +
  geom_point(data = d, aes(x = MaxDHW.mean, y = Mort.prop.nudge), alpha = 0.35, size = 1.2) +
  geom_ribbon(data = ga, aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = 0.2) +
  geom_line(data = ga, aes(x = MaxDHW.mean, y = pred), colour = "steelblue", linewidth = 1) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", title = "(A) LTMP: DHW Dose-Response") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11)

# Panel B: Prior Exposure alone (main effect)
hl_range <- seq(quantile(d$histmDHW6, 0.05, na.rm = TRUE), quantile(d$histmDHW6, 0.95, na.rm = TRUE), length.out = 100)
gb_main <- data.frame(histmDHW6 = hl_range, MaxDHW.mean = median(d$MaxDHW.mean, na.rm = TRUE),
  secc3m = median(d$secc3m, na.rm=T), cloudp_90 = median(d$cloudp_90, na.rm=T),
  winyear_sd = median(d$winyear_sd, na.rm=T), mcur_90 = median(d$mcur_90, na.rm=T),
  prop_acropora = median(d$prop_acropora, na.rm=T))
gb_main <- pci(m, gb_main)
pB2_main <- ggplot(gb_main, aes(x = histmDHW6, y = pred)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "purple", alpha = 0.15) +
  geom_line(colour = "purple", linewidth = 1) +
  labs(x = "Prior Severe Bleaching Count", y = "Mortality Proportion", title = "(B) Prior Exposure Main Effect",
       subtitle = "DHW held at median") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11)

# Panel C: DHW × Prior Exposure
hl <- quantile(d$histmDHW6, c(0.1, 0.5, 0.9), na.rm = TRUE)
gb <- expand.grid(MaxDHW.mean = seq(0, 20, length.out = 200), histmDHW6 = hl) %>%
  mutate(secc3m = median(d$secc3m, na.rm=T), cloudp_90 = median(d$cloudp_90, na.rm=T),
         winyear_sd = median(d$winyear_sd, na.rm=T), mcur_90 = median(d$mcur_90, na.rm=T),
         prop_acropora = median(d$prop_acropora, na.rm=T))
gb <- pci(m, gb) %>% mutate(lbl = factor(paste0(round(histmDHW6,1)," prior"), levels = paste0(round(hl,1)," prior")))
pC2 <- ggplot(gb, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Prior Exposure", fill = "Prior Exposure",
       title = "(C) DHW × Prior Exposure") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel D: DHW × Cloud Cover
cl <- quantile(d$cloudp_90, c(0.1, 0.5, 0.9), na.rm = TRUE)
gc <- expand.grid(MaxDHW.mean = seq(0, 20, length.out = 200), cloudp_90 = cl) %>%
  mutate(secc3m = median(d$secc3m, na.rm=T), histmDHW6 = median(d$histmDHW6, na.rm=T),
         winyear_sd = median(d$winyear_sd, na.rm=T), mcur_90 = median(d$mcur_90, na.rm=T),
         prop_acropora = median(d$prop_acropora, na.rm=T))
gc <- pci(m, gc) %>% mutate(lbl = factor(paste0(round(cloudp_90*100,0),"%"), levels = paste0(round(cl*100,0),"%")))
pD2 <- ggplot(gc, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Cloud Cover", fill = "Cloud Cover",
       title = "(D) DHW × Cloud Cover") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel E: DHW × % Acropora
al_raw <- quantile(d$prop_acropora, c(0.1, 0.5, 0.9), na.rm = TRUE)
ge_acrop <- expand.grid(MaxDHW.mean = seq(0, 20, length.out = 200), prop_acropora = al_raw) %>%
  mutate(secc3m = median(d$secc3m, na.rm=T), cloudp_90 = median(d$cloudp_90, na.rm=T),
         histmDHW6 = median(d$histmDHW6, na.rm=T), winyear_sd = median(d$winyear_sd, na.rm=T),
         mcur_90 = median(d$mcur_90, na.rm=T))
ge_acrop <- pci(m, ge_acrop) %>%
  mutate(lbl = factor(paste0(round(prop_acropora * 100, 0), "% Acrop."), levels = paste0(round(al_raw * 100, 0), "% Acrop.")))
pE2 <- ggplot(ge_acrop, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Acropora", fill = "Acropora",
       title = "(E) DHW × % Acropora") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel F: Worst vs Best
cl_10 <- quantile(d$cloudp_90, 0.10, na.rm = TRUE)
cl_90 <- quantile(d$cloudp_90, 0.90, na.rm = TRUE)
h_10  <- quantile(d$histmDHW6, 0.10, na.rm = TRUE)
h_90  <- quantile(d$histmDHW6, 0.90, na.rm = TRUE)
a_10  <- quantile(d$prop_acropora, 0.10, na.rm = TRUE)
a_90  <- quantile(d$prop_acropora, 0.90, na.rm = TRUE)

gf <- data.frame(MaxDHW.mean = rep(seq(0, 20, length.out = 200), 2),
  scenario = rep(c("Worst case", "Best case"), each = 200),
  cloudp_90 = rep(c(cl_10, cl_90), each = 200),
  histmDHW6 = rep(c(h_10, h_90), each = 200),
  prop_acropora = rep(c(a_90, a_10), each = 200),
  secc3m = median(d$secc3m, na.rm=T), winyear_sd = median(d$winyear_sd, na.rm=T), mcur_90 = median(d$mcur_90, na.rm=T))
gf <- pci(m, gf); gf$scenario <- factor(gf$scenario, levels = c("Worst case", "Best case"))
gfe <- gf %>% select(MaxDHW.mean, scenario, pred, pred_lo, pred_hi) %>%
  pivot_wider(id_cols = MaxDHW.mean, names_from = scenario, values_from = c(pred, pred_lo, pred_hi))

pF2 <- ggplot() + es + el +
  geom_ribbon(data = gfe, aes(x = MaxDHW.mean, ymin = `pred_Best case`, ymax = `pred_Worst case`), fill = "grey70", alpha = 0.25) +
  geom_ribbon(data = gf %>% filter(scenario == "Worst case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "firebrick", alpha = 0.12) +
  geom_ribbon(data = gf %>% filter(scenario == "Best case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "forestgreen", alpha = 0.12) +
  geom_line(data = gf, aes(x = MaxDHW.mean, y = pred, colour = scenario), linewidth = 1.2) +
  scale_colour_manual(values = c("Worst case" = "firebrick", "Best case" = "forestgreen")) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Scenario",
       title = "(F) Worst vs Best Case",
       subtitle = "Worst: low cloud (10th), no prior (10th), high Acropora (90th) | Best: high cloud (90th), prior (90th), low Acropora (10th)") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85), plot.subtitle = element_text(size = 8))

(pA2 + pB2_main) / (pC2 + pD2) / (pE2 + pF2)


## ----brms-ltmp-model----------------------------------------------------------
#| label: brms-ltmp-model
#| cache: true

# Scale LTMP predictors (same Gelman scaling)
scale_params_ltmp <- list(
  MaxDHW.mean = c(mean = mean(dat.ml.ltmp.clean$MaxDHW.mean, na.rm=T), sd2 = 2*sd(dat.ml.ltmp.clean$MaxDHW.mean, na.rm=T)),
  histmDHW6   = c(mean = mean(dat.ml.ltmp.clean$histmDHW6, na.rm=T), sd2 = 2*sd(dat.ml.ltmp.clean$histmDHW6, na.rm=T)),
  secc3m      = c(mean = mean(dat.ml.ltmp.clean$secc3m, na.rm=T), sd2 = 2*sd(dat.ml.ltmp.clean$secc3m, na.rm=T)),
  cloudp_90   = c(mean = mean(dat.ml.ltmp.clean$cloudp_90, na.rm=T), sd2 = 2*sd(dat.ml.ltmp.clean$cloudp_90, na.rm=T)),
  winyear_sd  = c(mean = mean(dat.ml.ltmp.clean$winyear_sd, na.rm=T), sd2 = 2*sd(dat.ml.ltmp.clean$winyear_sd, na.rm=T)),
  mcur_90     = c(mean = mean(dat.ml.ltmp.clean$mcur_90, na.rm=T), sd2 = 2*sd(dat.ml.ltmp.clean$mcur_90, na.rm=T)),
  prop_acropora = c(mean = mean(dat.ml.ltmp.clean$prop_acropora, na.rm=T), sd2 = 2*sd(dat.ml.ltmp.clean$prop_acropora, na.rm=T))
)

s2 <- function(x) (x - mean(x, na.rm=T)) / (2*sd(x, na.rm=T))
dat.ml.ltmp.clean <- dat.ml.ltmp.clean %>%
  mutate(MaxDHW_s = s2(MaxDHW.mean), histmDHW6_s = s2(histmDHW6), secc3m_s = s2(secc3m),
         cloudp_90_s = s2(cloudp_90), winyear_sd_s = s2(winyear_sd), mcur_90_s = s2(mcur_90),
         prop_acropora_s = s2(prop_acropora), obs_id = 1:n(),
         n_trials = 100L, y_binom = as.integer(round(Mort.prop.nudge * 100)))

cat("LTMP brms dataset:", nrow(dat.ml.ltmp.clean), "obs,", length(unique(dat.ml.ltmp.clean$ReefName)), "reefs\n")

# Beta model
cache_ltmp_beta <- "output/models/brms_ltmp_beta.rds"
if (file.exists(cache_ltmp_beta)) {
  fit_brms_ltmp <- readRDS(cache_ltmp_beta)
  cat("Loaded cached LTMP Beta model.\n")
} else {
  cat("Fitting LTMP Beta model...\n")
  fit_brms_ltmp <- brm(
    Mort.prop.nudge ~ MaxDHW_s * cloudp_90_s + MaxDHW_s * histmDHW6_s +
                      MaxDHW_s * mcur_90_s + MaxDHW_s * prop_acropora_s +
                      secc3m_s + winyear_sd_s + (1 | ReefName),
    data = dat.ml.ltmp.clean, family = Beta(link = "logit"),
    prior = c(prior(normal(-1, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.97)
  )
  saveRDS(fit_brms_ltmp, cache_ltmp_beta)
}

# Binomial-OLRE model
cache_ltmp_binom <- "output/models/brms_ltmp_binomial_olre.rds"
if (file.exists(cache_ltmp_binom)) {
  fit_brmsb_ltmp <- readRDS(cache_ltmp_binom)
  cat("Loaded cached LTMP Binomial-OLRE model.\n")
} else {
  cat("Fitting LTMP Binomial-OLRE model...\n")
  fit_brmsb_ltmp <- brm(
    y_binom | trials(n_trials) ~ MaxDHW_s * cloudp_90_s + MaxDHW_s * histmDHW6_s +
                      MaxDHW_s * mcur_90_s + MaxDHW_s * prop_acropora_s +
                      secc3m_s + winyear_sd_s + (1 | ReefName) + (1 | obs_id),
    data = dat.ml.ltmp.clean, family = binomial(link = "logit"),
    prior = c(prior(normal(-1, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(1), class = "sd")),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.97)
  )
  saveRDS(fit_brmsb_ltmp, cache_ltmp_binom)
}

cat("\n=== LTMP Beta Model ===\n")
print(summary(fit_brms_ltmp))
cat("\nBeta R²:"); print(bayes_R2(fit_brms_ltmp))
cat("\nBinomial-OLRE R²:"); print(bayes_R2(fit_brmsb_ltmp))


## ----brms-halfeye-ltmp--------------------------------------------------------
#| label: brms-halfeye-ltmp
#| fig-cap: "LTMP Benthic — Posterior coefficient distributions."
#| fig-width: 10
#| fig-height: 6

library(brms)
library(dplyr)
library(tidybayes)
draws_ltmp <- as_draws_df(fit_brms_ltmp) %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "term", values_to = "estimate") %>%
  filter(term != "b_Intercept") %>%
  mutate(term = gsub("b_", "", term),
         term = gsub("MaxDHW_s", "DHW", term),
         term = gsub("histmDHW6_s", "Prior Exposure", term),
         term = gsub("secc3m_s", "Secchi Depth", term),
         term = gsub("cloudp_90_s", "Cloud Cover", term),
         term = gsub("winyear_sd_s", "Winter SST SD", term),
         term = gsub("mcur_90_s", "Current Speed", term),
         term = gsub("prop_acropora_s", "Prop. Acropora", term),
         term = gsub(":", " × ", term))

ggplot(draws_ltmp, aes(x = estimate, y = reorder(term, estimate))) +
  stat_halfeye(.width = c(0.66, 0.95), point_interval = median_qi,
               fill = "darkorange", alpha = 0.7, normalize = "xy") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Posterior Estimate (logit scale)", y = NULL,
       title = "LTMP Benthic — Posterior Effect Distributions",
       subtitle = "Beta Model with (1|ReefName) random intercept") +
  theme_bw(base_size = 12)


## ----brms-dose-ltmp-----------------------------------------------------------
#| label: brms-dose-ltmp
#| fig-cap: "LTMP Benthic — 6-panel Bayesian dose-response summary."
#| fig-width: 12
#| fig-height: 14

library(brms)
library(dplyr)
d <- dat.ml.ltmp.clean; m <- fit_brms_ltmp; sp <- scale_params_ltmp
dhw_ext <- max(d$MaxDHW.mean, na.rm = TRUE)
dhw_max_s <- (20 - sp$MaxDHW.mean["mean"]) / sp$MaxDHW.mean["sd2"]

brms_pred_ltmp <- function(mod, nd) {
  pp <- posterior_epred(mod, newdata = nd, re_formula = NA)
  nd$pred    <- apply(pp, 2, median)
  nd$pred_lo <- apply(pp, 2, quantile, 0.05)
  nd$pred_hi <- apply(pp, 2, quantile, 0.95)
  nd$MaxDHW.mean <- nd$MaxDHW_s * sp$MaxDHW.mean["sd2"] + sp$MaxDHW.mean["mean"]
  nd
}

es <- annotate("rect", xmin = dhw_ext, xmax = 20, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.35)
el <- annotate("text", x = (dhw_ext + 20)/2, y = 0.05, label = "Extrapolation", fontface = "italic", colour = "grey40", size = 3)

# Panel A: DHW
ga <- data.frame(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 100),
  secc3m_s = 0, cloudp_90_s = 0, histmDHW6_s = 0, winyear_sd_s = 0, mcur_90_s = 0, prop_acropora_s = 0)
ga <- brms_pred_ltmp(m, ga)
pA2 <- ggplot() + es + el +
  geom_point(data = d, aes(x = MaxDHW.mean, y = Mort.prop.nudge), alpha = 0.35, size = 1.2) +
  geom_ribbon(data = ga, aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = 0.2) +
  geom_line(data = ga, aes(x = MaxDHW.mean, y = pred), colour = "steelblue", linewidth = 1) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", title = "(A) Bayesian LTMP: DHW Dose-Response") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11)

# Panel B: Prior Exposure alone (main effect)
hl_range_raw <- seq(quantile(d$histmDHW6, 0.05, na.rm = TRUE), quantile(d$histmDHW6, 0.95, na.rm = TRUE), length.out = 100)
hl_range_s <- (hl_range_raw - sp$histmDHW6["mean"]) / sp$histmDHW6["sd2"]
gb_main <- data.frame(histmDHW6_s = hl_range_s, MaxDHW_s = 0,
  secc3m_s = 0, cloudp_90_s = 0, winyear_sd_s = 0, mcur_90_s = 0, prop_acropora_s = 0)
gb_main <- brms_pred_ltmp(m, gb_main) %>%
  mutate(histmDHW6 = histmDHW6_s * sp$histmDHW6["sd2"] + sp$histmDHW6["mean"])
pB2_main <- ggplot(gb_main, aes(x = histmDHW6, y = pred)) +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), fill = "purple", alpha = 0.15) +
  geom_line(colour = "purple", linewidth = 1) +
  labs(x = "Prior Severe Bleaching Count", y = "Mortality Proportion", title = "(B) Prior Exposure Main Effect",
       subtitle = "DHW held at median") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11)

# Panel C: DHW × Prior Exposure
hl_raw <- quantile(d$histmDHW6, c(0.1, 0.5, 0.9), na.rm = TRUE)
hl_s <- (hl_raw - sp$histmDHW6["mean"]) / sp$histmDHW6["sd2"]
gb <- expand.grid(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 100), histmDHW6_s = hl_s) %>%
  mutate(secc3m_s = 0, cloudp_90_s = 0, winyear_sd_s = 0, mcur_90_s = 0, prop_acropora_s = 0)
gb <- brms_pred_ltmp(m, gb) %>%
  mutate(histmDHW6 = histmDHW6_s * sp$histmDHW6["sd2"] + sp$histmDHW6["mean"],
         lbl = factor(paste0(round(histmDHW6, 1), " prior"), levels = paste0(round(hl_raw, 1), " prior")))
pC2 <- ggplot(gb, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Prior Exposure", fill = "Prior Exposure",
       title = "(C) DHW × Prior Exposure") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel D: DHW × Cloud Cover
cl_raw <- quantile(d$cloudp_90, c(0.1, 0.5, 0.9), na.rm = TRUE)
cl_s <- (cl_raw - sp$cloudp_90["mean"]) / sp$cloudp_90["sd2"]
gc <- expand.grid(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 100), cloudp_90_s = cl_s) %>%
  mutate(secc3m_s = 0, histmDHW6_s = 0, winyear_sd_s = 0, mcur_90_s = 0, prop_acropora_s = 0)
gc <- brms_pred_ltmp(m, gc) %>%
  mutate(cloudp_90 = cloudp_90_s * sp$cloudp_90["sd2"] + sp$cloudp_90["mean"],
         lbl = factor(paste0(round(cloudp_90 * 100, 0), "%"), levels = paste0(round(cl_raw * 100, 0), "%")))
pD2 <- ggplot(gc, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Cloud Cover", fill = "Cloud Cover",
       title = "(D) DHW × Cloud Cover") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel E: DHW × % Acropora
al_raw <- quantile(d$prop_acropora, c(0.1, 0.5, 0.9), na.rm = TRUE)
al_s <- (al_raw - sp$prop_acropora["mean"]) / sp$prop_acropora["sd2"]
ge_acrop <- expand.grid(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 100), prop_acropora_s = al_s) %>%
  mutate(secc3m_s = 0, cloudp_90_s = 0, histmDHW6_s = 0, winyear_sd_s = 0, mcur_90_s = 0)
ge_acrop <- brms_pred_ltmp(m, ge_acrop) %>%
  mutate(prop_acropora = prop_acropora_s * sp$prop_acropora["sd2"] + sp$prop_acropora["mean"],
         lbl = factor(paste0(round(prop_acropora*100, 0), "% Acrop."), levels = paste0(round(al_raw*100, 0), "% Acrop.")))
pE2 <- ggplot(ge_acrop, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Acropora", fill = "Acropora",
       title = "(E) DHW × % Acropora") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel F: Worst vs Best
cl_10 <- (quantile(d$cloudp_90, 0.10, na.rm=T) - sp$cloudp_90["mean"]) / sp$cloudp_90["sd2"]
cl_90 <- (quantile(d$cloudp_90, 0.90, na.rm=T) - sp$cloudp_90["mean"]) / sp$cloudp_90["sd2"]
h_10  <- (quantile(d$histmDHW6, 0.10, na.rm=T) - sp$histmDHW6["mean"]) / sp$histmDHW6["sd2"]
h_90  <- (quantile(d$histmDHW6, 0.90, na.rm=T) - sp$histmDHW6["mean"]) / sp$histmDHW6["sd2"]
a_10  <- (quantile(d$prop_acropora, 0.10, na.rm=T) - sp$prop_acropora["mean"]) / sp$prop_acropora["sd2"]
a_90  <- (quantile(d$prop_acropora, 0.90, na.rm=T) - sp$prop_acropora["mean"]) / sp$prop_acropora["sd2"]

gf <- data.frame(MaxDHW_s = rep(seq(-1.5, dhw_max_s, length.out = 100), 2),
  scenario = rep(c("Worst case", "Best case"), each = 100),
  cloudp_90_s = rep(c(cl_10, cl_90), each = 100),
  histmDHW6_s = rep(c(h_10, h_90), each = 100),
  prop_acropora_s = rep(c(a_90, a_10), each = 100),
  secc3m_s = 0, winyear_sd_s = 0, mcur_90_s = 0)
gf <- brms_pred_ltmp(m, gf); gf$scenario <- factor(gf$scenario, levels = c("Worst case", "Best case"))
gfe <- gf %>% select(MaxDHW.mean, scenario, pred, pred_lo, pred_hi) %>%
  pivot_wider(id_cols = MaxDHW.mean, names_from = scenario, values_from = c(pred, pred_lo, pred_hi))

pF2 <- ggplot() + es + el +
  geom_ribbon(data = gfe, aes(x = MaxDHW.mean, ymin = `pred_Best case`, ymax = `pred_Worst case`), fill = "grey70", alpha = 0.25) +
  geom_ribbon(data = gf %>% filter(scenario == "Worst case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "firebrick", alpha = 0.12) +
  geom_ribbon(data = gf %>% filter(scenario == "Best case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "forestgreen", alpha = 0.12) +
  geom_line(data = gf, aes(x = MaxDHW.mean, y = pred, colour = scenario), linewidth = 1.2) +
  scale_colour_manual(values = c("Worst case" = "firebrick", "Best case" = "forestgreen")) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Scenario",
       title = "(F) Worst vs Best Case",
       subtitle = "Worst: low cloud (10th), no prior (10th), high Acropora (90th) | Best: high cloud (90th), prior (90th), low Acropora (10th)") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85), plot.subtitle = element_text(size = 8))

(pA2 + pB2_main) / (pC2 + pD2) / (pE2 + pF2)


## ----qbin-mmp-model-----------------------------------------------------------
#| label: qbin-mmp-model

dat.ml.mmp.clean <- dat.ml.mmp %>%
  filter(!is.na(mcur_90), !is.na(Rel.Change)) %>%
  left_join(acrop_comp %>% select(ReefName, report_year, depth, prop_acropora),
            by = c("ReefName", "report_year", "depth")) %>%
  mutate(prop_acropora = replace_na(prop_acropora, median(prop_acropora, na.rm = TRUE)))

N_mmp <- nrow(dat.ml.mmp.clean)
dat.ml.mmp.clean <- dat.ml.mmp.clean %>%
  mutate(
    Mort.prop = pmin(pmax(-Rel.Change, 0), 1),
    Mort.prop.nudge = (Mort.prop * (N_mmp - 1) + 0.5) / N_mmp
  )

cat("MMP Inshore model dataset:", N_mmp, "observations\n")
cat("  prop_acropora range:", round(range(dat.ml.mmp.clean$prop_acropora, na.rm = TRUE), 3), "\n")
cat("  Depth levels:", sort(unique(dat.ml.mmp.clean$depth)), "\n")

fit_qbin_mmp <- glm(
  Mort.prop.nudge ~ MaxDHW.mean * secc3m + MaxDHW.mean * cloudp_90 +
                    MaxDHW.mean * histmDHW6 + MaxDHW.mean * mcur_90 +
                    MaxDHW.mean * prop_acropora + winyear_sd + depth,
  data = dat.ml.mmp.clean,
  family = quasibinomial(link = "logit")
)
cat("MMP Quasibinomial GLM summary (with Acropora proportion + depth):\n")
print(summary(fit_qbin_mmp))


## ----qbin-mmp-effects---------------------------------------------------------
#| label: qbin-mmp-effects
#| fig-cap: "MMP Inshore — Coefficient estimates (logit scale) with 90% CI."
#| fig-width: 8
#| fig-height: 5

coef_mmp <- as.data.frame(coef(summary(fit_qbin_mmp)))
names(coef_mmp) <- c("Estimate", "SE", "t", "p")
coef_mmp$Term <- rownames(coef_mmp)
coef_mmp <- coef_mmp[coef_mmp$Term != "(Intercept)", ]
coef_mmp$lo <- coef_mmp$Estimate - z90 * coef_mmp$SE
coef_mmp$hi <- coef_mmp$Estimate + z90 * coef_mmp$SE
coef_mmp$Term <- gsub("MaxDHW.mean", "DHW", coef_mmp$Term)
coef_mmp$Term <- gsub("histmDHW6", "Prior Exposure (>6 DHW)", coef_mmp$Term)
coef_mmp$Term <- gsub("secc3m", "Secchi Depth", coef_mmp$Term)
coef_mmp$Term <- gsub("cloudp_90", "Cloud Cover", coef_mmp$Term)
coef_mmp$Term <- gsub("winyear_sd", "Winter SST SD", coef_mmp$Term)
coef_mmp$Term <- gsub("mcur_90", "Current Speed", coef_mmp$Term)
coef_mmp$Term <- gsub("prop_acropora", "Prop. Acropora", coef_mmp$Term)
coef_mmp$Term <- gsub("^depth$", "Survey Depth", coef_mmp$Term)
coef_mmp$Term <- gsub(":", " × ", coef_mmp$Term)
coef_mmp$sig <- ifelse(coef_mmp$lo > 0 | coef_mmp$hi < 0, "Significant", "Not significant")
coef_mmp$Term <- reorder(coef_mmp$Term, coef_mmp$Estimate)

ggplot(coef_mmp, aes(x = Estimate, y = Term, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.7, fatten = 3) +
  scale_colour_manual(values = c("Significant" = "black", "Not significant" = "grey60")) +
  labs(x = "Coefficient Estimate (logit scale)", y = NULL, colour = NULL,
       title = "MMP Inshore — Effect Sizes (90% CI)") +
  theme_bw(base_size = 12) + theme(legend.position = "bottom")


## ----qbin-mmp-dose------------------------------------------------------------
#| label: qbin-mmp-dose
#| fig-cap: "MMP Inshore — 4-panel dose-response summary."
#| fig-width: 12
#| fig-height: 10

dhw_ext <- max(dat.ml.mmp.clean$MaxDHW.mean, na.rm = TRUE)
es <- annotate("rect", xmin = dhw_ext, xmax = 20, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.35)
el <- annotate("text", x = (dhw_ext + 20)/2, y = 0.05, label = "Extrapolation", fontface = "italic", colour = "grey40", size = 3)
d <- dat.ml.mmp.clean; m <- fit_qbin_mmp

# Panel A: DHW dose-response
ga <- data.frame(MaxDHW.mean = seq(0, 20, length.out = 200), secc3m = median(d$secc3m, na.rm=T),
  cloudp_90 = median(d$cloudp_90, na.rm=T), histmDHW6 = median(d$histmDHW6, na.rm=T),
  winyear_sd = median(d$winyear_sd, na.rm=T), mcur_90 = median(d$mcur_90, na.rm=T),
  prop_acropora = median(d$prop_acropora, na.rm=T), depth = median(d$depth, na.rm=T))
ga <- pci(m, ga)
pA3 <- ggplot() + es + el +
  geom_point(data = d, aes(x = MaxDHW.mean, y = Mort.prop.nudge), alpha = 0.35, size = 1.2) +
  geom_ribbon(data = ga, aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = 0.2) +
  geom_line(data = ga, aes(x = MaxDHW.mean, y = pred), colour = "steelblue", linewidth = 1) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", title = "(A) MMP: DHW Dose-Response") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11)

# Panel B: DHW × Prop Acropora
al_raw <- quantile(d$prop_acropora, c(0.1, 0.5, 0.9), na.rm = TRUE)
gb <- expand.grid(MaxDHW.mean = seq(0, 20, length.out = 200), prop_acropora = al_raw) %>%
  mutate(secc3m = median(d$secc3m, na.rm=T), cloudp_90 = median(d$cloudp_90, na.rm=T),
         histmDHW6 = median(d$histmDHW6, na.rm=T), winyear_sd = median(d$winyear_sd, na.rm=T),
         mcur_90 = median(d$mcur_90, na.rm=T), depth = median(d$depth, na.rm=T))
gb <- pci(m, gb) %>% mutate(lbl = factor(paste0(round(prop_acropora*100,0),"% Acrop."), levels = paste0(round(al_raw*100,0),"% Acrop.")))
pB3 <- ggplot(gb, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Acropora", fill = "Acropora",
       title = "(B) DHW × Prop. Acropora") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel C: DHW × Cloud Cover
cl <- quantile(d$cloudp_90, c(0.1, 0.5, 0.9), na.rm = TRUE)
gc <- expand.grid(MaxDHW.mean = seq(0, 20, length.out = 200), cloudp_90 = cl) %>%
  mutate(secc3m = median(d$secc3m, na.rm=T), histmDHW6 = median(d$histmDHW6, na.rm=T),
         winyear_sd = median(d$winyear_sd, na.rm=T), mcur_90 = median(d$mcur_90, na.rm=T),
         prop_acropora = median(d$prop_acropora, na.rm=T), depth = median(d$depth, na.rm=T))
gc <- pci(m, gc) %>% mutate(lbl = factor(paste0(round(cloudp_90*100,0),"%"), levels = paste0(round(cl*100,0),"%")))
pC3 <- ggplot(gc, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Cloud Cover", fill = "Cloud Cover",
       title = "(C) DHW × Cloud Cover") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel D: Worst vs Best
cl_10 <- quantile(d$cloudp_90, 0.10, na.rm = TRUE)
cl_90 <- quantile(d$cloudp_90, 0.90, na.rm = TRUE)
a_10  <- quantile(d$prop_acropora, 0.10, na.rm = TRUE)
a_90  <- quantile(d$prop_acropora, 0.90, na.rm = TRUE)

gd <- data.frame(MaxDHW.mean = rep(seq(0, 20, length.out = 200), 2),
  scenario = rep(c("Worst case", "Best case"), each = 200),
  cloudp_90 = rep(c(cl_10, cl_90), each = 200),
  prop_acropora = rep(c(a_90, a_10), each = 200),
  secc3m = median(d$secc3m, na.rm=T), histmDHW6 = median(d$histmDHW6, na.rm=T),
  winyear_sd = median(d$winyear_sd, na.rm=T), mcur_90 = median(d$mcur_90, na.rm=T),
  depth = median(d$depth, na.rm=T))
gd <- pci(m, gd); gd$scenario <- factor(gd$scenario, levels = c("Worst case", "Best case"))
ge <- gd %>% select(MaxDHW.mean, scenario, pred, pred_lo, pred_hi) %>%
  pivot_wider(id_cols = MaxDHW.mean, names_from = scenario, values_from = c(pred, pred_lo, pred_hi))

pD3 <- ggplot() + es + el +
  geom_ribbon(data = ge, aes(x = MaxDHW.mean, ymin = `pred_Best case`, ymax = `pred_Worst case`), fill = "grey70", alpha = 0.25) +
  geom_ribbon(data = gd %>% filter(scenario == "Worst case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "firebrick", alpha = 0.12) +
  geom_ribbon(data = gd %>% filter(scenario == "Best case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "forestgreen", alpha = 0.12) +
  geom_line(data = gd, aes(x = MaxDHW.mean, y = pred, colour = scenario), linewidth = 1.2) +
  scale_colour_manual(values = c("Worst case" = "firebrick", "Best case" = "forestgreen")) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Scenario",
       title = "(D) Worst vs Best Case",
       subtitle = "Worst: low cloud (10th), high Acropora (90th) | Best: high cloud (90th), low Acropora (10th)") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85), plot.subtitle = element_text(size = 8))

(pA3 + pB3) / (pC3 + pD3)


## ----brms-mmp-model-----------------------------------------------------------
#| label: brms-mmp-model
#| cache: true

# Scale MMP predictors
scale_params_mmp <- list(
  MaxDHW.mean = c(mean = mean(dat.ml.mmp.clean$MaxDHW.mean, na.rm=T), sd2 = 2*sd(dat.ml.mmp.clean$MaxDHW.mean, na.rm=T)),
  cloudp_90   = c(mean = mean(dat.ml.mmp.clean$cloudp_90, na.rm=T), sd2 = 2*sd(dat.ml.mmp.clean$cloudp_90, na.rm=T)),
  secc3m      = c(mean = mean(dat.ml.mmp.clean$secc3m, na.rm=T), sd2 = 2*sd(dat.ml.mmp.clean$secc3m, na.rm=T)),
  prop_acropora = c(mean = mean(dat.ml.mmp.clean$prop_acropora, na.rm=T), sd2 = 2*sd(dat.ml.mmp.clean$prop_acropora, na.rm=T))
)

s2 <- function(x) (x - mean(x, na.rm=T)) / (2*sd(x, na.rm=T))
dat.ml.mmp.clean <- dat.ml.mmp.clean %>%
  mutate(MaxDHW_s = s2(MaxDHW.mean), cloudp_90_s = s2(cloudp_90),
         secc3m_s = s2(secc3m), prop_acropora_s = s2(prop_acropora),
         obs_id = 1:n(),
         n_trials = 100L, y_binom = as.integer(round(Mort.prop.nudge * 100)))

cat("MMP brms dataset:", nrow(dat.ml.mmp.clean), "obs,", length(unique(dat.ml.mmp.clean$ReefName)), "reefs\n")
cat("Depth levels:", sort(unique(dat.ml.mmp.clean$depth)), "\n")

# Parsimonious Beta model
cache_mmp_beta <- "output/models/brms_mmp_beta.rds"
if (file.exists(cache_mmp_beta)) {
  fit_brms_mmp <- readRDS(cache_mmp_beta)
  cat("Loaded cached MMP Beta model.\n")
} else {
  cat("Fitting MMP Beta model (parsimonious)...\n")
  fit_brms_mmp <- brm(
    Mort.prop.nudge ~ MaxDHW_s * cloudp_90_s + MaxDHW_s * prop_acropora_s +
                      secc3m_s + depth + (1 | ReefName),
    data = dat.ml.mmp.clean, family = Beta(link = "logit"),
    prior = c(prior(normal(-1, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(2), class = "sd")),  # stronger regularization for small N
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.98)
  )
  saveRDS(fit_brms_mmp, cache_mmp_beta)
}

# Parsimonious Binomial-OLRE model
cache_mmp_binom <- "output/models/brms_mmp_binomial_olre.rds"
if (file.exists(cache_mmp_binom)) {
  fit_brmsb_mmp <- readRDS(cache_mmp_binom)
  cat("Loaded cached MMP Binomial-OLRE model.\n")
} else {
  cat("Fitting MMP Binomial-OLRE model (parsimonious)...\n")
  fit_brmsb_mmp <- brm(
    y_binom | trials(n_trials) ~ MaxDHW_s * cloudp_90_s + MaxDHW_s * prop_acropora_s +
                      secc3m_s + depth + (1 | ReefName) + (1 | obs_id),
    data = dat.ml.mmp.clean, family = binomial(link = "logit"),
    prior = c(prior(normal(-1, 1.5), class = "Intercept"),
              prior(normal(0, 1), class = "b"),
              prior(exponential(2), class = "sd")),
    chains = 4, iter = 4000, warmup = 2000, cores = 4, seed = 42,
    control = list(adapt_delta = 0.98)
  )
  saveRDS(fit_brmsb_mmp, cache_mmp_binom)
}

cat("\n=== MMP Beta Model (Parsimonious) ===\n")
print(summary(fit_brms_mmp))
cat("\nBeta R²:"); print(bayes_R2(fit_brms_mmp))
cat("\nBinomial-OLRE R²:"); print(bayes_R2(fit_brmsb_mmp))


## ----brms-halfeye-mmp---------------------------------------------------------
#| label: brms-halfeye-mmp
#| fig-cap: "MMP Inshore — Posterior coefficient distributions."
#| fig-width: 10
#| fig-height: 5

library(brms)
library(dplyr)
library(tidybayes)
draws_mmp <- as_draws_df(fit_brms_mmp) %>%
  select(starts_with("b_")) %>%
  pivot_longer(everything(), names_to = "term", values_to = "estimate") %>%
  filter(term != "b_Intercept") %>%
  mutate(term = gsub("b_", "", term),
         term = gsub("MaxDHW_s", "DHW", term),
         term = gsub("secc3m_s", "Secchi Depth", term),
         term = gsub("cloudp_90_s", "Cloud Cover", term),
         term = gsub("prop_acropora_s", "Prop. Acropora", term),
         term = gsub("^depth$", "Survey Depth", term),
         term = gsub(":", " × ", term))

ggplot(draws_mmp, aes(x = estimate, y = reorder(term, estimate))) +
  stat_halfeye(.width = c(0.66, 0.95), point_interval = median_qi,
               fill = "forestgreen", alpha = 0.7, normalize = "xy") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(x = "Posterior Estimate (logit scale)", y = NULL,
       title = "MMP Inshore — Posterior Effect Distributions",
       subtitle = "Parsimonious Beta Model with (1|ReefName) random intercept") +
  theme_bw(base_size = 12)


## ----brms-dose-mmp------------------------------------------------------------
#| label: brms-dose-mmp
#| fig-cap: "MMP Inshore — 4-panel Bayesian dose-response summary."
#| fig-width: 12
#| fig-height: 10

library(brms)
library(dplyr)
d <- dat.ml.mmp.clean; m <- fit_brms_mmp; sp <- scale_params_mmp
dhw_ext <- max(d$MaxDHW.mean, na.rm = TRUE)
dhw_max_s <- (20 - sp$MaxDHW.mean["mean"]) / sp$MaxDHW.mean["sd2"]

brms_pred_mmp <- function(mod, nd) {
  pp <- posterior_epred(mod, newdata = nd, re_formula = NA)
  nd$pred    <- apply(pp, 2, median)
  nd$pred_lo <- apply(pp, 2, quantile, 0.05)
  nd$pred_hi <- apply(pp, 2, quantile, 0.95)
  nd$MaxDHW.mean <- nd$MaxDHW_s * sp$MaxDHW.mean["sd2"] + sp$MaxDHW.mean["mean"]
  nd
}

es <- annotate("rect", xmin = dhw_ext, xmax = 20, ymin = -Inf, ymax = Inf, fill = "grey80", alpha = 0.35)
el <- annotate("text", x = (dhw_ext + 20)/2, y = 0.05, label = "Extrapolation", fontface = "italic", colour = "grey40", size = 3)

# Panel A: DHW
ga <- data.frame(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 100),
  secc3m_s = 0, cloudp_90_s = 0, prop_acropora_s = 0, depth = median(d$depth, na.rm=T))
ga <- brms_pred_mmp(m, ga)
pA3 <- ggplot() + es + el +
  geom_point(data = d, aes(x = MaxDHW.mean, y = Mort.prop.nudge), alpha = 0.35, size = 1.2) +
  geom_ribbon(data = ga, aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "steelblue", alpha = 0.2) +
  geom_line(data = ga, aes(x = MaxDHW.mean, y = pred), colour = "steelblue", linewidth = 1) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", title = "(A) Bayesian MMP: DHW Dose-Response") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11)

# Panel B: DHW × Prop Acropora
al_raw <- quantile(d$prop_acropora, c(0.1, 0.5, 0.9), na.rm = TRUE)
al_s <- (al_raw - sp$prop_acropora["mean"]) / sp$prop_acropora["sd2"]
gb <- expand.grid(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 100), prop_acropora_s = al_s) %>%
  mutate(secc3m_s = 0, cloudp_90_s = 0, depth = median(d$depth, na.rm=T))
gb <- brms_pred_mmp(m, gb) %>%
  mutate(prop_acropora = prop_acropora_s * sp$prop_acropora["sd2"] + sp$prop_acropora["mean"],
         lbl = factor(paste0(round(prop_acropora*100, 0), "% Acrop."), levels = paste0(round(al_raw*100, 0), "% Acrop.")))
pB3 <- ggplot(gb, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Set1") + scale_fill_brewer(palette = "Set1") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Acropora", fill = "Acropora",
       title = "(B) DHW × Prop. Acropora") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel C: DHW × Cloud Cover
cl_raw <- quantile(d$cloudp_90, c(0.1, 0.5, 0.9), na.rm = TRUE)
cl_s <- (cl_raw - sp$cloudp_90["mean"]) / sp$cloudp_90["sd2"]
gc <- expand.grid(MaxDHW_s = seq(-1.5, dhw_max_s, length.out = 100), cloudp_90_s = cl_s) %>%
  mutate(secc3m_s = 0, prop_acropora_s = 0, depth = median(d$depth, na.rm=T))
gc <- brms_pred_mmp(m, gc) %>%
  mutate(cloudp_90 = cloudp_90_s * sp$cloudp_90["sd2"] + sp$cloudp_90["mean"],
         lbl = factor(paste0(round(cloudp_90 * 100, 0), "%"), levels = paste0(round(cl_raw * 100, 0), "%")))
pC3 <- ggplot(gc, aes(x = MaxDHW.mean, y = pred, colour = lbl, fill = lbl)) + es + el +
  geom_ribbon(aes(ymin = pred_lo, ymax = pred_hi), alpha = 0.12, colour = NA) + geom_line(linewidth = 1) +
  scale_colour_brewer(palette = "Dark2") + scale_fill_brewer(palette = "Dark2") +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Cloud Cover", fill = "Cloud Cover",
       title = "(C) DHW × Cloud Cover") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85))

# Panel D: Worst vs Best
cl_10 <- (quantile(d$cloudp_90, 0.10, na.rm=T) - sp$cloudp_90["mean"]) / sp$cloudp_90["sd2"]
cl_90 <- (quantile(d$cloudp_90, 0.90, na.rm=T) - sp$cloudp_90["mean"]) / sp$cloudp_90["sd2"]
a_10 <- (quantile(d$prop_acropora, 0.10, na.rm=T) - sp$prop_acropora["mean"]) / sp$prop_acropora["sd2"]
a_90 <- (quantile(d$prop_acropora, 0.90, na.rm=T) - sp$prop_acropora["mean"]) / sp$prop_acropora["sd2"]
gd <- data.frame(MaxDHW_s = rep(seq(-1.5, dhw_max_s, length.out = 100), 2),
  scenario = rep(c("Worst case", "Best case"), each = 100),
  cloudp_90_s = rep(c(cl_10, cl_90), each = 100),
  prop_acropora_s = rep(c(a_90, a_10), each = 100),
  secc3m_s = 0, depth = median(d$depth, na.rm=T))
gd <- brms_pred_mmp(m, gd); gd$scenario <- factor(gd$scenario, levels = c("Worst case", "Best case"))
ge <- gd %>% select(MaxDHW.mean, scenario, pred, pred_lo, pred_hi) %>%
  pivot_wider(id_cols = MaxDHW.mean, names_from = scenario, values_from = c(pred, pred_lo, pred_hi))
pD3 <- ggplot() + es + el +
  geom_ribbon(data = ge, aes(x = MaxDHW.mean, ymin = `pred_Best case`, ymax = `pred_Worst case`), fill = "grey70", alpha = 0.25) +
  geom_ribbon(data = gd %>% filter(scenario == "Worst case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "firebrick", alpha = 0.12) +
  geom_ribbon(data = gd %>% filter(scenario == "Best case"), aes(x = MaxDHW.mean, ymin = pred_lo, ymax = pred_hi), fill = "forestgreen", alpha = 0.12) +
  geom_line(data = gd, aes(x = MaxDHW.mean, y = pred, colour = scenario), linewidth = 1.2) +
  scale_colour_manual(values = c("Worst case" = "firebrick", "Best case" = "forestgreen")) +
  labs(x = "Max DHW (°C-weeks)", y = "Mortality Proportion", colour = "Scenario",
       title = "(D) Worst vs Best Case", 
       subtitle = "Worst: low cloud (10th), high Acropora (90th) | Best: high cloud (90th), low Acropora (10th)") +
  coord_cartesian(ylim = c(0, 1)) + theme_bw(base_size = 11) + theme(legend.position = c(0.35, 0.85), plot.subtitle = element_text(size = 8))

(pA3 + pB3) / (pC3 + pD3)


## ----brms-r2-calc-------------------------------------------------------------
#| label: brms-r2-calc

library(brms)
library(dplyr)
library(knitr)

# Load brms models
fit3_brms   <- readRDS("output/models/brms_fit3_full_interactive.rds")
fit3b_brms  <- readRDS("output/models/brms_fit3b_binomial_olre.rds")
fit_brms_ltmp  <- readRDS("output/models/brms_ltmp_beta.rds")
fit_brmsb_ltmp <- readRDS("output/models/brms_ltmp_binomial_olre.rds")
fit_brms_mmp   <- readRDS("output/models/brms_mmp_beta.rds")
fit_brmsb_mmp  <- readRDS("output/models/brms_mmp_binomial_olre.rds")

# Helpers
get_brms_r2 <- function(mod) {
  r2 <- bayes_R2(mod)
  paste0(round(r2[1], 3), " (95% CrI: ", round(r2[3], 3), " - ", round(r2[4], 3), ")")
}

get_qbin_r2 <- function(fit) {
  paste0(round(1 - fit$deviance / fit$null.deviance, 3))
}

r2_df <- data.frame(
  `Dataset` = c("Manta Tow (Broadscale, 6–9 m)", "LTMP Benthic (Offshore, 5 m)", "MMP Inshore (Inshore, 2–5 m)"),
  `Bayesian Beta R²` = c(get_brms_r2(fit3_brms), get_brms_r2(fit_brms_ltmp), get_brms_r2(fit_brms_mmp)),
  `Bayesian Binomial-OLRE R²` = c(get_brms_r2(fit3b_brms), get_brms_r2(fit_brmsb_ltmp), get_brms_r2(fit_brmsb_mmp)),
  `Frequentist Quasibinomial McFadden R²` = c(get_qbin_r2(fit_qbin), get_qbin_r2(fit_qbin_ltmp), get_qbin_r2(fit_qbin_mmp)),
  check.names = FALSE
)

kable(r2_df, caption = "Variance Explained (R²) Comparison across Datasets and Model Families")

