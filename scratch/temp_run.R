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

ggplot(data = newdata, aes(y = fit, x = DHW)) +
  geom_point(data = dat.bleach, aes(y = Mort.bin)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) +
  geom_line() +
  geom_line(data = newdata2, aes(y = fit, x = DHW), colour = "red") +
  theme_classic() +
  labs(x = "DHW (°C-weeks)", y = "Mortality proportion",
       caption = "Black = binomial GLM, Red = beta regression")

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

datatable(
  dat.ml.all %>%
    select(ReefName, DHWYear, survey_type, Abs.Change, Rel.Change,
           mcur_90, cloudp_90, histmDHW6) %>%
    slice_head(n = 50) %>%
    mutate(across(where(is.numeric), ~ round(., 3))),
  caption = "Sample of Combined ML dataset with Cheung predictors"
)

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

print('=== scale_params ===')
print(scale_params)

