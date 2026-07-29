library(dplyr)
library(ggplot2)
library(patchwork)
library(DT)

# Load cached data
metrics_df <- readRDS("output/models/brt_manta_metrics.rds")
var_imp <- readRDS("output/models/brt_manta_var_imp.rds")
pdp_summary <- readRDS("output/models/brt_manta_pdp_boot.rds")
grid_dhw_sec <- readRDS("output/models/brt_manta_pdp_2way_sec.rds")
grid_dhw_hist <- readRDS("output/models/brt_manta_pdp_2way_hist.rds")

# 1. 1D Bootstrapped PDPs with 80% & 95% CIs
p_pdp_1d <- ggplot(pdp_summary, aes(x = x, y = median_yhat)) +
  geom_ribbon(aes(ymin = lo_95, ymax = hi_95, fill = "95% Bootstrap CI"), alpha = 0.2) +
  geom_ribbon(aes(ymin = lo_80, ymax = hi_80, fill = "80% Bootstrap CI"), alpha = 0.35) +
  geom_line(colour = "#004d4d", linewidth = 1.1) +
  scale_fill_manual(name = "Confidence Interval", values = c("95% Bootstrap CI" = "#008080", "80% Bootstrap CI" = "#008080")) +
  facet_wrap(~pred_label, scales = "free_x", ncol = 3) +
  labs(
    x = "Predictor Value",
    y = "Partial Dependence (Mortality Proportion)",
    title = "Boosted Regression Trees (BRT) — 1D Partial Dependence Plots",
    subtitle = "Marginal predictor effects with 80% (darker ribbon) and 95% (lighter ribbon) bootstrapped CIs (B = 100)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"), legend.position = "bottom")

ggsave("scratch/test_brt_pdp_1d_80_95.png", p_pdp_1d, width = 10, height = 6.5)

# 2. 2-Way Partial Dependence Interaction Heatmaps
p_2way_sec <- ggplot(grid_dhw_sec, aes(x = MaxDHW.mean, y = secc3m, z = pred)) +
  geom_tile(aes(fill = pred)) +
  geom_contour(colour = "white", alpha = 0.4, linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", name = "Predicted\nMortality") +
  labs(
    x = "Max DHW (°C-weeks)", y = "Secchi Depth (m)",
    title = "(A) Thermal Stress × Water Clarity",
    subtitle = "DHW × Secchi Depth Interaction Surface"
  ) +
  theme_bw(base_size = 11)

p_2way_hist <- ggplot(grid_dhw_hist, aes(x = MaxDHW.mean, y = histmDHW6, z = pred)) +
  geom_tile(aes(fill = pred)) +
  geom_contour(colour = "white", alpha = 0.4, linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", name = "Predicted\nMortality") +
  labs(
    x = "Max DHW (°C-weeks)", y = "Prior Stress Count (>6 DHW)",
    title = "(B) Thermal Stress × Prior Exposure",
    subtitle = "DHW × Prior Thermal Stress Interaction Surface"
  ) +
  theme_bw(base_size = 11)

p_2way_combined <- p_2way_sec + p_2way_hist + plot_layout(guides = "collect") +
  plot_annotation(
    title = "Boosted Regression Trees (BRT) — 2-Way Partial Dependence Interaction Surfaces",
    subtitle = "Non-linear 2D response surfaces capturing environmental modulation of bleaching mortality"
  )

ggsave("scratch/test_brt_pdp_2way.png", p_2way_combined, width = 11, height = 5)

cat("=== Test plots generated successfully! ===\n")
