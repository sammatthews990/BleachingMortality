library(dplyr)
library(gbm)
library(ggplot2)
library(patchwork)
library(purrr)

# Load cached BRT objects
brt_fit <- readRDS("output/models/brt_manta_fit.rds")
var_imp <- readRDS("output/models/brt_manta_var_imp.rds")
pdp_summary <- readRDS("output/models/brt_manta_pdp_boot.rds")

cat("=== 1. Variable Importance ===\n")
print(var_imp)

# Plot Variable Importance
p_vi <- ggplot(var_imp, aes(x = reorder(Variable_Label, rel.inf), y = rel.inf)) +
  geom_col(fill = "#1f77b4", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = sprintf("%.1f%%", rel.inf)), hjust = -0.15, size = 3.8, fontface = "bold") +
  coord_flip(ylim = c(0, max(var_imp$rel.inf) * 1.15)) +
  labs(
    x = NULL, y = "Relative Influence (%)",
    title = "Boosted Regression Trees (BRT) — Variable Importance",
    subtitle = "Relative contribution of environmental predictors to coral bleaching mortality (Manta Tow, N=110)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("scratch/test_brt_vi.png", p_vi, width = 8, height = 4.5)

# Plot Bootstrapped PDPs
p_pdp <- ggplot(pdp_summary, aes(x = x, y = median_yhat)) +
  geom_ribbon(aes(ymin = lo_95, ymax = hi_95), fill = "#008080", alpha = 0.2) +
  geom_ribbon(aes(ymin = lo_90, ymax = hi_90), fill = "#008080", alpha = 0.3) +
  geom_line(colour = "#004d4d", linewidth = 1.1) +
  facet_wrap(~pred_label, scales = "free_x", ncol = 3) +
  labs(
    x = "Predictor Value",
    y = "Partial Dependence (Mortality Proportion)",
    title = "Boosted Regression Trees (BRT) — Bootstrapped Partial Dependence Plots",
    subtitle = "Marginal effects with 90% and 95% bootstrapped confidence intervals (B = 100)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"))

ggsave("scratch/test_brt_pdp_final.png", p_pdp, width = 10, height = 6)
cat("=== Plots generated successfully! ===\n")
