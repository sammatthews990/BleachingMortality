library(dplyr)
library(ggplot2)

sindy_coefs <- readRDS("output/models/sindy_manta_coefs.rds")
grid_sindy <- readRDS("output/models/sindy_manta_dhw_grid.rds")
pdp_summary <- readRDS("output/models/brt_manta_pdp_boot.rds")
pdp_dhw <- pdp_summary %>% filter(predictor == "MaxDHW.mean")

cat("SINDy Discovered Coefficients:\n")
print(sindy_coefs)

p_comp <- ggplot() +
  geom_ribbon(data = pdp_dhw, aes(x = x, ymin = lo_95, ymax = hi_95), fill = "#008080", alpha = 0.2) +
  geom_line(data = pdp_dhw, aes(x = x, y = median_yhat, colour = "BRT Non-Parametric PDP"), linewidth = 1.2) +
  geom_line(data = grid_sindy, aes(x = DHW, y = pred, colour = "SINDy Parsimonious Law"), linewidth = 1.2, linetype = "dotdash") +
  scale_colour_manual(name = "Model Methodology", values = c("BRT Non-Parametric PDP" = "#004d4d", "SINDy Parsimonious Law" = "#e66101")) +
  labs(
    x = "Max DHW (°C-weeks)", y = "Predicted Mortality Proportion",
    title = "SINDy Law Discovery vs Non-Parametric BRT Response",
    subtitle = "Closed-form sparse polynomial equation vs machine learning boosted tree marginal curve"
  ) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("scratch/test_sindy_comp.png", p_comp, width = 8.5, height = 5.5)
cat("=== Test plot generated! ===\n")
