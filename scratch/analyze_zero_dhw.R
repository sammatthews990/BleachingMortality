# Analyze zero DHW for LTMP Benthic and MMP Inshore models with correct predictor sets
library(brms)
library(dplyr)

print("=== CHECKING LTMP BENTHIC MODELS ===")
m_ltmp_olre <- readRDS("output/models/brms_ltmp_binomial_olre.rds")
m_ltmp_beta <- readRDS("output/models/brms_ltmp_beta.rds")

d_ltmp <- m_ltmp_olre$data
print("LTMP Benthic Data MaxDHW_s range:")
print(range(d_ltmp$MaxDHW_s))
print("LTMP Benthic Data histmDHW6_s range:")
print(range(d_ltmp$histmDHW6_s))

# Predict at MaxDHW_s at minimum (which corresponds to 0 raw DHW)
nd_ltmp <- data.frame(
  MaxDHW_s = min(d_ltmp$MaxDHW_s),
  histmDHW6_s = c(min(d_ltmp$histmDHW6_s), median(d_ltmp$histmDHW6_s), max(d_ltmp$histmDHW6_s)),
  secc3m_s = 0,
  cloudp_90_s = 0,
  winyear_sd_s = 0,
  mcur_90_s = 0,
  prop_acropora_s = 0,
  n_trials = 100
)
pp_ltmp_olre <- posterior_epred(m_ltmp_olre, newdata = nd_ltmp, re_formula = NA)
nd_ltmp$pred_olre <- apply(pp_ltmp_olre, 2, median) / 100

nd_ltmp_beta <- data.frame(
  MaxDHW_s = min(d_ltmp$MaxDHW_s),
  histmDHW6_s = c(min(d_ltmp$histmDHW6_s), median(d_ltmp$histmDHW6_s), max(d_ltmp$histmDHW6_s)),
  secc3m_s = 0,
  cloudp_90_s = 0,
  winyear_sd_s = 0,
  mcur_90_s = 0,
  prop_acropora_s = 0
)
pp_ltmp_beta <- posterior_epred(m_ltmp_beta, newdata = nd_ltmp_beta, re_formula = NA)
nd_ltmp$pred_beta <- apply(pp_ltmp_beta, 2, median)

print("LTMP Benthic predictions at minimum MaxDHW_s (0 raw DHW):")
print(nd_ltmp %>% select(MaxDHW_s, histmDHW6_s, pred_olre, pred_beta))


print("=== CHECKING MMP INSHORE MODELS ===")
m_mmp_olre <- readRDS("output/models/brms_mmp_binomial_olre.rds")
m_mmp_beta <- readRDS("output/models/brms_mmp_beta.rds")

d_mmp <- m_mmp_olre$data
print("MMP Inshore Data MaxDHW_s range:")
print(range(d_mmp$MaxDHW_s))

nd_mmp <- data.frame(
  MaxDHW_s = min(d_mmp$MaxDHW_s),
  cloudp_90_s = 0,
  prop_acropora_s = c(min(d_mmp$prop_acropora_s), median(d_mmp$prop_acropora_s), max(d_mmp$prop_acropora_s)),
  secc3m_s = 0,
  depth = 5, # default depth
  n_trials = 100
)
pp_mmp_olre <- posterior_epred(m_mmp_olre, newdata = nd_mmp, re_formula = NA)
nd_mmp$pred_olre <- apply(pp_mmp_olre, 2, median) / 100

nd_mmp_beta <- data.frame(
  MaxDHW_s = min(d_mmp$MaxDHW_s),
  cloudp_90_s = 0,
  prop_acropora_s = c(min(d_mmp$prop_acropora_s), median(d_mmp$prop_acropora_s), max(d_mmp$prop_acropora_s)),
  secc3m_s = 0,
  depth = 5
)
pp_mmp_beta <- posterior_epred(m_mmp_beta, newdata = nd_mmp_beta, re_formula = NA)
nd_mmp$pred_beta <- apply(pp_mmp_beta, 2, median)

print("MMP Inshore predictions at minimum MaxDHW_s (0 raw DHW):")
print(nd_mmp %>% select(MaxDHW_s, prop_acropora_s, pred_olre, pred_beta))
