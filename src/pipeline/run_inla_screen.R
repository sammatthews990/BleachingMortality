# Explicit entry point for the full spatial INLA candidate screen.

Sys.setenv(INLA_ST_RUN = '1')
source('src/models/fit_inla_spatiotemporal_screen.R')
