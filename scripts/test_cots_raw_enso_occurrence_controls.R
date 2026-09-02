# Operational screens for uncompressed COTS pressure/compact timing,
# continuous La Nina state in mortality occurrence, and named reef controls.
# Historical disturbance labels supervise the COTS cause model only.
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(INLA); library(readr)
  library(stringr); library(tidyr)
})
Sys.setenv(INLA_ST_RUN = '0')
source('scripts/model_registry.R')
source('scripts/model_diagnostics.R')
root <- normalizePath('.', winslash = '/', mustWork = TRUE)
out_dir <- file.path(root, 'output', 'cots_raw_enso_occurrence')
fig_dir <- file.path(root, 'output', 'fig')
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

metric_summary_local <- function(rows) rows |> summarise(
  n = n(),
  rmse = sqrt(mean((observed_mortality - predicted_mortality)^2)),
  mae = mean(abs(observed_mortality - predicted_mortality)),
  predictive_r2 = 1 - sum((observed_mortality-predicted_mortality)^2) /
    sum((observed_mortality-mean(observed_mortality))^2),
  severe_rmse = sqrt(mean((observed_mortality[observed_mortality >= .5] -
    predicted_mortality[observed_mortality >= .5])^2)),
  occurrence_brier = mean((observed_occurrence-predicted_occurrence)^2),
  false_extreme_rate = mean(
    predicted_mortality[observed_mortality < .3] >= .3),
  .groups = 'drop')
max_or_na <- function(x) if (all(is.na(x))) NA_real_ else max(x,na.rm=TRUE)

selected <- read_csv(file.path(root, 'output',
  'cots_timing_cover_cyclone_soft_gate', 'cv_predictions.csv'),
  show_col_types = FALSE) |>
  filter(candidate == 'baseline_cots_logistic20_5_cyclone') |>
  mutate(source_observation_id=as.character(source_observation_id),
         fold=as.character(fold))
joint_context <- read_csv(file.path(root, 'output', 'explanatory_event_dhw',
  'event_dhw_brt_data.csv'), show_col_types=FALSE) |>
  mutate(source_observation_id=as.character(source_observation_id)) |>
  distinct(source_observation_id, .keep_all=TRUE)
reef_locations <- joint_context |> distinct(ReefID, lon, lat)
cots_hindcast <- read_csv(file.path(root, 'data', 'gbrPredsAdj_20262408.csv'),
  show_col_types=FALSE) |>
  transmute(ReefName=reefName, event_year=as.integer(year),
            cots_outbreak_probability=as.numeric(outbrProb))

# Only current and preceding years are used; no future interval pressure enters
# the compact candidate.
rrn_history <- read_csv(file.path(root, 'data', 'processed',
  'rrn_pressure_reef_year.csv'), show_col_types=FALSE) |>
  transmute(ReefID=LABEL_ID, event_year=as.integer(event_year),
            cots_density_year=coalesce(
              pmax(as.numeric(cot_idwmeanpertow), 0), 0)) |>
  arrange(ReefID, event_year) |> group_by(ReefID) |>
  mutate(cots_peak_to_date=cummax(coalesce(cots_density_year,0)),
    new_peak=cots_density_year >= cots_peak_to_date,
    peak_marker=if_else(new_peak,event_year,-Inf),
    cots_peak_year=cummax(peak_marker),
    cots_years_since_peak=pmin(event_year-cots_peak_year,10),
    cots_years_since_peak=if_else(is.finite(cots_years_since_peak),
                                  cots_years_since_peak,10),
    cots_density_excess_raw=pmax(cots_density_year-.22,0)) |> ungroup()

annual <- read_csv(file.path(root, 'data', 'processed',
  'annual_coral_transitions.csv'), show_col_types=FALSE) |>
  mutate(source_observation_id=as.character(source_observation_id)) |>
  left_join(cots_hindcast, by=c('ReefName','event_year'),
            relationship='many-to-one') |>
  left_join(rrn_history, by=c('ReefID','event_year'),
            relationship='many-to-one') |>
  left_join(reef_locations, by='ReefID', relationship='many-to-one') |>
  mutate(relative_loss=pmin(pmax(-cover_change_pp/(100*pmax(pre_cover,.02)),0),.999),
    absolute_cover_loss=pmin(pmax(-cover_change_pp/100,0),.999),
    cots_label=coalesce(disturbance_has_cots,FALSE) |
      str_detect(str_to_lower(coalesce(disturbance_text,'')),
                 'cots|crown-of-thorns'),
    cots_positive_relative=cots_label & relative_loss>.02,
    cots_positive_absolute=cots_label & absolute_cover_loss>.02,
    across(c(cots_outbreak_probability,cots_density_year,cots_peak_to_date,
             cots_density_excess_raw), ~coalesce(.x,0)),
    cots_years_since_peak=coalesce(cots_years_since_peak,10),
    acropora=prop_acropora_pre,
    is_manta=as.numeric(programme_key=='manta'),
    is_mmp=as.numeric(programme_key=='mmp'))

assessment_features <- selected |> distinct(source_observation_id,scheme,fold,
  .keep_all=TRUE) |>
  select(-any_of(c('cots_density_year','cots_peak_to_date',
                   'cots_density_excess_raw','cots_years_since_peak'))) |>
  left_join(rrn_history |> select(ReefID,event_year,cots_density_year,
      cots_peak_to_date,cots_density_excess_raw,cots_years_since_peak),
    by=c('ReefID','event_year'), relationship='many-to-one') |>
  mutate(across(c(cots_density_year,cots_peak_to_date,cots_density_excess_raw),
                ~coalesce(.x,0)),
    cots_years_since_peak=coalesce(cots_years_since_peak,10),
    cots_interval_excess_raw=pmax(coalesce(cots_interval_max,0)-.22,0),
    is_manta=as.numeric(programme_key=='manta'),
    is_mmp=as.numeric(programme_key=='mmp'))

prepare_features <- function(training,assessment,features) {
  for (feature in features) {
    tr <- as.numeric(training[[feature]]); av <- as.numeric(assessment[[feature]])
    replacement <- median(tr[is.finite(tr)],na.rm=TRUE)
    if (!is.finite(replacement)) replacement <- 0
    tr[!is.finite(tr)] <- replacement; av[!is.finite(av)] <- replacement
    centre <- mean(tr); spread <- sd(tr)
    if (!is.finite(spread) || spread<1e-8) spread <- 1
    training[[paste0(feature,'_z')]] <- (tr-centre)/spread
    assessment[[paste0(feature,'_z')]] <- (av-centre)/spread
  }
  list(training=training,assessment=assessment)
}

fit_cots_hurdle <- function(training,assessment,features,
                            loss_scale=c('relative','absolute')) {
  loss_scale <- match.arg(loss_scale)
  z <- prepare_features(training,assessment,features)
  training <- z$training; assessment <- z$assessment
  terms <- paste0(features,'_z')
  form <- as.formula(paste('response ~ 1 +',paste(terms,collapse=' + ')))
  pos <- paste0('cots_positive_',loss_scale)
  loss <- if (loss_scale=='relative') 'relative_loss' else 'absolute_cover_loss'
  od <- bind_rows(
    training |> transmute(response=as.numeric(.data[[pos]]),across(all_of(terms))),
    assessment |> transmute(response=NA_real_,across(all_of(terms))))
  of <- inla(form,family='binomial',data=od,
    control.predictor=list(compute=TRUE,link=1),verbose=FALSE)
  oi <- seq.int(nrow(training)+1L,nrow(od)); occurrence <- of$summary.fitted.values$mean[oi]
  pt <- training |> filter(.data[[pos]]) |>
    transmute(response=pmin(pmax(.data[[loss]],.001),.999),across(all_of(terms)))
  md <- bind_rows(pt,assessment |> transmute(response=NA_real_,across(all_of(terms))))
  mf <- inla(form,family='beta',data=md,
    control.predictor=list(compute=TRUE,link=1),verbose=FALSE)
  mi <- seq.int(nrow(pt)+1L,nrow(md)); magnitude <- mf$summary.fitted.values$mean[mi]
  relmag <- if (loss_scale=='absolute')
    pmin(magnitude/pmax(assessment$pre_cover,.02),.999) else magnitude
  list(prediction=occurrence*relmag,occurrence=occurrence,magnitude=relmag,
    fixed=bind_rows(
      as_tibble(of$summary.fixed,rownames='term') |> mutate(component='occurrence'),
      as_tibble(mf$summary.fixed,rownames='term') |>
        mutate(component=paste0(loss_scale,'_magnitude'))))
}

common_features <- c('cots_outbreak_probability','pre_cover','acropora',
  'lon','lat','is_manta','is_mmp')
cots_specs <- tribble(
  ~candidate,~pressure_source,~loss_scale,~include_timing,
  'cots_raw_interval_relative','interval','relative',FALSE,
  'cots_raw_current_timing_relative','current','relative',TRUE,
  'cots_raw_current_timing_absolute','current','absolute',TRUE)
cots_predictions <- tibble(); cots_fixed <- tibble()

for (scheme_name in unique(selected$scheme)) {
  for (fold_name in unique(selected$fold[selected$scheme==scheme_name])) {
    assessment <- assessment_features |>
      filter(scheme==scheme_name,fold==fold_name)
    if (scheme_name=='leave_one_event_out') {
      held <- as.integer(fold_name)
      training <- annual |>
        filter(!(baseline_report_year<=held & report_year>=held))
    } else {
      held_reefs <- unique(assessment$ReefID)
      training <- annual |> filter(!ReefID %in% held_reefs)
    }
    for (j in seq_len(nrow(cots_specs))) {
      spec <- cots_specs[j,]; tr <- training; av <- assessment
      tr$pressure_raw <- tr$cots_density_excess_raw
      av$pressure_raw <- if (spec$pressure_source=='interval')
        av$cots_interval_excess_raw else av$cots_density_excess_raw
      features <- c('pressure_raw',common_features)
      if (isTRUE(spec$include_timing)) features <-
        c('pressure_raw','cots_years_since_peak',common_features)
      fitted <- fit_cots_hurdle(tr,av,features,spec$loss_scale)
      activation <- if (spec$pressure_source=='interval')
        as.numeric(coalesce(av$cots_interval_max,0)>.22) else
        as.numeric(av$cots_density_year>.22 |
          (av$cots_peak_to_date>.22 & av$cots_years_since_peak<=2))
      cots_predictions <- bind_rows(cots_predictions,av |>
        mutate(cots_prediction_candidate=fitted$prediction,
          cots_occurrence_candidate=fitted$occurrence,
          cots_activation_candidate=activation,
          predicted_mortality=1-(1-thermal_prediction)*
            (1-cots_prediction_candidate*activation)*
            (1-cyclone_prediction*cyclone_activation),
          predicted_occurrence=1-(1-thermal_occurrence)*
            (1-cots_occurrence_candidate*activation)*
            (1-cyclone_occurrence*cyclone_activation),
          residual=observed_mortality-predicted_mortality,
          candidate=spec$candidate))
      cots_fixed <- bind_rows(cots_fixed,fitted$fixed |>
        mutate(candidate=spec$candidate,scheme=scheme_name,fold=fold_name))
    }
  }
}

all_cots_predictions <- bind_rows(
  selected |> mutate(candidate='current_log_interval_selected'),
  cots_predictions)
cots_comparison <- all_cots_predictions |> group_by(candidate,scheme) |>
  metric_summary_local() |> arrange(scheme,rmse)
cots_event_metrics <- all_cots_predictions |>
  filter(scheme=='leave_one_event_out') |> group_by(candidate,event_year) |>
  metric_summary_local()
base_event <- cots_comparison |> filter(
  candidate=='current_log_interval_selected',scheme=='leave_one_event_out')
base_reef_rmse <- cots_comparison |> filter(
  candidate=='current_log_interval_selected',scheme=='reef_blocked_5fold') |>
  pull(rmse)
eligible_cots <- cots_comparison |> filter(scheme=='leave_one_event_out') |>
  left_join(cots_comparison |> filter(scheme=='reef_blocked_5fold') |>
    select(candidate,reef_rmse=rmse),by='candidate') |>
  filter(rmse<base_event$rmse,severe_rmse<=base_event$severe_rmse+.005,
         reef_rmse<=base_reef_rmse+.002) |> arrange(rmse)
selected_cots <- if (nrow(eligible_cots)) eligible_cots$candidate[[1]] else
  'current_log_interval_selected'
write_csv(all_cots_predictions,file.path(out_dir,'cots_cv_predictions.csv'))
write_csv(cots_comparison,file.path(out_dir,'cots_model_comparison.csv'))
write_csv(cots_event_metrics,file.path(out_dir,'cots_event_metrics.csv'))
write_csv(cots_fixed,file.path(out_dir,'cots_fixed_effects.csv'))
write_lines(selected_cots,file.path(out_dir,'selected_cots_candidate.txt'))

# Full-data coefficient ledger for the promoted raw interval candidate. The
# annual training rows use same-year pressure; the prediction rows use the
# available interval maximum, matching the validation screen.
full_av <- assessment_features |> filter(scheme=='leave_one_event_out') |>
  distinct(source_observation_id,.keep_all=TRUE)
full_tr <- annual
full_tr$pressure_raw <- full_tr$cots_density_excess_raw
full_av$pressure_raw <- full_av$cots_interval_excess_raw
full_raw_cots <- fit_cots_hurdle(full_tr,full_av,
  c('pressure_raw',common_features),'relative')
old_causes <- read_csv(file.path(root,'output','cause_aware_competing_hazards',
  'cause_fixed_effects.csv'),show_col_types=FALSE)
selected_cause_fixed <- bind_rows(
  full_raw_cots$fixed |>
    mutate(component=paste0('cots_',component)),
  old_causes |> filter(str_detect(component,'^cyclone_')))
write_csv(selected_cause_fixed,
          file.path(out_dir,'selected_cause_fixed_effects.csv'))
compact_correlation <- annual |>
  select(cots_outbreak_probability,cots_density_excess_raw,
         cots_years_since_peak,pre_cover,acropora) |>
  mutate(across(everything(),as.numeric)) |> cor(use='pairwise.complete.obs',
  method='spearman') |> as.data.frame() |> as_tibble(rownames='feature')
write_csv(compact_correlation,file.path(out_dir,'cots_compact_correlations.csv'))
cots_focal <- all_cots_predictions |>
  filter(scheme=='leave_one_event_out',
    ReefID %in% c('21-556','18-065','14-126','20-104','16-015',
                  '17-034','18-112','18-120')) |>
  group_by(candidate,ReefID,ReefName,event_year) |>
  summarise(observed_mortality=mean(observed_mortality),
    predicted_mortality=mean(predicted_mortality,na.rm=TRUE),
    cots_density_year=max(cots_density_year,na.rm=TRUE),
    cots_interval_max=max(cots_interval_max,na.rm=TRUE),
    years_since_peak=max_or_na(cots_years_since_peak),.groups='drop')
write_csv(cots_focal,file.path(out_dir,'cots_focal_predictions.csv'))

# Initial-forecast event state. The WQC delta is a coloured-water/runoff proxy,
# not measured river discharge. RONI and SOI are alternatives, never cofit.
base_for_occurrence <- all_cots_predictions |> filter(candidate==selected_cots) |>
  select(-candidate) |>
  left_join(joint_context |> select(source_observation_id,survey_date,depth,
    observed_pre_cover,prop_acropora_pre,log_coastal_rain30,wqc_prior10_delta,
    wqc_prior10_percentile,wqc_freqcc12,cloudp_90,mcur_90,secc3m,ann_maxdhw,
    applied_dhw_uplift,cot_interval_idw_max,cots_outbreak_probability,
    disturbance_text),by='source_observation_id',relationship='many-to-one') |>
  left_join(read_csv(file.path(root,'data','processed','enso_event_context.csv'),
      show_col_types=FALSE) |>
    select(event_year,roni_bleaching_summer_mean,soi_dec_mar_mean),
    by='event_year',relationship='many-to-one') |>
  group_by(scheme,ReefID) |>
  mutate(rainfall_event_anomaly=log_coastal_rain30-
    median(log_coastal_rain30,na.rm=TRUE)) |> ungroup() |>
  mutate(runoff_proxy_anomaly=wqc_prior10_delta,
    base_conditional_magnitude=pmin(
      predicted_mortality/pmax(predicted_occurrence,.01),.999))

fit_occurrence_calibrator <- function(training,assessment,features) {
  z <- prepare_features(training,assessment,features)
  training <- z$training; assessment <- z$assessment
  terms <- paste0(features,'_z')
  form <- as.formula(paste('response ~ 1 +',paste(terms,collapse=' + '),
                           '+ offset(base_logit)'))
  combined <- bind_rows(
    training |> transmute(response=observed_occurrence,
      base_logit=qlogis(pmin(pmax(predicted_occurrence,.001),.999)),
      across(all_of(terms))),
    assessment |> transmute(response=NA_real_,
      base_logit=qlogis(pmin(pmax(predicted_occurrence,.001),.999)),
      across(all_of(terms))))
  fit <- inla(form,family='binomial',data=combined,
    control.fixed=list(mean=0,prec=4,mean.intercept=0,prec.intercept=1),
    control.predictor=list(compute=TRUE,link=1),verbose=FALSE)
  index <- seq.int(nrow(training)+1L,nrow(combined))
  list(occurrence=fit$summary.fitted.values$mean[index],
       fixed=as_tibble(fit$summary.fixed,rownames='term'))
}

occurrence_specs <- list(
  occurrence_offset_only=character(),
  occurrence_roni_rain_wqc=c('roni_bleaching_summer_mean',
    'rainfall_event_anomaly','runoff_proxy_anomaly'),
  occurrence_soi_rain_wqc=c('soi_dec_mar_mean','rainfall_event_anomaly',
    'runoff_proxy_anomaly'))
occurrence_predictions <- tibble(); occurrence_fixed <- tibble()
for (scheme_name in unique(base_for_occurrence$scheme)) {
  sr <- base_for_occurrence |> filter(scheme==scheme_name)
  for (fold_name in unique(sr$fold)) {
    av <- sr |> filter(fold==fold_name); tr <- sr |> filter(fold!=fold_name)
    for (name in names(occurrence_specs)) {
      features <- occurrence_specs[[name]]
      if (!length(features)) { pred <- av$predicted_occurrence; fixed <- tibble() }
      else { fit <- fit_occurrence_calibrator(tr,av,features); pred <- fit$occurrence
        fixed <- fit$fixed |> mutate(candidate=name,scheme=scheme_name,fold=fold_name) }
      occurrence_predictions <- bind_rows(occurrence_predictions,av |>
        mutate(predicted_occurrence_initial=pred,
          predicted_mortality_initial=pred*base_conditional_magnitude,
          residual_initial=observed_mortality-predicted_mortality_initial,
          occurrence_candidate=name))
      occurrence_fixed <- bind_rows(occurrence_fixed,fixed)
    }
  }
}
occurrence_scored <- occurrence_predictions |>
  mutate(predicted_occurrence=predicted_occurrence_initial,
         predicted_mortality=predicted_mortality_initial)
occurrence_comparison <- occurrence_scored |>
  group_by(occurrence_candidate,scheme) |> metric_summary_local() |>
  arrange(scheme,rmse)
occurrence_event_metrics <- occurrence_scored |>
  filter(scheme=='leave_one_event_out') |>
  group_by(occurrence_candidate,event_year) |>
  summarise(metric_summary_local(pick(everything())),
    observed_mean=mean(observed_mortality),
    predicted_mean=mean(predicted_mortality),
    observed_occurrence_mean=mean(observed_occurrence),
    predicted_occurrence_mean=mean(predicted_occurrence),.groups='drop')
offset_metric <- occurrence_comparison |> filter(
  occurrence_candidate=='occurrence_offset_only',scheme=='leave_one_event_out')
eligible_occurrence <- occurrence_comparison |>
  filter(scheme=='leave_one_event_out',
    occurrence_candidate!='occurrence_offset_only',rmse<offset_metric$rmse,
    severe_rmse<=offset_metric$severe_rmse+.005) |> arrange(rmse)
selected_occurrence <- if (nrow(eligible_occurrence))
  eligible_occurrence$occurrence_candidate[[1]] else 'occurrence_offset_only'
write_csv(occurrence_predictions,file.path(out_dir,'occurrence_cv_predictions.csv'))
write_csv(occurrence_comparison,file.path(out_dir,'occurrence_model_comparison.csv'))
write_csv(occurrence_event_metrics,file.path(out_dir,'occurrence_event_metrics.csv'))
write_csv(occurrence_fixed,file.path(out_dir,'occurrence_fixed_effects.csv'))
write_lines(selected_occurrence,file.path(out_dir,'selected_occurrence_candidate.txt'))

# Sequential nowcast: earliest 20% of surveyed reefs form a shrunk prevalence
# signal and are excluded from scoring, so the update does not predict itself.
loeo <- base_for_occurrence |> filter(scheme=='leave_one_event_out') |>
  mutate(survey_date=as.Date(survey_date))
early_lookup <- loeo |> group_by(event_year,ReefID) |>
  summarise(first_survey=min(survey_date,na.rm=TRUE),
    observed_occurrence=max(observed_occurrence),
    predicted_occurrence=mean(predicted_occurrence),.groups='drop') |>
  group_by(event_year) |> arrange(first_survey,ReefID,.by_group=TRUE) |>
  mutate(early_n=max(3L,ceiling(.2*n())),is_early_reef=row_number()<=early_n) |>
  ungroup()
early_summary <- early_lookup |> group_by(event_year) |>
  summarise(early_n=sum(is_early_reef),
    early_observed_prevalence=(sum(observed_occurrence[is_early_reef])+1)/
      (sum(is_early_reef)+2),
    early_predicted_prevalence=(sum(predicted_occurrence[is_early_reef])+1)/
      (sum(is_early_reef)+2),
    early_logit_delta=qlogis(early_observed_prevalence)-
      qlogis(early_predicted_prevalence),
    early_reefs=paste(ReefID[is_early_reef],collapse=';'),.groups='drop')
nowcast_rows <- loeo |>
  left_join(early_summary,by='event_year') |>
  left_join(early_lookup |> select(event_year,ReefID,is_early_reef),
            by=c('event_year','ReefID')) |> filter(!is_early_reef)
nowcast_features <- unique(c(occurrence_specs[[selected_occurrence]],
                             'early_logit_delta'))
nowcast_predictions <- tibble()
for (fold_name in unique(nowcast_rows$fold)) {
  av <- nowcast_rows |> filter(fold==fold_name)
  tr <- nowcast_rows |> filter(fold!=fold_name)
  fit <- fit_occurrence_calibrator(tr,av,nowcast_features)
  nowcast_predictions <- bind_rows(nowcast_predictions,av |>
    mutate(predicted_occurrence_nowcast=fit$occurrence,
      predicted_mortality_nowcast=fit$occurrence*base_conditional_magnitude,
      residual_nowcast=observed_mortality-predicted_mortality_nowcast))
}
nowcast_comparison <- bind_rows(
  nowcast_rows |> mutate(candidate='Initial operational forecast'),
  nowcast_predictions |>
    mutate(predicted_occurrence=predicted_occurrence_nowcast,
      predicted_mortality=predicted_mortality_nowcast,
      candidate='Sequential early-prevalence update')) |>
  group_by(candidate) |> metric_summary_local()
nowcast_event_metrics <- bind_rows(
  nowcast_rows |> mutate(candidate='Initial operational forecast'),
  nowcast_predictions |>
    mutate(predicted_occurrence=predicted_occurrence_nowcast,
      predicted_mortality=predicted_mortality_nowcast,
      candidate='Sequential early-prevalence update')) |>
  group_by(candidate,event_year) |> metric_summary_local()
write_csv(early_summary,file.path(out_dir,'early_prevalence_summary.csv'))
write_csv(nowcast_predictions,file.path(out_dir,'nowcast_cv_predictions.csv'))
write_csv(nowcast_comparison,file.path(out_dir,'nowcast_model_comparison.csv'))
write_csv(nowcast_event_metrics,file.path(out_dir,'nowcast_event_metrics.csv'))

# Named 2020/2022 controls: model inputs, survey support and the operational
# local-first logger correction available for that reef-event.
audit_ids <- c('21-556','18-065','14-126','20-104','16-015','11-049',
  '11-162','18-051','21-139','21-588','20-351a','20-351b','21-060')
survey_support <- read_csv(file.path(root,'data','processed',
  'annual_coral_transitions.csv'),show_col_types=FALSE) |>
  transmute(source_observation_id=as.character(source_observation_id),
    baseline_survey_date,survey_date,interval_days,interval_years,
    pre_cover_survey=pre_cover,post_cover_survey=post_cover,cover_change_pp,
    depth_survey=depth,survey_sample_type,
    disturbance_text_survey=disturbance_text)
logger_support <- read_csv(file.path(root,'data','processed',
  'noaa_dhw_correction_layer_local_first_validation.csv'),
  show_col_types=FALSE) |>
  select(ReefID,event_year,local_first_correction,correction_source,
    correction_sd,nearest_local_logger_km,effective_local_loggers,
    direct_logger_site)
selected_initial <- occurrence_predictions |>
  filter(scheme=='leave_one_event_out',
         occurrence_candidate==selected_occurrence) |>
  mutate(predicted_mortality=predicted_mortality_initial,
    predicted_occurrence=predicted_occurrence_initial,
    residual=observed_mortality-predicted_mortality_initial)
named_audit <- selected_initial |>
  filter(ReefID %in% audit_ids,event_year %in% c(2020L,2022L)) |>
  left_join(survey_support,by='source_observation_id',relationship='many-to-one') |>
  left_join(logger_support,by=c('ReefID','event_year'),
            relationship='many-to-one') |>
  mutate(control_group=case_when(
    ReefID %in% c('21-556','18-065','14-126','20-104','16-015') ~
      'Leading La Nina residual / logger audit',
    ReefID %in% c('11-049','11-162','18-051') ~ 'High-DHW low-loss control',
    TRUE ~ '2022 low-DHW occurrence-floor control'),
    logger_evidence=case_when(
      str_detect(coalesce(correction_source,''),'direct') ~
        'Direct logger correction',
      str_detect(coalesce(correction_source,''),'local') ~
        'Nearby logger interpolation',
      TRUE ~ 'No supported local correction'),
    survey_support_flag=case_when(
      pre_cover_survey<.08 ~ 'Low starting cover',
      interval_years>1.5 ~ 'Long survey interval',
      TRUE ~ 'Standard annual support')) |>
  arrange(control_group,event_year,desc(abs(residual)))
write_csv(named_audit,file.path(out_dir,'named_reef_audit.csv'))

# Report/paper figures.
cots_plot_data <- cots_comparison |> filter(scheme=='leave_one_event_out') |>
  select(candidate,rmse,predictive_r2,severe_rmse) |>
  pivot_longer(-candidate,names_to='metric',values_to='value') |>
  mutate(candidate=recode(candidate,
    current_log_interval_selected='Current log interval',
    cots_raw_interval_relative='Raw interval density',
    cots_raw_current_timing_relative='Raw current + timing',
    cots_raw_current_timing_absolute='Raw current + timing (absolute loss)'),
    metric=recode(metric,rmse='RMSE',predictive_r2='Predictive R²',
                  severe_rmse='Severe RMSE'))
p_cots <- ggplot(cots_plot_data,aes(value,candidate,colour=candidate))+
  geom_point(size=3)+facet_wrap(~metric,scales='free_x')+
  labs(title='Uncompressed COTS pressure and compact timing screen',
    subtitle='Event-held-out validation; lower RMSE and higher R² are better',
    x=NULL,y=NULL)+theme_minimal(base_size=12)+theme(legend.position='none')
event_plot_data <- occurrence_event_metrics |>
  mutate(candidate=recode(occurrence_candidate,
    occurrence_offset_only='Current occurrence',
    occurrence_roni_rain_wqc='RONI + rainfall/WQC',
    occurrence_soi_rain_wqc='SOI + rainfall/WQC'))
p_event <- ggplot(event_plot_data,aes(factor(event_year),predicted_mean,
  colour=candidate,group=candidate))+geom_line()+geom_point(size=2.5)+
  geom_point(aes(y=observed_mean),colour='black',shape=4,size=3,stroke=1.2)+
  labs(title='Does event state lower the La Nina mortality floor?',
    subtitle='Black crosses are observed event means; colours are held-out predictions',
    x='Bleaching event',y='Mean relative mortality',colour=NULL)+
  theme_minimal(base_size=12)+theme(legend.position='bottom')
p_audit <- named_audit |> group_by(control_group,ReefID,ReefName,event_year) |>
  summarise(residual=mean(residual),.groups='drop') |>
  mutate(label=paste0(str_remove(ReefName,' \\([^()]+\\)$'),' ',event_year)) |>
  ggplot(aes(residual,reorder(label,residual),colour=control_group))+
  geom_vline(xintercept=0,linetype=2,colour='grey50')+
  geom_segment(aes(x=0,xend=residual,yend=reorder(label,residual)),linewidth=.7)+
  geom_point(size=2.5)+labs(title='Named La Nina reef controls',
    subtitle='Positive residuals are underpredictions; negative are overpredictions',
    x='Observed minus predicted relative mortality',y=NULL,colour=NULL)+
  theme_minimal(base_size=11)+theme(legend.position='bottom')
figures <- list(
  list(plot=p_cots,stem='Fig-INLA-16_cots_raw_timing_screen_v001',w=11,h=5.5),
  list(plot=p_event,stem='Fig-INLA-17_lanina_occurrence_state_v001',w=10,h=6),
  list(plot=p_audit,stem='Fig-INLA-18_lanina_named_reef_audit_v001',w=11,h=7))
for (item in figures) {
  ggsave(file.path(fig_dir,paste0(item$stem,'.png')),item$plot,
         width=item$w,height=item$h,dpi=300)
  ggsave(file.path(fig_dir,paste0(item$stem,'.pdf')),item$plot,
         width=item$w,height=item$h)
  ggsave(file.path(out_dir,paste0(item$stem,'.png')),item$plot,
         width=item$w,height=item$h,dpi=220)
}
save_figure_bundle(p_cots,'Fig-INLA-16_cots_raw_timing_screen',cots_plot_data,
  'Held-out comparison of log-compressed COTS intensity, raw interval density, current-year temporal state and an absolute-cover-loss biomass candidate.',
  'Tests whether retaining the original density contrast and compact temporal support improves transfer without redundant interactions.',
  'The selected interval maximum is retrospectively timed and raw density slightly worsens the Gannett severe positive control.',
  'cots_raw_interval_logistic20_5_cyclone','INLA composite',
  'model_comparison','operational_candidate_test',root,TRUE,
  code_source='scripts/test_cots_raw_enso_occurrence_controls.R',width=11,height=5.5)
save_figure_bundle(p_event,'Fig-INLA-17_lanina_occurrence_state',event_plot_data,
  'Observed and held-out predicted event means for the current occurrence layer and continuous RONI/SOI plus rainfall/WQC candidates.',
  'Shows that the tested climate-state corrections improve 2022 but worsen 2020 and suppress 2024.',
  'Only five events are available and WQC is a coloured-water/runoff proxy rather than measured discharge.',
  'cots_raw_interval_logistic20_5_cyclone','INLA occurrence calibration',
  'event_occurrence_state','operational_candidate_test',root,TRUE,
  code_source='scripts/test_cots_raw_enso_occurrence_controls.R',width=10,height=6)
save_figure_bundle(p_audit,'Fig-INLA-18_lanina_named_reef_audit',named_audit,
  'Held-out residuals for named 2020 and 2022 positive and negative reef controls.',
  'Connects the leading La Nina misses to local logger correction, composition, cooling and survey support.',
  'The audit is hypothesis-generating and repeated programme rows can describe the same reef-event.',
  'cots_raw_interval_logistic20_5_cyclone','INLA composite',
  'named_reef_residual_audit','operational_diagnostic',root,TRUE,
  code_source='scripts/test_cots_raw_enso_occurrence_controls.R',width=11,height=7)
write_figure_readme(root)
message('Selected COTS: ',selected_cots)
message('Selected initial occurrence: ',selected_occurrence)
