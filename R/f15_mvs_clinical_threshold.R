#!/usr/bin/env Rscript
# f15: MVS construction sensitivity.
# (1) Clinically defined count: general obesity, central obesity and raised
#     blood pressure (0-3), with WHO/NCEP thresholds in all cohorts and
#     Chinese criteria as a CHARLS sensitivity; also a binary >=2 criteria.
# (2) BMI-free continuous score: mean of within-cohort z-scores of waist, SBP
#     and DBP (all three required), re-standardised in the eligible sample.
# All models use the primary Model-3 interaction formula and risk-set rules.

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f15_mvs_clinical_threshold")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
eligible <- inputs$eligible

ACCUMULATION <- TRANSITION_MAP[TRANSITION_MAP$transition %in% c("no_cmd_to_single_cmd", "no_cmd_to_cmm", "single_cmd_to_cmm"), ]

flag <- function(value, cut) ifelse(is.na(value), NA_integer_, as.integer(value >= cut))
raised_bp <- function(sbp, dbp) {
  hi <- (sbp >= 140 & !is.na(sbp)) | (dbp >= 90 & !is.na(dbp))
  ifelse(hi, 1L, ifelse(!is.na(sbp) & !is.na(dbp), 0L, NA_integer_))
}

score_rows <- list(); prevalence_rows <- list()
for (cohort_name in COHORTS) {
  x <- eligible[eligible$cohort == cohort_name, , drop = FALSE]
  male <- x$sex == 1
  sets <- list(who = list(bmi = 30, waist_m = 102, waist_f = 88))
  if (cohort_name == "CHARLS") sets$china <- list(bmi = 28, waist_m = 90, waist_f = 85)
  out <- data.frame(cohort = cohort_name, person_id = x$person_id,
                    clin_count_china = NA_integer_, clin_ge2_china = NA_integer_)
  for (nm in names(sets)) {
    s <- sets[[nm]]
    ob <- flag(x$bmi, s$bmi)
    co <- ifelse(is.na(x$waist_cm) | is.na(male), NA_integer_, as.integer(x$waist_cm >= ifelse(male, s$waist_m, s$waist_f)))
    bp <- raised_bp(x$sbp, x$dbp)
    cnt <- ob + co + bp
    out[[paste0("clin_count_", nm)]] <- cnt
    out[[paste0("clin_ge2_", nm)]] <- as.integer(cnt >= 2)
    prevalence_rows[[length(prevalence_rows) + 1]] <- data.frame(
      cohort = cohort_name, criteria = nm, n_eligible = nrow(x), n_scored = sum(!is.na(cnt)),
      general_obesity_pct = 100 * mean(ob[!is.na(cnt)]), central_obesity_pct = 100 * mean(co[!is.na(cnt)]),
      raised_bp_pct = 100 * mean(bp[!is.na(cnt)]),
      count0_pct = 100 * mean(cnt == 0, na.rm = TRUE), count1_pct = 100 * mean(cnt == 1, na.rm = TRUE),
      count2_pct = 100 * mean(cnt == 2, na.rm = TRUE), count3_pct = 100 * mean(cnt == 3, na.rm = TRUE))
  }
  zw <- zscore(as.numeric(x$waist_cm)); zs <- zscore(as.numeric(x$sbp)); zd <- zscore(as.numeric(x$dbp))
  out$bmi_free_per_sd <- zscore(ifelse(is.finite(zw) & is.finite(zs) & is.finite(zd), (zw + zs + zd) / 3, NA_real_))
  score_rows[[length(score_rows) + 1]] <- out
}
scores <- do.call(rbind, score_rows)
write_csv(do.call(rbind, prevalence_rows), "clinical_criteria_prevalence.csv")

idx <- match(paste(transitions$cohort, transitions$person_id), paste(scores$cohort, scores$person_id))
variants <- c("corrected_per_sd", "bmi_free_per_sd", "clin_count_who", "clin_ge2_who", "clin_count_china", "clin_ge2_china")
for (v in setdiff(variants, "corrected_per_sd")) transitions[[v]] <- scores[[v]][idx]

formula_interaction <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

rows <- list(); diag_rows <- list()
for (cohort_name in COHORTS) for (i in seq_len(nrow(ACCUMULATION))) for (v in variants) {
  if (grepl("china", v) && cohort_name != "CHARLS") next
  tr <- ACCUMULATION$transition[i]; fs <- ACCUMULATION$from_state[i]
  x <- prepare_model3_riskset(transitions, cohort_name, tr, fs, exposure_var = v)
  z <- fit_capture(function() coxph(formula_interaction, data = x, ties = "efron"))
  diag_rows[[length(diag_rows) + 1]] <- data.frame(cohort = cohort_name, transition = tr, variant = v,
    risk_intervals = nrow(x), events = sum(x$event_target), unique_persons = length(unique(x$person_id)),
    status = if (is.null(z$fit)) "error" else "fit", warning_or_error = z$warning)
  if (is.null(z$fit)) next
  b <- coef(z$fit)[["exposure"]]; se <- sqrt(vcov(z$fit)["exposure", "exposure"])
  rows[[length(rows) + 1]] <- data.frame(cohort = cohort_name, transition = tr, variant = v,
    unit = c(corrected_per_sd = "per SD", bmi_free_per_sd = "per SD", clin_count_who = "per criterion",
             clin_ge2_who = ">=2 vs <2 criteria", clin_count_china = "per criterion",
             clin_ge2_china = ">=2 vs <2 criteria")[[v]],
    events = sum(x$event_target), unique_persons = length(unique(x$person_id)),
    beta = b, robust_se = se, hazard_ratio = exp(b), ci_low = exp(b - 1.96 * se), ci_high = exp(b + 1.96 * se),
    p_value = 2 * pnorm(abs(b / se), lower.tail = FALSE))
}
write_csv(do.call(rbind, rows), "mvs_threshold_and_bmi_free_cox_results.csv")
write_csv(do.call(rbind, diag_rows), "mvs_threshold_and_bmi_free_cox_diagnostics.csv")
write_session_log("f15_mvs_clinical_threshold", c(
  "WHO/NCEP criteria: BMI>=30; waist>=102 (men)/>=88 (women) cm; SBP>=140 or DBP>=90 mmHg.",
  "CHARLS China criteria: BMI>=28; waist>=90/>=85 cm; same BP rule."))
message("f15 complete: ", output_dir)
