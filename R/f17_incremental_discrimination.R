#!/usr/bin/env Rscript
# f17: incremental discrimination of the MVS.
# For each disease-accumulation transition, Harrell's C of the Model-3
# covariate set plus the social proxy (base) is compared with the same model
# plus the per-SD MVS (main effects; no product term, so the comparison is a
# single added variable). Variance of C and of the difference uses the
# infinitesimal jackknife clustered by person (survival::concordance).
# Discrimination is in-sample; no external risk model is available in all
# three cohorts, so this measures added information beyond these covariates only.

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f17_incremental_discrimination")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
ACCUMULATION <- TRANSITION_MAP[TRANSITION_MAP$transition %in% c("no_cmd_to_single_cmd", "no_cmd_to_cmm", "single_cmd_to_cmm"), ]

f_base <- Surv(tstart, tstop, event_target) ~ baseline_social_isolated_2item + baseline_age + baseline_sex +
  factor(baseline_education_level) + factor(baseline_wealth_tertile) + baseline_current_smoking +
  baseline_current_alcohol + baseline_cancer + baseline_lung_disease
f_age_sex <- Surv(tstart, tstop, event_target) ~ baseline_age + baseline_sex
add_mvs <- function(f) update(f, . ~ . + exposure)

rows <- list()
for (cohort_name in COHORTS) for (i in seq_len(nrow(ACCUMULATION))) {
  tr <- ACCUMULATION$transition[i]; fs <- ACCUMULATION$from_state[i]
  x <- prepare_model3_riskset(transitions, cohort_name, tr, fs)
  for (base_name in c("model3_covariates", "age_sex_only")) {
    fb <- if (base_name == "model3_covariates") f_base else f_age_sex
    z0 <- fit_capture(function() coxph(fb, data = x, ties = "efron"))
    z1 <- fit_capture(function() coxph(add_mvs(fb), data = x, ties = "efron"))
    if (is.null(z0$fit) || is.null(z1$fit)) next
    cc <- concordance(z0$fit, z1$fit, cluster = x$person_id)
    est <- coef(cc); v <- vcov(cc)
    d <- est[2] - est[1]; se_d <- sqrt(v[1, 1] + v[2, 2] - 2 * v[1, 2])
    rows[[length(rows) + 1]] <- data.frame(
      cohort = cohort_name, transition = tr, base_model = base_name,
      events = sum(x$event_target), unique_persons = length(unique(x$person_id)),
      c_base = est[1], c_base_se = sqrt(v[1, 1]), c_with_mvs = est[2], c_with_mvs_se = sqrt(v[2, 2]),
      delta_c = d, delta_c_ci_low = d - 1.96 * se_d, delta_c_ci_high = d + 1.96 * se_d,
      p_delta = 2 * pnorm(abs(d / se_d), lower.tail = FALSE),
      warnings = paste(c(z0$warning, z1$warning)[nzchar(c(z0$warning, z1$warning))], collapse = " | "),
      row.names = NULL)
  }
}
write_csv(do.call(rbind, rows), "mvs_incremental_concordance.csv")
write_session_log("f17_incremental_discrimination", c(
  "Harrell's C (counting-process form), person-clustered infinitesimal-jackknife variance.",
  "Base models: Model-3 covariates + proxy; and age + sex only."))
message("f17 complete: ", output_dir)
