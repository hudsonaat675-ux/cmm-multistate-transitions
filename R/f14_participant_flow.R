#!/usr/bin/env Rscript
# f14: full participant flow for Fig. 1a, from all
# respondents at the selected baseline wave to the primary complete-case
# sample, with mutually exclusive, sequential exclusion reasons. Uses the
# frozen conventions in _common.R for the pool, the MVS and the Model-3 sample.

options(stringsAsFactors = FALSE, width = 180)
source("_common.R")

output_dir <- make_output_dir("f14_participant_flow")
write_csv <- make_write_csv(output_dir)

PERSON_WAVE_DIR <- env_path("CMM_PERSON_WAVE_DIR")
BASE_WAVE <- c(CHARLS = 1, ELSA = 4, HRS = 10)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
baseline <- inputs$baseline
eligible <- inputs$eligible

model3_persons <- function(cohort_name) {
  unique(unlist(lapply(seq_len(nrow(TRANSITION_MAP)), function(i) {
    prepare_model3_riskset(transitions, cohort_name, TRANSITION_MAP$transition[i], TRANSITION_MAP$from_state[i])$person_id
  })))
}

flow_rows <- list(); reason_rows <- list(); covariate_rows <- list()
add_step <- function(cohort, order, step, n) {
  flow_rows[[length(flow_rows) + 1]] <<- data.frame(cohort = cohort, step_order = order, step = step, n = n)
}
add_reason <- function(cohort, after_step, reason, n) {
  reason_rows[[length(reason_rows) + 1]] <<- data.frame(cohort = cohort, excluded_before_step = after_step, reason = reason, n = n)
}

for (cohort_name in COHORTS) {
  p <- read.csv(file.path(PERSON_WAVE_DIR, paste0(tolower(cohort_name), "_person_wave.csv")), check.names = FALSE)
  p <- p[p$wave == BASE_WAVE[[cohort_name]], , drop = FALSE]
  b <- baseline[baseline$cohort == cohort_name & baseline$preferred_baseline_wave == 1, , drop = FALSE]
  b <- b[!duplicated(b$person_id), , drop = FALSE]

  # Step 0 -> 1: respondents at the baseline wave not carried into the candidate file.
  add_step(cohort_name, 0, sprintf("Respondents at baseline wave %d (harmonised person-wave file)", BASE_WAVE[[cohort_name]]), nrow(p))
  dropped <- p[!p$person_id %in% b$person_id, , drop = FALSE]
  if (nrow(dropped)) {
    why <- ifelse(is.na(dropped$cmd_count), "CMD status (diabetes/heart disease/stroke) missing at baseline",
           ifelse(is.na(dropped$age), "age missing", "other (not carried into baseline-candidate file)"))
    for (w in unique(why)) add_reason(cohort_name, 1, w, sum(why == w))
  }
  add_step(cohort_name, 1, "Baseline candidates with known CMD status", nrow(b))

  # Step 1 -> 2: age.
  add_reason(cohort_name, 2, "age < 45 years", sum(b$age < 45, na.rm = TRUE))
  add_reason(cohort_name, 2, "age missing", sum(is.na(b$age)))
  b2 <- b[b$eligible_age45 %in% 1, , drop = FALSE]
  add_step(cohort_name, 2, "Aged >= 45 years", nrow(b2))

  # Step 2 -> 3: prevalent CMM.
  add_reason(cohort_name, 3, "CMM (>= 2 of diabetes, heart disease, stroke) at baseline", sum(!b2$eligible_no_cmm %in% 1))
  b3 <- b2[b2$eligible_no_cmm %in% 1, , drop = FALSE]
  add_step(cohort_name, 3, "Free of CMM at baseline", nrow(b3))

  # Step 3 -> 4: follow-up.
  add_reason(cohort_name, 4, "no later interview or death record", sum(!b3$analysis_eligible_no_cmm %in% 1))
  pool <- eligible[eligible$cohort == cohort_name, , drop = FALSE]
  stopifnot(nrow(pool) == sum(b3$analysis_eligible_no_cmm %in% 1))
  add_step(cohort_name, 4, "Analysis-eligible no-CMM pool", nrow(pool))

  # Step 4 -> 5: exposure modules (MVS >= 3 of 4 components; two-item proxy).
  mvs_ok <- is.finite(pool$corrected_per_sd)
  proxy_ok <- !is.na(pool$social_isolated_2item)
  add_reason(cohort_name, 5, "MVS unavailable (< 3 of 4 components) and proxy unavailable", sum(!mvs_ok & !proxy_ok))
  add_reason(cohort_name, 5, "MVS unavailable (< 3 of 4 components), proxy available", sum(!mvs_ok & proxy_ok))
  add_reason(cohort_name, 5, "proxy unavailable, MVS available", sum(mvs_ok & !proxy_ok))
  exposed <- pool[mvs_ok & proxy_ok, , drop = FALSE]
  add_step(cohort_name, 5, "Exposure-available (MVS and social-connection proxy)", nrow(exposed))

  # Step 5 -> 6: Model-3 complete case (covariates + at least one valid interval).
  cc <- model3_persons(cohort_name)
  covar_base <- c(age = "age", sex = "sex", education = "education_level", wealth = "wealth_tertile",
                  smoking = "current_smoking", alcohol = "current_alcohol", cancer = "baseline_cancer",
                  lung_disease = "baseline_lung_disease")
  miss <- sapply(covar_base, function(v) is.na(exposed[[v]]))
  any_miss <- rowSums(miss) > 0
  in_cc <- exposed$person_id %in% cc
  add_reason(cohort_name, 6, "missing >= 1 Model-3 covariate", sum(any_miss & !in_cc))
  add_reason(cohort_name, 6, "no valid follow-up interval in transition data", sum(!any_miss & !in_cc))
  add_reason(cohort_name, 6, "check: complete covariates but excluded / missing covariate but included", sum(any_miss & in_cc))
  for (nm in names(covar_base)) covariate_rows[[length(covariate_rows) + 1]] <- data.frame(
    cohort = cohort_name, covariate = nm, n_missing_among_exposure_available = sum(miss[, nm]),
    n_missing_only_this = sum(miss[, nm] & rowSums(miss) == 1))
  add_step(cohort_name, 6, "Primary complete-case sample (Model 3)", length(cc))
  stopifnot(all(cc %in% exposed$person_id))
}

flow <- do.call(rbind, flow_rows)
reasons <- do.call(rbind, reason_rows)
reasons <- reasons[reasons$n > 0 | grepl("^check", reasons$reason), , drop = FALSE]
write_csv(flow, "participant_flow_steps.csv")
write_csv(reasons, "participant_flow_exclusion_reasons.csv")
write_csv(do.call(rbind, covariate_rows), "participant_flow_covariate_missingness.csv")

# Wide layout for the figure: one column per cohort.
wide <- reshape(flow[, c("cohort", "step_order", "step", "n")], idvar = c("step_order", "step"), timevar = "cohort", direction = "wide")
write_csv(wide[order(wide$step_order), ], "participant_flow_wide.csv")

write_session_log("f14_participant_flow", c(
  "Full participant flow from baseline-wave respondents to the Model-3 sample.",
  paste("Person-wave directory:", PERSON_WAVE_DIR)
))
message("f14 complete: ", output_dir)
