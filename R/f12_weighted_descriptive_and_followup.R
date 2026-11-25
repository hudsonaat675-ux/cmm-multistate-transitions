# f12: (a) descriptive weighted Cox point estimates.
#        HRS: physical-measure weights (column baseline_pm_weight of the analysis
#             dataset; source HRS Cross-Wave Tracker, MPMWGTR).
#        CHARLS: Wave-1 cross-sectional individual weights with non-response
#             adjustment (ind_weight_ad2, official 2011 file weight.dta).
#        Reported as descriptive point estimates only, because the survey design
#        (strata and sampling units) is not used.
#      (b) follow-up summary per cohort.
#
# Run from the R/ folder:  Rscript f12_weighted_descriptive_and_followup.R

options(stringsAsFactors = FALSE, width = 170)
library(survival)
library(haven)      # read_dta(), to read the official Stata file
source("_common.R")

output_dir <- make_output_dir("f12_weighted_followup")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions

# ---- CHARLS Wave-1 weights from the official 2011 file weight.dta ----------
cw <- read_dta(env_path("CMM_CHARLS_2011_WEIGHT_DTA"))
cw <- as.data.frame(cw)
names(cw) <- tolower(names(cw))

# The 2011 identifier has 11 characters (9-character household + 2-character
# person). Later releases use 12 characters, with "0" inserted after the
# household part. The analysis dataset uses the 12-character form, stored as a
# whole number without leading zeros after read.csv().
id_text <- function(x) {
  sprintf("%.0f", suppressWarnings(as.numeric(as.character(x))))
}
id_2011 <- as.character(cw$id)
id_release <- paste0(substr(id_2011, 1, 9), "0", substr(id_2011, 10, 11))
cw$key <- paste("CHARLS", id_text(id_release), sep = "__")
cw$w_ad2 <- suppressWarnings(as.numeric(cw$ind_weight_ad2))
tr_key <- paste(transitions$cohort, id_text(transitions$person_id), sep = "__")
transitions$charls_w <- cw$w_ad2[match(tr_key, cw$key)]

wt_summary <- data.frame(
  cohort = c("CHARLS", "HRS"),
  weight = c("ind_weight_ad2 (CHARLS 2011 weight.dta; cross-sectional, non-response adjusted)",
             "baseline_pm_weight (analysis dataset; HRS Tracker MPMWGTR)"),
  positive_rows = c(sum(is.finite(transitions$charls_w) & transitions$charls_w > 0 & transitions$cohort == "CHARLS"),
                    sum(is.finite(transitions$baseline_pm_weight) & transitions$baseline_pm_weight > 0 & transitions$cohort == "HRS")),
  stringsAsFactors = FALSE
)
write_csv(wt_summary, "weight_source_summary.csv")

cox_formula <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

runs <- list(
  list(cohort = "CHARLS", weight_col = "charls_w", label = "charls_wave1_weighted"),
  list(cohort = "HRS", weight_col = "baseline_pm_weight", label = "hrs_pm_weighted")
)
res_rows <- list(); diag_rows <- list()
for (spec in runs) {
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    tr <- TRANSITION_MAP$transition[i]; fs <- TRANSITION_MAP$from_state[i]
    x <- prepare_model3_riskset(transitions, spec$cohort, tr, fs)
    x$analysis_weight <- as.numeric(x[[spec$weight_col]])
    x <- x[is.finite(x$analysis_weight) & x$analysis_weight > 0, , drop = FALSE]
    z <- fit_capture(function() coxph(cox_formula, data = x, weights = analysis_weight, ties = "efron", x = TRUE))
    diag_rows[[length(diag_rows) + 1]] <- data.frame(
      analysis = spec$label, cohort = spec$cohort, transition = tr,
      risk_intervals = nrow(x), events = sum(x$event_target),
      unique_persons = length(unique(x$person_id)),
      status = if (is.null(z$fit)) "error" else "fit", warning_or_error = z$warning,
      stringsAsFactors = FALSE
    )
    if (is.null(z$fit)) next
    term_map <- c(exposure = "mvs", baseline_social_isolated_2item = "proxy", exposure_x_proxy = "mvs_x_proxy")
    b <- coef(z$fit)[names(term_map)]
    res_rows[[length(res_rows) + 1]] <- data.frame(
      analysis = spec$label, cohort = spec$cohort, transition = tr,
      term = unname(term_map), weighted_hr_point = exp(unname(b)),
      inference = "descriptive point estimate only (survey design not independently reproducible)",
      risk_intervals = nrow(x), events = sum(x$event_target), stringsAsFactors = FALSE
    )
  }
}
write_csv(do.call(rbind, res_rows), "weighted_descriptive_point_estimates.csv")
write_csv(do.call(rbind, diag_rows), "weighted_descriptive_diagnostics.csv")

# ---- follow-up summary (Model 3 persons) ------------------------------------
fu_rows <- list()
for (cohort_name in COHORTS) {
  persons <- unique(unlist(lapply(seq_len(nrow(TRANSITION_MAP)), function(i) {
    x <- prepare_model3_riskset(transitions, cohort_name, TRANSITION_MAP$transition[i], TRANSITION_MAP$from_state[i])
    unique(x$person_id)
  })))
  valid <- transitions$cohort == cohort_name & is.finite(transitions$tstart) &
    is.finite(transitions$tstop) & transitions$tstop > transitions$tstart &
    transitions$person_id %in% persons
  fu <- tapply(transitions$tstop[valid], transitions$person_id[valid], max)
  q <- quantile(fu, c(.25, .5, .75), names = FALSE)
  fu_rows[[length(fu_rows) + 1]] <- data.frame(
    cohort = cohort_name, persons = length(fu),
    median_followup_years = q[2], iqr_low = q[1], iqr_high = q[3], stringsAsFactors = FALSE
  )
}
write_csv(do.call(rbind, fu_rows), "followup_summary_model3_persons.csv")

write_session_log("f12_weighted_descriptive_and_followup", c(
  "CHARLS weights: ind_weight_ad2 from the official 2011 weight.dta; HRS physical-measure weights from the analysis dataset.",
  "Weighted results are descriptive point estimates only."
))
message("f12 complete: ", output_dir)
