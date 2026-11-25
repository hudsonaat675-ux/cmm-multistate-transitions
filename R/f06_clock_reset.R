#!/usr/bin/env Rscript
# f06: clock-forward (primary) versus clock-reset (sensitivity) time scales.
# Logic unchanged from the pre-analysis; inherits the half-month rule via the
# common loader. Baseline prevalent single-CMD persons have unknown pre-baseline
# sojourn time; reset time therefore starts at study baseline (documented).

options(stringsAsFactors = FALSE, width = 160)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f06_clock_reset")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions

valid <- is.finite(transitions$tstart) & is.finite(transitions$tstop) & transitions$tstop > transitions$tstart
group_key <- paste(transitions$cohort, transitions$person_id, transitions$from_state, sep = "__")
state_entry <- ave(ifelse(valid, transitions$tstart, NA_real_), group_key,
                   FUN = function(z) if (all(is.na(z))) NA_real_ else min(z, na.rm = TRUE))
transitions$reset_tstart <- transitions$tstart - state_entry
transitions$reset_tstop <- transitions$tstop - state_entry

fit_formula <- function(start_name, stop_name) {
  as.formula(paste0(
    "Surv(", start_name, ", ", stop_name, ", event_target) ~ exposure + baseline_social_isolated_2item + ",
    "exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) + ",
    "factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol + ",
    "baseline_cancer + baseline_lung_disease + cluster(person_id)"
  ))
}

result_rows <- list(); diagnostic_rows <- list()
for (cohort_name in COHORTS) {
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    tr <- TRANSITION_MAP$transition[i]
    from <- TRANSITION_MAP$from_state[i]
    x <- transitions[transitions$cohort == cohort_name & transitions$from_state == from & valid, , drop = FALSE]
    x$event_target <- as.integer(x$event == 1 & x$transition == tr)
    x$exposure <- as.numeric(x$corrected_per_sd)
    x$exposure_x_proxy <- x$exposure * as.numeric(x$baseline_social_isolated_2item)
    required <- c("person_id", "tstart", "tstop", "reset_tstart", "reset_tstop", "event_target",
                  "exposure", "exposure_x_proxy", MODEL3_COVARS)
    x <- x[complete.cases(x[, required, drop = FALSE]), , drop = FALSE]
    for (time_scale in c("clock_forward", "clock_reset")) {
      start_name <- if (time_scale == "clock_forward") "tstart" else "reset_tstart"
      stop_name <- if (time_scale == "clock_forward") "tstop" else "reset_tstop"
      z <- fit_capture(function() coxph(fit_formula(start_name, stop_name), data = x, ties = "efron", x = TRUE))
      diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = tr, time_scale = time_scale,
        risk_intervals = nrow(x), events = sum(x$event_target), unique_persons = length(unique(x$person_id)),
        status = if (is.null(z$fit)) "error" else "fit", warning_or_error = z$warning, stringsAsFactors = FALSE
      )
      if (is.null(z$fit)) next
      term_map <- c(exposure = "mvs", baseline_social_isolated_2item = "proxy", exposure_x_proxy = "mvs_x_proxy")
      b <- coef(z$fit)[names(term_map)]
      se <- sqrt(diag(vcov(z$fit)))[names(term_map)]
      result_rows[[length(result_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = tr, time_scale = time_scale,
        term = unname(term_map), beta = unname(b), robust_se = unname(se),
        hazard_ratio = exp(unname(b)), ci_low = exp(unname(b) - 1.96 * unname(se)),
        ci_high = exp(unname(b) + 1.96 * unname(se)),
        p_value = 2 * pnorm(abs(unname(b / se)), lower.tail = FALSE), stringsAsFactors = FALSE
      )
    }
  }
}

results <- do.call(rbind, result_rows)
write_csv(results, "clock_forward_vs_reset_results.csv")
write_csv(do.call(rbind, diagnostic_rows), "clock_forward_vs_reset_diagnostics.csv")

wide <- reshape(
  results[results$term == "mvs", c("cohort", "transition", "time_scale", "hazard_ratio", "ci_low", "ci_high", "p_value")],
  idvar = c("cohort", "transition"), timevar = "time_scale", direction = "wide"
)
write_csv(wide, "mvs_clock_forward_vs_reset_comparison.csv")

write_session_log("f06_clock_reset", c(
  "Clock-forward is the actual primary time scale; clock-reset is sensitivity (eMethods statement to be corrected accordingly)."
))
message("f06 complete: ", output_dir)
