#!/usr/bin/env Rscript
# f00: definitive risk-set/event-count ledger under the frozen conventions,

options(stringsAsFactors = FALSE, width = 160)
source("_common.R")

output_dir <- make_output_dir("f00_riskset_counts")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
eligible <- inputs$eligible

manifest <- data.frame(
  item = c("analysis_id", "global_seed", "transitions_md5", "baseline_md5",
           "half_month_rule_fixed_intervals", "nonpositive_nondeath_dropped",
           "r_version"),
  value = c(ANALYSIS_ID, GLOBAL_SEED, inputs$input_md5["transitions"], inputs$input_md5["baseline"],
            inputs$n_half_month_fixed, inputs$n_nonpositive_dropped, R.version.string),
  stringsAsFactors = FALSE
)
write_csv(manifest, "analysis_id_manifest.csv")

# Disposition of the previously zero-length death intervals, by cohort/transition.
fixed <- transitions[transitions$half_month_rule_applied == 1, , drop = FALSE]
if (nrow(fixed)) {
  disposition <- aggregate(list(n_intervals = fixed$event),
                           by = list(cohort = fixed$cohort, transition = fixed$transition), FUN = length)
} else {
  disposition <- data.frame(cohort = character(), transition = character(), n_intervals = integer())
}
write_csv(disposition, "half_month_rule_disposition.csv")

# Count ladder: raw -> valid-time -> Model 3 complete case, per cohort x transition.
count_rows <- list()
for (cohort_name in COHORTS) {
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    tr <- TRANSITION_MAP$transition[i]
    fs <- TRANSITION_MAP$from_state[i]
    x0 <- transitions[transitions$cohort == cohort_name & transitions$from_state == fs, , drop = FALSE]
    x0$event_target <- as.integer(x0$event == 1 & x0$transition == tr)
    valid <- is.finite(x0$tstart) & is.finite(x0$tstop) & x0$tstop > x0$tstart
    x1 <- x0[valid, , drop = FALSE]
    x2 <- prepare_model3_riskset(transitions, cohort_name, tr, fs)
    count_rows[[length(count_rows) + 1]] <- data.frame(
      cohort = cohort_name, transition = tr, from_state = fs,
      raw_intervals = nrow(x0), raw_events = sum(x0$event_target, na.rm = TRUE),
      valid_time_intervals = nrow(x1), valid_time_events = sum(x1$event_target, na.rm = TRUE),
      model3_intervals = nrow(x2), model3_events = sum(x2$event_target),
      model3_unique_persons = length(unique(x2$person_id)),
      stringsAsFactors = FALSE
    )
  }
}
count_ladder <- do.call(rbind, count_rows)
write_csv(count_ladder, "riskset_event_count_ladder.csv")

# Number of participants under each counting rule.
person_rows <- list()
for (cohort_name in COHORTS) {
  model3_persons <- unique(unlist(lapply(seq_len(nrow(TRANSITION_MAP)), function(i) {
    x <- prepare_model3_riskset(transitions, cohort_name, TRANSITION_MAP$transition[i], TRANSITION_MAP$from_state[i])
    unique(x$person_id)
  })))
  valid_rows <- transitions$cohort == cohort_name & is.finite(transitions$tstart) &
    is.finite(transitions$tstop) & transitions$tstop > transitions$tstart
  person_rows[[length(person_rows) + 1]] <- data.frame(
    cohort = cohort_name,
    eligible_persons = sum(eligible$cohort == cohort_name),
    valid_transition_persons = length(unique(transitions$person_id[valid_rows])),
    model3_persons = length(model3_persons),
    stringsAsFactors = FALSE
  )
}
write_csv(do.call(rbind, person_rows), "person_count_unification.csv")

# MVS scale audit: the corrected per-SD score must have SD = 1 in eligible sample.
scale_rows <- do.call(rbind, lapply(COHORTS, function(cohort_name) {
  z <- inputs$scores$corrected_per_sd[inputs$scores$cohort == cohort_name]
  raw <- eligible$mvs_anthro_bp_z[eligible$cohort == cohort_name]
  data.frame(cohort = cohort_name,
             raw_composite_sd = sd(as.numeric(raw), na.rm = TRUE),
             corrected_per_sd_mean = mean(z, na.rm = TRUE),
             corrected_per_sd_sd = sd(z, na.rm = TRUE),
             n_scored = sum(is.finite(z)), stringsAsFactors = FALSE)
}))
write_csv(scale_rows, "mvs_scale_audit.csv")

write_session_log("f00_riskset_and_counts", c(
  paste("Half-month-rule fixed intervals:", inputs$n_half_month_fixed),
  paste("Non-death nonpositive intervals dropped:", inputs$n_nonpositive_dropped)
))
message("f00 complete: ", output_dir)
