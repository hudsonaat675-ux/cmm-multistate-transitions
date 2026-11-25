#!/usr/bin/env Rscript
# f11: descriptive Aalen-Johansen state occupancy and CMD-free years.
#   - Death-inclusive occupancy for CHARLS and HRS (Figure 3 main panels).
#   - ELSA disease-state occupancy with death treated as censoring, clearly
#     labelled: ELSA leaves the death-inclusive cross-cohort comparison.
#   - CMD-free years over 8 years by cohort-specific corrected-MVS median group
#     (area under the no-CMD occupancy curve), for Fig S4 / Table S12 and the
#     Results sentence.

options(stringsAsFactors = FALSE, width = 170)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f11_occupancy_cmdfree")
write_csv <- make_write_csv(output_dir)

HORIZON_YEARS <- 8

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
scores <- inputs$scores

valid <- is.finite(transitions$tstart) & is.finite(transitions$tstop) & transitions$tstop > transitions$tstart
tv <- transitions[valid, , drop = FALSE]

# Multi-state counting-process outcome: state entered at tstop (censor = stay).
state_levels <- c("censor", "single_cmd", "cmm", "death")
build_ms <- function(x, censor_death = FALSE) {
  to_state <- ifelse(x$event == 1, x$to_state, "censor")
  if (censor_death) to_state[to_state == "death"] <- "censor"
  x$ms_state <- factor(ifelse(to_state == "no_cmd", "censor", to_state), levels = state_levels)
  x
}

# Occupancy from the baseline no-CMM population (starting states no_cmd or
# single_cmd handled via initial state in survfit through istate).
occupancy_rows <- list()
cmdfree_rows <- list()

for (cohort_name in COHORTS) {
  x <- tv[tv$cohort == cohort_name, , drop = FALSE]
  x <- x[order(x$person_id, x$tstart), , drop = FALSE]
  censor_death <- cohort_name == "ELSA"
  x <- build_ms(x, censor_death = censor_death)
  x$istate <- factor(x$from_state, levels = c("no_cmd", "single_cmd", "cmm"))

  fit <- survfit(Surv(tstart, tstop, ms_state) ~ 1, data = x, id = person_id, istate = istate)
  pst <- fit$pstate
  states <- fit$states
  keep_times <- fit$time <= HORIZON_YEARS
  occ <- data.frame(
    cohort = cohort_name, time_years = c(0, fit$time[keep_times]),
    stringsAsFactors = FALSE
  )
  p0 <- if (!is.null(fit$p0)) fit$p0 else pst[1, ]
  for (j in seq_along(states)) {
    nm <- states[j]
    nm_label <- if (nm %in% c("(s0)", "no_cmd")) "no_cmd" else nm
    occ[[paste0("occupancy_", nm_label)]] <- c(unname(p0[j]), pst[keep_times, j])
  }
  occ$death_handling <- if (censor_death) "death censored (incomplete ascertainment; not death-inclusive)" else "death included as absorbing state"
  occupancy_rows[[length(occupancy_rows) + 1]] <- occ

  # CMD-free years by corrected-MVS median group. Two population conventions:
  #   baseline_no_cmm_pool = whole baseline population including single-CMD
  #     starters (0 no-CMD years by construction) - the submitted S12 convention;
  #   baseline_no_cmd_starters = restricted to persons starting CMD-free.
  cutoff <- median(scores$corrected_per_sd[scores$cohort == cohort_name], na.rm = TRUE)
  for (population in c("baseline_no_cmm_pool", "baseline_no_cmd_starters")) {
    xg <- x[is.finite(x$corrected_per_sd), , drop = FALSE]
    if (population == "baseline_no_cmd_starters") xg <- xg[xg$baseline_state == "no_cmd", , drop = FALSE]
    xg$mvs_group <- ifelse(xg$corrected_per_sd >= cutoff, "high_mvs", "low_mvs")
    for (grp in c("low_mvs", "high_mvs")) {
      xx <- xg[xg$mvs_group == grp, , drop = FALSE]
      if (!nrow(xx)) next
      fit_g <- survfit(Surv(tstart, tstop, ms_state) ~ 1, data = xx, id = person_id, istate = xx$istate)
      s0_col <- which(fit_g$states %in% c("(s0)", "no_cmd"))
      if (length(s0_col) != 1) next
      keep_g <- fit_g$time <= HORIZON_YEARS
      share_start_no_cmd <- mean(tapply(xx$baseline_state, xx$person_id, function(z) z[1]) == "no_cmd")
      times <- c(0, fit_g$time[keep_g], HORIZON_YEARS)
      p0_path <- fit_g$pstate[keep_g, s0_col]
      p0_start <- if (population == "baseline_no_cmm_pool") share_start_no_cmd else 1
      p0 <- c(p0_start, p0_path, if (length(p0_path)) p0_path[length(p0_path)] else p0_start)
      auc <- sum(diff(times) * p0[-length(p0)])
      cmdfree_rows[[length(cmdfree_rows) + 1]] <- data.frame(
        cohort = cohort_name, population = population, mvs_group = grp,
        horizon_years = HORIZON_YEARS,
        persons = length(unique(xx$person_id)),
        cmd_free_years = auc,
        death_handling = if (censor_death) "death censored" else "death absorbing",
        stringsAsFactors = FALSE
      )
    }
  }
}

occupancy <- do.call(rbind, lapply(occupancy_rows, function(z) {
  cols <- c("cohort", "time_years", "occupancy_no_cmd", "occupancy_single_cmd",
            "occupancy_cmm", "occupancy_death", "death_handling")
  for (cc in cols) if (!cc %in% names(z)) z[[cc]] <- NA_real_
  z[, cols, drop = FALSE]
}))
write_csv(occupancy, "state_occupancy_aalen_johansen.csv")

cmdfree <- do.call(rbind, cmdfree_rows)
# Reshape and merge on cohort AND population: the two population conventions are
# separate rows for the same cohort, so a cohort-only key silently recycles the
# pool difference onto the no-CMD-starter rows.
wide <- reshape(cmdfree[, c("cohort", "population", "mvs_group", "cmd_free_years")],
                idvar = c("cohort", "population"), timevar = "mvs_group", direction = "wide")
wide$difference_low_minus_high <- wide$cmd_free_years.low_mvs - wide$cmd_free_years.high_mvs
cmdfree_out <- merge(cmdfree, wide[, c("cohort", "population", "difference_low_minus_high")],
                     by = c("cohort", "population"))
write_csv(cmdfree_out, "cmd_free_years_by_mvs_group.csv")

write_session_log("f11_occupancy_cmdfree", c(
  "ELSA occupancy: death censored and labelled; excluded from death-inclusive cross-cohort comparison.",
  "CMD-free years: area under no-CMD occupancy to 8 years, by cohort-specific corrected-MVS median group."
))
message("f11 complete: ", output_dir)
