#!/usr/bin/env Rscript
# f08: feasibility of a panel-observed continuous-time multistate model (msm),
# States: 1 no CMD, 2 single CMD,
# 3 CMM, 4 death (exact death times). Direct 1->3 is not modelled as an
# instantaneous intensity; interval censoring accounts for observed 1->3 jumps
# between visits. Covariates kept minimal (corrected per-SD MVS, proxy, baseline
# age/10, sex) for identifiability. Whatever the outcome, this script writes a
# feasibility record usable in the response letter.

options(stringsAsFactors = FALSE, width = 170)
if (!requireNamespace("msm", quietly = TRUE)) stop("Package msm required")
suppressPackageStartupMessages(library(msm))
source("_common.R")

output_dir <- make_output_dir("f08_msm_feasibility")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions

state_code_map <- c(no_cmd = 1, single_cmd = 2, cmm = 3, death = 4)

# Reconstruct panel observations per person: baseline state at first valid
# tstart, then the destination state at each interval's tstop.
valid <- is.finite(transitions$tstart) & is.finite(transitions$tstop) & transitions$tstop > transitions$tstart
tv <- transitions[valid, , drop = FALSE]
tv <- tv[order(tv$cohort, tv$person_id, tv$tstart), , drop = FALSE]

panel_parts <- list()
for (cohort_name in COHORTS) {
  x <- tv[tv$cohort == cohort_name, , drop = FALSE]
  first_row <- !duplicated(x$person_id)
  base_obs <- data.frame(
    cohort = cohort_name, person_id = x$person_id[first_row],
    time = x$tstart[first_row], state = unname(state_code_map[x$from_state[first_row]]),
    stringsAsFactors = FALSE
  )
  follow_obs <- data.frame(
    cohort = cohort_name, person_id = x$person_id,
    time = x$tstop, state = unname(state_code_map[x$to_state]),
    stringsAsFactors = FALSE
  )
  covars <- data.frame(
    cohort = cohort_name, person_id = x$person_id[first_row],
    mvs = as.numeric(x$corrected_per_sd[first_row]),
    proxy = as.numeric(x$baseline_social_isolated_2item[first_row]),
    age10 = as.numeric(x$baseline_age[first_row]) / 10,
    sex = as.numeric(x$baseline_sex[first_row]),
    stringsAsFactors = FALSE
  )
  obs <- rbind(base_obs, follow_obs)
  obs <- merge(obs, covars, by = c("cohort", "person_id"))
  obs <- obs[complete.cases(obs[, c("time", "state", "mvs", "proxy", "age10", "sex")]), , drop = FALSE]
  obs <- obs[order(obs$person_id, obs$time), , drop = FALSE]
  # msm requires strictly increasing times within subject; drop duplicated times
  # (state at a repeated time is the later record, which reflects the update).
  keep <- !duplicated(obs[, c("person_id", "time")], fromLast = TRUE)
  obs <- obs[keep, , drop = FALSE]
  # Keep only subjects with at least two observations.
  n_obs <- table(obs$person_id)
  obs <- obs[obs$person_id %in% names(n_obs)[n_obs >= 2], , drop = FALSE]
  panel_parts[[cohort_name]] <- obs
}

Q <- rbind(
  c(0, 0.1, 0,   0.05),
  c(0, 0,   0.1, 0.05),
  c(0, 0,   0,   0.1),
  c(0, 0,   0,   0)
)
dimnames(Q) <- list(c("no_cmd", "single_cmd", "cmm", "death"), c("no_cmd", "single_cmd", "cmm", "death"))

feasibility_rows <- list()
hr_rows <- list()

for (cohort_name in COHORTS) {
  obs <- panel_parts[[cohort_name]]
  elapsed <- NA_real_
  status <- "not_attempted"
  note <- ""
  minus2ll <- NA_real_
  fit <- NULL
  t0 <- proc.time()["elapsed"]
  # The logged "Absorbing - absorbing transition at observation k" warning is a
  # property of the input panel, not of the fit, so it is diagnosed here from the
  # data rather than by intercepting warnings. A calling handler around msm is
  # not viable: the optimiser can raise a warning per iteration, and muffling
  # each one by restart costs far more than R's own "keep the first 50" default,
  # to the point where the ELSA fit no longer finishes.
  ord <- obs[order(obs$person_id, obs$time), c("person_id", "state", "time")]
  same_person <- c(FALSE, ord$person_id[-1] == ord$person_id[-nrow(ord)])
  prev_absorbing <- c(FALSE, ord$state[-nrow(ord)] == 4)
  absorbing_followed <- same_person & prev_absorbing
  n_absorbing_followed <- sum(absorbing_followed)

  attempt <- tryCatch({
    qinit <- crudeinits.msm(state ~ time, subject = person_id, data = obs, qmatrix = Q)
    msm(state ~ time, subject = person_id, data = obs,
        qmatrix = qinit, deathexact = 4,
        covariates = ~ mvs + proxy + age10 + sex,
        control = list(maxit = 500, fnscale = 40000, reltol = 1e-9))
  }, error = function(e) e)
  elapsed <- unname(proc.time()["elapsed"] - t0)
  if (inherits(attempt, "error")) {
    status <- "error"
    note <- conditionMessage(attempt)
  } else {
    fit <- attempt
    conv <- fit$opt$convergence
    minus2ll <- fit$minus2loglik
    status <- if (!is.null(conv) && conv == 0) "converged" else paste0("optimizer_code_", conv)
    fnd <- tryCatch(all(is.finite(sqrt(diag(fit$covmat)))), error = function(e) FALSE)
    if (!fnd) {
      status <- paste(status, "nonfinite_se", sep = "+")
      note <- "covariance matrix not positive definite or SEs non-finite"
    }
  }
  feasibility_rows[[length(feasibility_rows) + 1]] <- data.frame(
    cohort = cohort_name, subjects = length(unique(obs$person_id)), observations = nrow(obs),
    status = status, minus2loglik = minus2ll, elapsed_seconds = round(elapsed, 1),
    elsa_death_note = if (cohort_name == "ELSA") "death ascertainment incomplete; mortality intensities not interpretable" else "",
    observations_after_absorbing_state = n_absorbing_followed,
    absorbing_state_note = if (n_absorbing_followed == 0)
      "no observation follows the absorbing death state; the msm absorbing-to-absorbing message comes from its internal deathexact handling and does not indicate duplicated death records"
      else "panel contains observations after the absorbing state; check the death-record construction",
    note = note, stringsAsFactors = FALSE
  )
  if (!is.null(fit) && grepl("^converged", status)) {
    hr <- tryCatch(hazard.msm(fit), error = function(e) NULL)
    if (!is.null(hr)) {
      for (cov_name in names(hr)) {
        tab <- hr[[cov_name]]
        hr_rows[[length(hr_rows) + 1]] <- data.frame(
          cohort = cohort_name, covariate = cov_name, transition = rownames(tab),
          hazard_ratio = tab[, "HR"], ci_low = tab[, "L"], ci_high = tab[, "U"],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

write_csv(do.call(rbind, feasibility_rows), "msm_feasibility_summary.csv")
if (length(hr_rows)) write_csv(do.call(rbind, hr_rows), "msm_hazard_ratios.csv")

write_session_log("f08_msm_feasibility", c(
  "States: 1 no CMD, 2 single CMD, 3 CMM, 4 death (deathexact).",
  "Direct 1->3 handled by interval censoring, not an instantaneous intensity.",
  "Covariates minimal for identifiability: corrected per-SD MVS, proxy, baseline age/10, sex."
))
message("f08 complete: ", output_dir)
