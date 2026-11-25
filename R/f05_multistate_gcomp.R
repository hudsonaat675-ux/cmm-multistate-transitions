#!/usr/bin/env Rscript
# f05: common-target-population 5-year multistate g-computation (product
# integral over all six cause-specific hazards) with RD/RR/RERI.
# Formal-run changes: person-level bootstrap B = 500 parallelised with
# per-replicate deterministic seeds (reproducible regardless of core count);
# ordinal SES retained in the transition models (documented).

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(parallel))
source("_common.R")

output_dir <- make_output_dir("f05_multistate_gcomp")
write_csv <- make_write_csv(output_dir)

BOOTSTRAP_B <- as.integer(Sys.getenv("FORMAL_BOOTSTRAP_B", "500"))
CORES <- as.integer(Sys.getenv("FORMAL_BOOTSTRAP_CORES", "6"))
# mclapply() cannot use several cores on Windows. Each bootstrap sample sets its
# own random seed, so the results are the same with one core or several.
if (.Platform$OS.type == "windows") {
  CORES <- 1L
}
horizon <- 5

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
eligible <- inputs$eligible
scores <- inputs$scores

scenario_map <- data.frame(
  scenario = c("R00", "R10", "R01", "R11"),
  high_mvs = c(0, 1, 0, 1),
  proxy = c(0, 0, 1, 1),
  stringsAsFactors = FALSE
)

model_covars <- c(
  "baseline_age", "baseline_sex", "baseline_education_level", "baseline_wealth_tertile",
  "baseline_current_smoking", "baseline_current_alcohol", "baseline_cancer",
  "baseline_lung_disease", "baseline_social_isolated_2item", "corrected_per_sd"
)

cutoff_rows <- list()
for (cohort_name in COHORTS) {
  cutoff <- median(scores$corrected_per_sd[scores$cohort == cohort_name], na.rm = TRUE)
  cutoff_rows[[length(cutoff_rows) + 1]] <- data.frame(
    cohort = cohort_name, corrected_mvs_median_cutoff = cutoff, stringsAsFactors = FALSE
  )
  transitions$high_mvs[transitions$cohort == cohort_name] <- as.integer(
    transitions$corrected_per_sd[transitions$cohort == cohort_name] >= cutoff
  )
}
write_csv(do.call(rbind, cutoff_rows), "binary_mvs_cutoffs.csv")

valid_complete <- transitions[
  is.finite(transitions$tstart) & is.finite(transitions$tstop) & transitions$tstop > transitions$tstart &
    complete.cases(transitions[, c("person_id", model_covars), drop = FALSE]),
  , drop = FALSE
]
selected_keys <- unique(paste(valid_complete$cohort, valid_complete$person_id, sep = "__"))

target <- data.frame(
  cohort = eligible$cohort, person_id = eligible$person_id, baseline_state = eligible$state,
  baseline_age = eligible$age, baseline_sex = eligible$sex,
  baseline_education_level = eligible$education_level, baseline_wealth_tertile = eligible$wealth_tertile,
  baseline_current_smoking = eligible$current_smoking, baseline_current_alcohol = eligible$current_alcohol,
  baseline_cancer = eligible$baseline_cancer, baseline_lung_disease = eligible$baseline_lung_disease,
  baseline_social_isolated_2item = eligible$social_isolated_2item,
  corrected_per_sd = eligible$corrected_per_sd,
  stringsAsFactors = FALSE
)
target$key <- paste(target$cohort, target$person_id, sep = "__")
target <- target[target$key %in% selected_keys & complete.cases(target[, c("baseline_state", model_covars), drop = FALSE]), , drop = FALSE]
for (i in seq_along(COHORTS)) {
  cutoff <- cutoff_rows[[i]]$corrected_mvs_median_cutoff
  sel <- target$cohort == COHORTS[i]
  target$high_mvs[sel] <- as.integer(target$corrected_per_sd[sel] >= cutoff)
}
target$exposure_x_proxy <- target$high_mvs * target$baseline_social_isolated_2item

target_counts <- as.data.frame(table(target$cohort, target$baseline_state), stringsAsFactors = FALSE)
names(target_counts) <- c("cohort", "baseline_state", "target_persons")
target_counts$mortality_boundary <- ifelse(
  target_counts$cohort == "ELSA",
  "death ascertainment incomplete; excluded from main absolute-probability displays",
  "included in main absolute-probability displays"
)
write_csv(target_counts, "common_target_population_counts.csv")

gcomp_formula <- Surv(tstart, tstop, event_target) ~ high_mvs + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + baseline_education_level +
  baseline_wealth_tertile + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

fit_transition_models <- function(data, cohort_name) {
  fits <- vector("list", nrow(TRANSITION_MAP))
  diagnostics <- list()
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    tr <- TRANSITION_MAP$transition[i]
    from <- TRANSITION_MAP$from_state[i]
    x <- data[
      data$cohort == cohort_name & data$from_state == from &
        is.finite(data$tstart) & is.finite(data$tstop) & data$tstop > data$tstart,
      , drop = FALSE
    ]
    x$event_target <- as.integer(x$event == 1 & x$transition == tr)
    x$exposure_x_proxy <- x$high_mvs * x$baseline_social_isolated_2item
    required <- c("person_id", "tstart", "tstop", "event_target", "high_mvs", "exposure_x_proxy", model_covars[-10])
    x <- x[complete.cases(x[, unique(required), drop = FALSE]), , drop = FALSE]
    warning_text <- character(); error_text <- NULL
    fit <- NULL
    if (sum(x$event_target) > 0) {
      fit <- tryCatch(
        withCallingHandlers(
          coxph(gcomp_formula, data = x, ties = "efron", x = TRUE),
          warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") }
        ),
        error = function(e) { error_text <<- conditionMessage(e); NULL }
      )
    }
    if (sum(x$event_target) > 0 && is.null(fit)) {
      return(list(fits = NULL, diagnostics = data.frame(
        cohort = cohort_name, transition = tr, risk_intervals = nrow(x), events = sum(x$event_target),
        status = "error", warning_or_error = paste(unique(c(warning_text, error_text)), collapse = " | "),
        stringsAsFactors = FALSE
      )))
    }
    fits[i] <- list(fit)
    diagnostics[[length(diagnostics) + 1]] <- data.frame(
      cohort = cohort_name, transition = tr, risk_intervals = nrow(x), events = sum(x$event_target),
      status = if (sum(x$event_target) == 0) "zero_event_hazard_set_to_zero" else "fit",
      warning_or_error = paste(unique(c(warning_text, error_text)), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
  names(fits) <- TRANSITION_MAP$transition
  list(fits = fits, diagnostics = do.call(rbind, diagnostics))
}

hazard_increment_matrix <- function(fits, horizon) {
  increment_lists <- lapply(fits, function(fit) {
    if (is.null(fit)) return(data.frame(time = numeric(), increment = numeric()))
    bh <- basehaz(fit, centered = FALSE)
    bh <- bh[is.finite(bh$time) & bh$time <= horizon, , drop = FALSE]
    if (!nrow(bh)) return(data.frame(time = numeric(), increment = numeric()))
    data.frame(time = bh$time, increment = c(bh$hazard[1], diff(bh$hazard)))
  })
  times <- sort(unique(unlist(lapply(increment_lists, function(z) z$time))))
  dh <- matrix(0, nrow = length(times), ncol = length(fits), dimnames = list(NULL, names(fits)))
  for (j in seq_along(increment_lists)) {
    z <- increment_lists[[j]]
    if (nrow(z)) dh[match(z$time, times), j] <- z$increment
  }
  list(times = times, increments = dh)
}

predict_rate_multipliers <- function(fits, newdata) {
  out <- matrix(1, nrow = nrow(newdata), ncol = length(fits), dimnames = list(NULL, names(fits)))
  for (j in seq_along(fits)) {
    if (!is.null(fits[[j]])) {
      out[, j] <- exp(predict(fits[[j]], newdata = newdata, type = "lp", reference = "zero"))
    }
  }
  out
}

# Product integral over the six cause-specific hazards.
#
# Each observed event time carries a piecewise-constant generator. A single
# one-step update per event time cannot represent the no CMD -> single CMD ->
# CMM path inside one increment, and the Breslow increments are not uniformly
# small (the largest CHARLS row increment is 3.42), so that approximation is not
# adequate here. The step is instead solved exactly by uniformisation
# (Jensen's method): with Lambda the largest row rate in the step and
# R = I + Q / Lambda a stochastic matrix, exp(Q) = sum_k Pois(k; Lambda) R^k,
# a nonnegative series with no cancellation. Terms are taken until the remaining
# Poisson mass is below UNIFORMISATION_TOLERANCE, so the transition matrix is
# exact to that tolerance rather than to first order in the step size.
#
# Cumulative incidences are read off an expanded, still-Markov state space that
# records how each absorbing state was reached:
#   1 no CMD | 2 single CMD | 3a CMM entered directly from no CMD
#   3b CMM entered from single CMD | 4a death from no CMD | 4b death from single
#   CMD | 4ca death from CMM entered directly | 4cb death from CMM via single CMD
# so every cause-specific cumulative incidence is a sum of expanded-state
# occupancies and needs no separate numerical integration.
UNIFORMISATION_TOLERANCE <- as.numeric(Sys.getenv("FORMAL_GCOMP_UNIF_TOL", "1e-13"))

poisson_terms_needed <- function(lambda, tol) {
  if (!is.finite(lambda) || lambda <= 0) return(0L)
  k <- max(10L, as.integer(ceiling(lambda + 8 * sqrt(lambda) + 12)))
  while (stats::ppois(k, lambda, lower.tail = FALSE) > tol) k <- k + 5L
  k
}

propagate_multistate <- function(dh, multipliers, initial_state) {
  n <- nrow(multipliers)
  S <- c("no_cmd", "single_cmd", "cmm_direct", "cmm_via_single",
         "death_from_no", "death_from_single", "death_cmm_direct", "death_cmm_via_single")
  v <- matrix(0, nrow = n, ncol = 8, dimnames = list(NULL, S))
  v[, if (initial_state == "no_cmd") 1L else 2L] <- 1
  max_row_hazard_increment <- 0
  max_poisson_terms <- 0L
  for (k in seq_len(nrow(dh))) {
    h <- sweep(multipliers, 2, dh[k, ], `*`)
    a12 <- h[, 1]; a13 <- h[, 2]; a14 <- h[, 3]; a23 <- h[, 4]; a24 <- h[, 5]; a34 <- h[, 6]
    r1 <- a12 + a13 + a14; r2 <- a23 + a24; r3 <- a34
    lambda <- max(r1, r2, r3)
    max_row_hazard_increment <- max(max_row_hazard_increment, lambda)
    if (!is.finite(lambda) || lambda <= 0) next
    K <- poisson_terms_needed(lambda, UNIFORMISATION_TOLERANCE)
    max_poisson_terms <- max(max_poisson_terms, K)
    # Sub-stochastic one-jump kernel R = I + Q / lambda (all entries >= 0).
    s12 <- a12 / lambda; s13 <- a13 / lambda; s14 <- a14 / lambda
    s23 <- a23 / lambda; s24 <- a24 / lambda; s34 <- a34 / lambda
    d1 <- 1 - r1 / lambda; d2 <- 1 - r2 / lambda; d3 <- 1 - r3 / lambda
    w <- stats::dpois(0:K, lambda)
    acc <- v * w[1]
    u <- v
    for (j in seq_len(K)) {
      u <- cbind(
        u[, 1] * d1,
        u[, 1] * s12 + u[, 2] * d2,
        u[, 1] * s13 + u[, 3] * d3,
        u[, 2] * s23 + u[, 4] * d3,
        u[, 1] * s14 + u[, 5],
        u[, 2] * s24 + u[, 6],
        u[, 3] * s34 + u[, 7],
        u[, 4] * s34 + u[, 8]
      )
      acc <- acc + u * w[j + 1]
    }
    # Residual Poisson mass (< UNIFORMISATION_TOLERANCE) left in the last iterate.
    v <- acc + u * (1 - sum(w))
    colnames(v) <- S
  }
  p_no <- v[, 1]; p_single <- v[, 2]
  p_cmm <- v[, 3] + v[, 4]
  p_death <- v[, 5] + v[, 6] + v[, 7] + v[, 8]
  cif13 <- v[, 3] + v[, 7]        # entered CMM directly from no CMD
  cif23 <- v[, 4] + v[, 8]        # entered CMM from single CMD
  cif14 <- v[, 5]; cif24 <- v[, 6]
  cif34 <- v[, 7] + v[, 8]
  # Cumulative incidence of the no CMD -> single CMD transition: everyone who
  # ever occupied single CMD is either still there, has moved on to CMM, or has
  # died from it. Baseline single-CMD starters never occupy no CMD, so the
  # transition cannot occur for them.
  cif12 <- if (initial_state == "no_cmd") p_single + cif23 + cif24 else rep(0, n)
  data.frame(
    probability_no_cmd = p_no, probability_single_cmd = p_single,
    probability_cmm_occupancy = p_cmm, probability_death_occupancy = p_death,
    cif_no_to_single = cif12, cif_direct_no_to_cmm = cif13,
    cif_single_to_cmm = cif23, cif_total_cmm_entry = cif13 + cif23,
    cif_death_from_no = cif14, cif_death_from_single = cif24, cif_death_from_cmm = cif34,
    probability_row_sum = p_no + p_single + p_cmm + p_death,
    min_probability_diagnostic = min(v),
    max_row_hazard_increment_diagnostic = max_row_hazard_increment,
    max_poisson_terms_diagnostic = max_poisson_terms,
    uniformisation_tolerance_diagnostic = UNIFORMISATION_TOLERANCE,
    stringsAsFactors = FALSE
  )
}

scenario_probabilities <- function(fits, target_data, cohort_name, replicate_id = 0) {
  hz <- hazard_increment_matrix(fits, horizon)
  result_rows <- list(); diagnostic_rows <- list()
  for (initial_state in c("no_cmd", "single_cmd")) {
    x0 <- target_data[target_data$baseline_state == initial_state, , drop = FALSE]
    if (!nrow(x0)) next
    for (s in seq_len(nrow(scenario_map))) {
      x <- x0
      x$high_mvs <- scenario_map$high_mvs[s]
      x$baseline_social_isolated_2item <- scenario_map$proxy[s]
      x$exposure_x_proxy <- x$high_mvs * x$baseline_social_isolated_2item
      mult <- predict_rate_multipliers(fits, x)
      pred <- propagate_multistate(hz$increments, mult, initial_state)
      metrics <- setdiff(names(pred), c("min_probability_diagnostic", "max_row_hazard_increment_diagnostic",
                                        "max_poisson_terms_diagnostic", "uniformisation_tolerance_diagnostic",
                                        "probability_row_sum"))
      for (metric in metrics) {
        result_rows[[length(result_rows) + 1]] <- data.frame(
          replicate = replicate_id, cohort = cohort_name, initial_state = initial_state,
          scenario = scenario_map$scenario[s], metric = metric,
          probability = mean(pred[[metric]]), per_1000 = 1000 * mean(pred[[metric]]),
          target_persons = nrow(x), stringsAsFactors = FALSE
        )
      }
      diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
        replicate = replicate_id, cohort = cohort_name, initial_state = initial_state,
        scenario = scenario_map$scenario[s], event_times_to_horizon = length(hz$times),
        min_state_probability = min(pred$min_probability_diagnostic),
        max_abs_probability_row_sum_error = max(abs(pred$probability_row_sum - 1)),
        max_row_hazard_increment = max(pred$max_row_hazard_increment_diagnostic),
        uniformisation_tolerance = UNIFORMISATION_TOLERANCE,
        max_poisson_terms = max(pred$max_poisson_terms_diagnostic),
        stringsAsFactors = FALSE
      )
    }
  }
  list(probabilities = do.call(rbind, result_rows), diagnostics = do.call(rbind, diagnostic_rows))
}

interaction_metrics <- function(probabilities, replicate_id = 0) {
  outcome_spec <- data.frame(
    initial_state = c("no_cmd", "no_cmd", "single_cmd"),
    metric = c("cif_direct_no_to_cmm", "cif_total_cmm_entry", "cif_single_to_cmm"),
    outcome = c("direct no CMD to CMM", "any CMM entry from baseline no CMD", "single CMD to CMM"),
    stringsAsFactors = FALSE
  )
  rows <- list()
  for (cohort_name in unique(probabilities$cohort)) {
    for (i in seq_len(nrow(outcome_spec))) {
      x <- probabilities[
        probabilities$cohort == cohort_name & probabilities$initial_state == outcome_spec$initial_state[i] &
          probabilities$metric == outcome_spec$metric[i],
        , drop = FALSE
      ]
      risk <- setNames(x$probability, x$scenario)
      if (!all(c("R00", "R10", "R01", "R11") %in% names(risk)) || risk["R00"] <= 0) next
      rr <- risk / risk["R00"]
      rows[[length(rows) + 1]] <- data.frame(
        replicate = replicate_id, cohort = cohort_name,
        initial_state = outcome_spec$initial_state[i], outcome = outcome_spec$outcome[i],
        risk_R00_per_1000 = 1000 * risk["R00"], risk_R10_per_1000 = 1000 * risk["R10"],
        risk_R01_per_1000 = 1000 * risk["R01"], risk_R11_per_1000 = 1000 * risk["R11"],
        rd_R10_vs_R00_per_1000 = 1000 * (risk["R10"] - risk["R00"]),
        rd_R01_vs_R00_per_1000 = 1000 * (risk["R01"] - risk["R00"]),
        rd_R11_vs_R00_per_1000 = 1000 * (risk["R11"] - risk["R00"]),
        rr_R10_vs_R00 = rr["R10"], rr_R01_vs_R00 = rr["R01"], rr_R11_vs_R00 = rr["R11"],
        additive_joint_excess_per_1000 = 1000 * (risk["R11"] - risk["R10"] - risk["R01"] + risk["R00"]),
        reri = rr["R11"] - rr["R10"] - rr["R01"] + 1,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

# ---- Point estimates --------------------------------------------------------
point_probability_rows <- list(); point_diagnostic_rows <- list(); model_diagnostic_rows <- list()
for (cohort_name in COHORTS) {
  zfit <- fit_transition_models(transitions, cohort_name)
  model_diagnostic_rows[[length(model_diagnostic_rows) + 1]] <- zfit$diagnostics
  if (is.null(zfit$fits)) next
  zpred <- scenario_probabilities(zfit$fits, target[target$cohort == cohort_name, , drop = FALSE], cohort_name, 0)
  point_probability_rows[[length(point_probability_rows) + 1]] <- zpred$probabilities
  point_diagnostic_rows[[length(point_diagnostic_rows) + 1]] <- zpred$diagnostics
}
point_probabilities <- do.call(rbind, point_probability_rows)
point_interaction <- interaction_metrics(point_probabilities, 0)
write_csv(do.call(rbind, model_diagnostic_rows), "gcomp_transition_model_diagnostics.csv")
write_csv(point_probabilities, "five_year_multistate_scenario_probabilities_point.csv")
write_csv(do.call(rbind, point_diagnostic_rows), "five_year_multistate_probability_diagnostics.csv")
write_csv(point_interaction, "five_year_common_population_rd_rr_reri_point.csv")

# ---- Bootstrap (parallel, deterministic per-replicate seeds) ----------------
run_one_bootstrap <- function(b, cohort_name, cohort_target, cohort_transitions, split_rows, seed_base) {
  set.seed(seed_base + b)
  original_ids <- cohort_target$person_id
  sampled_index <- sample(seq_along(original_ids), length(original_ids), replace = TRUE)
  sampled_ids <- original_ids[sampled_index]
  boot_target <- cohort_target[sampled_index, , drop = FALSE]
  boot_target$person_id <- paste0("boot_", seq_len(nrow(boot_target)))
  boot_parts <- vector("list", length(sampled_ids))
  for (j in seq_along(sampled_ids)) {
    rows <- split_rows[[as.character(sampled_ids[j])]]
    piece <- cohort_transitions[rows, , drop = FALSE]
    piece$person_id <- paste0("boot_", j)
    boot_parts[[j]] <- piece
  }
  boot_transitions <- do.call(rbind, boot_parts)
  rownames(boot_transitions) <- NULL
  zfit <- fit_transition_models(boot_transitions, cohort_name)
  if (is.null(zfit$fits)) {
    return(list(status = data.frame(cohort = cohort_name, replicate = b, status = "model_error", n_transition_models_with_warning = NA_integer_,
                                    error = paste(zfit$diagnostics$warning_or_error, collapse = " | "),
                                    stringsAsFactors = FALSE),
                probabilities = NULL, interaction = NULL))
  }
  zpred <- tryCatch(scenario_probabilities(zfit$fits, boot_target, cohort_name, b), error = function(e) e)
  if (inherits(zpred, "error")) {
    return(list(status = data.frame(cohort = cohort_name, replicate = b, status = "prediction_error", n_transition_models_with_warning = NA_integer_,
                                    error = conditionMessage(zpred), stringsAsFactors = FALSE),
                probabilities = NULL, interaction = NULL))
  }
  # A replicate that fits is not necessarily a clean fit: separation warnings on
  # the exposure terms are recorded so that "500/500 successful" cannot be read
  # as "500/500 without warnings".
  warn <- zfit$diagnostics$warning_or_error
  warn <- warn[nzchar(warn)]
  list(status = data.frame(
         cohort = cohort_name, replicate = b,
         status = if (length(warn)) "fit_with_warning" else "fit",
         n_transition_models_with_warning = length(warn),
         error = paste(unique(warn), collapse = " | "), stringsAsFactors = FALSE),
       probabilities = zpred$probabilities, interaction = interaction_metrics(zpred$probabilities, b))
}

bootstrap_probability_rows <- list(); bootstrap_interaction_rows <- list(); bootstrap_status_rows <- list()
if (BOOTSTRAP_B > 0) {
  transition_key_all <- paste(transitions$cohort, transitions$person_id, sep = "__")
  for (ci in seq_along(COHORTS)) {
    cohort_name <- COHORTS[ci]
    cohort_target <- target[target$cohort == cohort_name, , drop = FALSE]
    cohort_transitions <- transitions[transitions$cohort == cohort_name & transition_key_all %in% cohort_target$key, , drop = FALSE]
    split_rows <- split(seq_len(nrow(cohort_transitions)), cohort_transitions$person_id)
    # Integer-safe deterministic seed base (fits in 32-bit): cohort block + date tag.
    seed_base <- ci * 10000000L + 260828L
    message("Bootstrap ", cohort_name, ": B = ", BOOTSTRAP_B, " on ", CORES, " cores...")
    results <- mclapply(
      seq_len(BOOTSTRAP_B),
      function(b) run_one_bootstrap(b, cohort_name, cohort_target, cohort_transitions, split_rows, seed_base),
      mc.cores = CORES, mc.preschedule = TRUE
    )
    for (res in results) {
      if (inherits(res, "try-error") || is.null(res) || is.null(res$status)) {
        bootstrap_status_rows[[length(bootstrap_status_rows) + 1]] <- data.frame(
          cohort = cohort_name, replicate = NA_integer_, status = "worker_error",
          error = paste(as.character(res), collapse = " | "), stringsAsFactors = FALSE)
        next
      }
      bootstrap_status_rows[[length(bootstrap_status_rows) + 1]] <- res$status
      if (!is.null(res$probabilities)) bootstrap_probability_rows[[length(bootstrap_probability_rows) + 1]] <- res$probabilities
      if (!is.null(res$interaction)) bootstrap_interaction_rows[[length(bootstrap_interaction_rows) + 1]] <- res$interaction
    }
  }
  write_csv(do.call(rbind, bootstrap_status_rows), "person_bootstrap_status.csv")

  if (length(bootstrap_probability_rows)) {
    boot_prob <- do.call(rbind, bootstrap_probability_rows)
    prob_ci <- aggregate(
      probability ~ cohort + initial_state + scenario + metric, data = boot_prob,
      FUN = function(z) c(low = quantile(z, .025, na.rm = TRUE), median = median(z, na.rm = TRUE),
                          high = quantile(z, .975, na.rm = TRUE), successful = sum(is.finite(z)))
    )
    prob_ci <- cbind(prob_ci[c("cohort", "initial_state", "scenario", "metric")], as.data.frame(prob_ci$probability))
    point_merge <- point_probabilities[, c("cohort", "initial_state", "scenario", "metric", "probability", "per_1000", "target_persons")]
    prob_ci <- merge(point_merge, prob_ci, by = c("cohort", "initial_state", "scenario", "metric"), all.x = TRUE)
    prob_ci$ci_low_per_1000 <- 1000 * prob_ci$low
    prob_ci$ci_high_per_1000 <- 1000 * prob_ci$high
    write_csv(prob_ci, "five_year_multistate_scenario_probabilities_bootstrap_ci.csv")
  }
  if (length(bootstrap_interaction_rows)) {
    boot_int <- do.call(rbind, bootstrap_interaction_rows)
    value_vars <- c(
      "risk_R00_per_1000", "risk_R10_per_1000", "risk_R01_per_1000", "risk_R11_per_1000",
      "rd_R10_vs_R00_per_1000", "rd_R01_vs_R00_per_1000", "rd_R11_vs_R00_per_1000",
      "rr_R10_vs_R00", "rr_R01_vs_R00", "rr_R11_vs_R00",
      "additive_joint_excess_per_1000", "reri"
    )
    ci_rows <- list()
    for (cohort_name in unique(boot_int$cohort)) {
      for (outcome_name in unique(boot_int$outcome[boot_int$cohort == cohort_name])) {
        x <- boot_int[boot_int$cohort == cohort_name & boot_int$outcome == outcome_name, , drop = FALSE]
        point <- point_interaction[point_interaction$cohort == cohort_name & point_interaction$outcome == outcome_name, , drop = FALSE]
        for (v in value_vars) {
          q <- quantile(x[[v]], c(.025, .5, .975), na.rm = TRUE, names = FALSE)
          ci_rows[[length(ci_rows) + 1]] <- data.frame(
            cohort = cohort_name, outcome = outcome_name, estimand = v,
            point_estimate = point[[v]], bootstrap_ci_low = q[1], bootstrap_median = q[2], bootstrap_ci_high = q[3],
            successful_replicates = sum(is.finite(x[[v]])), requested_replicates = BOOTSTRAP_B,
            elsa_display_note = if (cohort_name == "ELSA") "excluded from main absolute-probability displays" else "",
            stringsAsFactors = FALSE
          )
        }
      }
    }
    write_csv(do.call(rbind, ci_rows), "five_year_common_population_rd_rr_reri_bootstrap_ci.csv")
  }
}

write_session_log("f05_multistate_gcomp", c(
  paste("Bootstrap replicates:", BOOTSTRAP_B, "on", CORES, "cores; per-replicate seeds = cohort_index*1e7 + 260828 + b (core-count independent)."),
  "Transition models use ordinal SES: eliminates documented zero-event education-level separation; MVS estimates insensitive (f04).",
  paste0("Estimator: piecewise-exponential product integral over all six cause-specific hazards; ",
         "each Breslow increment solved exactly by uniformisation on an expanded state space that records ",
         "how each absorbing state was reached (Poisson truncation tolerance ", UNIFORMISATION_TOLERANCE,
         "; diagnostics in five_year_multistate_probability_diagnostics.csv)."),
  "ELSA: mortality ascertainment incomplete; absolute probabilities excluded from main displays."
))
message("f05 complete: ", output_dir)
