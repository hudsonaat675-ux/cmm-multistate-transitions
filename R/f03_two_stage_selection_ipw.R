#!/usr/bin/env Rscript
# f03: two-stage selection analysis.
#   Stage A: exposure-module availability (MVS with >=3 components AND social
#            proxy observed) among all eligible persons. Descriptive rates plus
#            a demographics-only logistic model whose discrimination is reported;
#            no weighting across this boundary (HRS availability is dominated by
#            the EFTF random half-sample design). Estimand statement: analyses
#            are conditional on the exposure-observed population.
#   Stage B: complete-case selection among exposure-observed persons, modelled
#            with substantive covariates only (no near-deterministic missingness
#            indicators). IPW re-weights complete cases to the exposure-observed
#            population, with predefined truncation (P1-99 primary, P2.5-97.5
#            sensitivity), full weight distributions, ESS, overlap, and SMD
#            balance. IPW-weighted Model 3 Cox fits follow.

options(stringsAsFactors = FALSE, width = 170)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f03_selection_ipw")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
eligible <- inputs$eligible

ess <- function(w) {
  w <- w[is.finite(w) & w > 0]
  if (!length(w)) return(NA_real_)
  sum(w)^2 / sum(w^2)
}
weighted_mean <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(ok) || sum(w[ok]) <= 0) return(NA_real_)
  sum(w[ok] * x[ok]) / sum(w[ok])
}
weighted_var <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  if (sum(ok) < 2 || sum(w[ok]) <= 0) return(NA_real_)
  m <- weighted_mean(x[ok], w[ok])
  sum(w[ok] * (x[ok] - m)^2) / sum(w[ok])
}
weight_summary <- function(w, cohort_name, version) {
  q <- quantile(w, c(0, .01, .05, .5, .95, .99, 1), na.rm = TRUE, names = FALSE)
  data.frame(cohort = cohort_name, weight_version = version, n = sum(is.finite(w)),
             mean = mean(w, na.rm = TRUE), sd = sd(w, na.rm = TRUE),
             min = q[1], p01 = q[2], p05 = q[3], median = q[4], p95 = q[5], p99 = q[6], max = q[7],
             cv = sd(w, na.rm = TRUE) / mean(w, na.rm = TRUE),
             effective_sample_size = ess(w), stringsAsFactors = FALSE)
}
cstat <- function(y, p) {
  ok <- is.finite(p) & !is.na(y)
  y <- y[ok]; p <- p[ok]
  r <- rank(p)
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ---- Selection definitions -------------------------------------------------

eligible$exposure_observed <- as.integer(
  is.finite(as.numeric(eligible$mvs_anthro_bp_z)) & !is.na(eligible$social_isolated_2item)
)

model3_required_transition <- c("baseline_mvs_anthro_bp_z", "baseline_social_isolated_2item", MODEL3_COVARS[!MODEL3_COVARS %in% "baseline_social_isolated_2item"])
valid_model3 <- transitions[
  is.finite(transitions$tstart) & is.finite(transitions$tstop) & transitions$tstop > transitions$tstart &
    complete.cases(transitions[, unique(model3_required_transition), drop = FALSE]),
  , drop = FALSE
]
selected_keys <- unique(paste(valid_model3$cohort, valid_model3$person_id, sep = "__"))
eligible$selected_model3 <- as.integer(eligible$key %in% selected_keys)

stage_counts <- do.call(rbind, lapply(split(eligible, eligible$cohort), function(x) {
  data.frame(
    cohort = unique(x$cohort), eligible_n = nrow(x),
    exposure_observed_n = sum(x$exposure_observed),
    exposure_observed_percent = 100 * mean(x$exposure_observed),
    selected_model3_n = sum(x$selected_model3),
    selected_percent_of_eligible = 100 * mean(x$selected_model3),
    selected_percent_of_exposure_observed = 100 * sum(x$selected_model3) / max(1, sum(x$exposure_observed)),
    stringsAsFactors = FALSE
  )
}))
write_csv(stage_counts, "two_stage_selection_counts.csv")

# Availability of key variables among selected vs excluded.
availability_vars <- c(
  "social_participation", "living_alone", "social_isolated_2item",
  "bmi", "waist_cm", "sbp", "dbp", "mvs_anthro_bp_z",
  "education_level", "wealth_tertile", "current_smoking", "current_alcohol",
  "baseline_cancer", "baseline_lung_disease", "cesd_score"
)
availability_rows <- list()
for (cohort_name in COHORTS) {
  x <- eligible[eligible$cohort == cohort_name, , drop = FALSE]
  for (v in availability_vars) {
    availability_rows[[length(availability_rows) + 1]] <- data.frame(
      cohort = cohort_name, variable = v, eligible_n = nrow(x),
      observed_percent = 100 * mean(!is.na(x[[v]])),
      selected_observed_percent = 100 * mean(!is.na(x[[v]][x$selected_model3 == 1])),
      excluded_observed_percent = 100 * mean(!is.na(x[[v]][x$selected_model3 == 0])),
      stringsAsFactors = FALSE
    )
  }
}
write_csv(do.call(rbind, availability_rows), "baseline_variable_availability.csv")

# ---- Stage A: module availability (design-dominated; no weights) ------------

stage_a_rows <- list()
stage_a_fit_rows <- list()
for (cohort_name in COHORTS) {
  x <- eligible[eligible$cohort == cohort_name, , drop = FALSE]
  fit <- glm(exposure_observed ~ splines::ns(age, df = 3) + sex + factor(state_code),
             family = binomial(), data = x)
  p <- predict(fit, type = "response")
  stage_a_rows[[length(stage_a_rows) + 1]] <- data.frame(
    cohort = cohort_name, eligible_n = nrow(x),
    exposure_observed_percent = 100 * mean(x$exposure_observed),
    demographics_model_c_statistic = cstat(x$exposure_observed, p),
    interpretation = if (cohort_name == "HRS")
      "availability dominated by EFTF random half-sample design; demographics discriminate weakly"
    else "module availability high; residual unavailability described, not weighted",
    stringsAsFactors = FALSE
  )
  coefs <- summary(fit)$coefficients
  stage_a_fit_rows[[length(stage_a_fit_rows) + 1]] <- data.frame(
    cohort = cohort_name, term = rownames(coefs), beta = coefs[, 1],
    standard_error = coefs[, 2], p_value = coefs[, 4], stringsAsFactors = FALSE
  )
}
write_csv(do.call(rbind, stage_a_rows), "stage_a_module_availability_summary.csv")
write_csv(do.call(rbind, stage_a_fit_rows), "stage_a_demographics_model_coefficients.csv")

# ---- Stage B: complete-case selection among exposure-observed ---------------

prep_stage_b <- function(x) {
  impute_fraction <- list()
  for (v in c("living_alone", "bmi", "cesd_score")) {
    value <- as.numeric(x[[v]])
    missing <- !is.finite(value)
    med <- median(value[!missing], na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    value[missing] <- med
    x[[paste0(v, "_imp")]] <- value
    impute_fraction[[v]] <- mean(missing)
  }
  attr(x, "impute_fraction") <- impute_fraction
  x
}

stage_b_formula <- selected_model3 ~ splines::ns(age, df = 3) + sex + factor(state_code) +
  living_alone_imp + bmi_imp + cesd_score_imp

balance_vars <- c("age", "sex", "state_code", "living_alone", "bmi", "cesd_score",
                  "education_level", "wealth_tertile", "current_smoking", "current_alcohol",
                  "baseline_cancer", "baseline_lung_disease", "social_isolated_2item",
                  "mvs_anthro_bp_z")

# Reported both as the ordinal/numeric coding above and as level indicators.
CATEGORICAL_BALANCE_VARS <- c("state_code", "education_level", "wealth_tertile")

stage_b_model_rows <- list()
stage_b_impute_rows <- list()
weight_rows <- list()
overlap_rows <- list()
balance_rows <- list()
weighted_person_parts <- list()

for (cohort_name in COHORTS) {
  x <- eligible[eligible$cohort == cohort_name & eligible$exposure_observed == 1, , drop = FALSE]
  x <- prep_stage_b(x)
  frac <- attr(x, "impute_fraction")
  stage_b_impute_rows[[length(stage_b_impute_rows) + 1]] <- data.frame(
    cohort = cohort_name, variable = names(frac),
    median_imputed_percent = 100 * unlist(frac), stringsAsFactors = FALSE
  )

  fit <- glm(stage_b_formula, family = binomial(), data = x)
  p <- pmin(pmax(predict(fit, type = "response"), 1e-6), 1 - 1e-6)
  pbar <- mean(x$selected_model3)
  x$selection_probability <- p
  x$ipw_untruncated <- ifelse(x$selected_model3 == 1, pbar / p, NA_real_)
  selected_weight <- x$ipw_untruncated[x$selected_model3 == 1]
  q01 <- quantile(selected_weight, c(.01, .99), names = FALSE)
  q025 <- quantile(selected_weight, c(.025, .975), names = FALSE)
  x$ipw_p01_p99 <- ifelse(x$selected_model3 == 1, pmin(pmax(x$ipw_untruncated, q01[1]), q01[2]), NA_real_)
  x$ipw_p025_p975 <- ifelse(x$selected_model3 == 1, pmin(pmax(x$ipw_untruncated, q025[1]), q025[2]), NA_real_)

  coefs <- summary(fit)$coefficients
  stage_b_model_rows[[length(stage_b_model_rows) + 1]] <- data.frame(
    cohort = cohort_name, term = rownames(coefs), beta = coefs[, 1], standard_error = coefs[, 2],
    odds_ratio = exp(coefs[, 1]), ci_low = exp(coefs[, 1] - 1.96 * coefs[, 2]),
    ci_high = exp(coefs[, 1] + 1.96 * coefs[, 2]), p_value = coefs[, 4],
    exposure_observed_n = nrow(x), selected_n = sum(x$selected_model3),
    c_statistic = cstat(x$selected_model3, p), stringsAsFactors = FALSE
  )

  for (version in c("ipw_untruncated", "ipw_p01_p99", "ipw_p025_p975")) {
    weight_rows[[length(weight_rows) + 1]] <- weight_summary(x[[version]][x$selected_model3 == 1], cohort_name, version)
  }

  for (group_name in c("exposure_observed", "selected", "not_selected")) {
    keep <- switch(group_name,
                   exposure_observed = rep(TRUE, nrow(x)),
                   selected = x$selected_model3 == 1,
                   not_selected = x$selected_model3 == 0)
    pp <- p[keep]
    qq <- quantile(pp, c(0, .01, .05, .5, .95, .99, 1), na.rm = TRUE, names = FALSE)
    overlap_rows[[length(overlap_rows) + 1]] <- data.frame(
      cohort = cohort_name, group = group_name, n = sum(keep),
      min = qq[1], p01 = qq[2], p05 = qq[3], median = qq[4], p95 = qq[5], p99 = qq[6], max = qq[7],
      probability_below_0_05_percent = 100 * mean(pp < .05),
      probability_below_0_10_percent = 100 * mean(pp < .10), stringsAsFactors = FALSE
    )
  }

  # Balance: complete cases (weighted) versus the exposure-observed target.
  # Multi-category covariates (baseline state, education, wealth) are expanded
  # into level indicators as well: an ordinal mean can sit near the target while
  # individual categories are imbalanced in opposite directions.
  target_w <- rep(1, nrow(x))
  balance_terms <- lapply(balance_vars, function(v) list(label = v, value = as.numeric(x[[v]])))
  for (v in CATEGORICAL_BALANCE_VARS) {
    lv <- sort(unique(x[[v]][!is.na(x[[v]])]))
    for (l in lv) balance_terms[[length(balance_terms) + 1]] <- list(
      label = sprintf("%s [level %s]", v, l), value = as.numeric(x[[v]] == l))
  }
  for (bt in balance_terms) {
    v <- bt$label
    value <- bt$value
    target_mean <- weighted_mean(value, target_w)
    target_var <- weighted_var(value, target_w)
    for (version in c("unweighted_selected", "ipw_untruncated", "ipw_p01_p99", "ipw_p025_p975")) {
      keep <- x$selected_model3 == 1
      w <- if (version == "unweighted_selected") rep(1, sum(keep)) else x[[version]][keep]
      m <- weighted_mean(value[keep], w)
      vv <- weighted_var(value[keep], w)
      denom <- sqrt((target_var + vv) / 2)
      balance_rows[[length(balance_rows) + 1]] <- data.frame(
        cohort = cohort_name, variable = v, comparison = version,
        exposure_observed_target_mean = target_mean, selected_mean = m,
        standardized_mean_difference = if (is.finite(denom) && denom > 0) (m - target_mean) / denom else NA_real_,
        absolute_smd = if (is.finite(denom) && denom > 0) abs((m - target_mean) / denom) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  weighted_person_parts[[length(weighted_person_parts) + 1]] <- x[x$selected_model3 == 1, c(
    "cohort", "person_id", "selection_probability", "ipw_untruncated", "ipw_p01_p99", "ipw_p025_p975"
  ), drop = FALSE]
}

write_csv(do.call(rbind, stage_b_model_rows), "stage_b_selection_model_coefficients.csv")
write_csv(do.call(rbind, stage_b_impute_rows), "stage_b_predictor_median_imputation.csv")
write_csv(do.call(rbind, weight_rows), "selection_weight_distribution_and_ess.csv")
write_csv(do.call(rbind, overlap_rows), "selection_probability_overlap.csv")
balance <- do.call(rbind, balance_rows)
write_csv(balance, "selection_balance_smd.csv")
balance_summary <- aggregate(
  absolute_smd ~ cohort + comparison, data = balance,
  FUN = function(z) c(max = max(z, na.rm = TRUE), mean = mean(z, na.rm = TRUE), over_0_1 = sum(z > .1, na.rm = TRUE))
)
balance_summary <- cbind(balance_summary[c("cohort", "comparison")], as.data.frame(balance_summary$absolute_smd))
write_csv(balance_summary, "selection_balance_summary.csv")

# Level-wise diagnostic reported separately so the "all |SMD| below x" sentence in
# the manuscript can be checked against category indicators, not ordinal means.
balance$term_type <- ifelse(grepl("\\[level ", balance$variable), "category indicator", "numeric/ordinal coding")
balance_by_type <- aggregate(
  absolute_smd ~ cohort + comparison + term_type, data = balance,
  FUN = function(z) c(max = max(z, na.rm = TRUE), mean = mean(z, na.rm = TRUE),
                      over_0_03 = sum(z > .03, na.rm = TRUE), over_0_1 = sum(z > .1, na.rm = TRUE)))
balance_by_type <- cbind(balance_by_type[c("cohort", "comparison", "term_type")],
                         as.data.frame(balance_by_type$absolute_smd))
write_csv(balance_by_type, "selection_balance_summary_by_term_type.csv")

# ---- IPW-weighted Model 3 Cox fits -----------------------------------------

person_weights <- do.call(rbind, weighted_person_parts)
pw_key <- paste(person_weights$cohort, person_weights$person_id, sep = "__")
tr_key <- paste(transitions$cohort, transitions$person_id, sep = "__")
pw_index <- match(tr_key, pw_key)
for (v in c("ipw_untruncated", "ipw_p01_p99", "ipw_p025_p975")) {
  transitions[[v]] <- person_weights[[v]][pw_index]
}

cox_formula <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex +
  factor(baseline_education_level) + factor(baseline_wealth_tertile) +
  baseline_current_smoking + baseline_current_alcohol + baseline_cancer +
  baseline_lung_disease + cluster(person_id)

ipw_result_rows <- list()
ipw_diagnostic_rows <- list()
for (cohort_name in COHORTS) {
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    transition_name <- TRANSITION_MAP$transition[i]
    from_state_name <- TRANSITION_MAP$from_state[i]
    for (version in c("unweighted", "ipw_untruncated", "ipw_p01_p99", "ipw_p025_p975")) {
      x <- prepare_model3_riskset(transitions, cohort_name, transition_name, from_state_name)
      x$analysis_weight <- if (version == "unweighted") 1 else x[[version]]
      x <- x[is.finite(x$analysis_weight), , drop = FALSE]
      z <- fit_capture(function() coxph(cox_formula, data = x, weights = analysis_weight, ties = "efron", x = TRUE))
      ipw_diagnostic_rows[[length(ipw_diagnostic_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = transition_name, weight_version = version,
        risk_intervals = nrow(x), events = sum(x$event_target), unique_persons = length(unique(x$person_id)),
        status = if (is.null(z$fit)) "error" else "fit", warning_or_error = z$warning, stringsAsFactors = FALSE
      )
      if (is.null(z$fit)) next
      term_map <- c(exposure = "mvs", baseline_social_isolated_2item = "proxy", exposure_x_proxy = "mvs_x_proxy")
      b <- coef(z$fit)[names(term_map)]
      se <- sqrt(diag(vcov(z$fit)))[names(term_map)]
      ipw_result_rows[[length(ipw_result_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = transition_name, weight_version = version,
        term = unname(term_map), beta = unname(b), robust_se = unname(se),
        hazard_ratio = exp(unname(b)), ci_low = exp(unname(b) - 1.96 * unname(se)),
        ci_high = exp(unname(b) + 1.96 * unname(se)),
        p_value = 2 * pnorm(abs(unname(b / se)), lower.tail = FALSE),
        risk_intervals = nrow(x), events = sum(x$event_target),
        unique_persons = length(unique(x$person_id)), stringsAsFactors = FALSE
      )
    }
  }
}
write_csv(do.call(rbind, ipw_result_rows), "selection_ipw_cox_results.csv")
write_csv(do.call(rbind, ipw_diagnostic_rows), "selection_ipw_cox_diagnostics.csv")

write_session_log("f03_two_stage_selection_ipw", c(
  "Stage A: exposure-module availability described, not weighted (estimand conditional on exposure-observed population).",
  "Stage B: complete-case selection among exposure-observed persons; predictors are substantive covariates only, no missingness indicators.",
  "Truncation predefined: P1-P99 primary, P2.5-P97.5 sensitivity."
))
message("f03 complete: ", output_dir)
