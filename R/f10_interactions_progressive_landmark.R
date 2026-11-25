#!/usr/bin/env Rscript
# f10: remaining manuscript-support models on the corrected per-SD scale:
#   (a) progressive adjustment Models 1/2/3 (S11, Results robustness sentence)
#   (b) MVS x age (per 10 years) and MVS x sex interaction models (S13)
#   (c) landmark sensitivity excluding the first 2 years of follow-up (S14)
#   (d) fixed-effect vs DerSimonian-Laird random-effects pooling of Model 3
#       (S15), alongside the REML/modified-HK table from f04.

options(stringsAsFactors = FALSE, width = 170)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f10_interactions_progressive")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions

formula_m1 <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + cluster(person_id)
formula_m2 <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + cluster(person_id)
formula_m3 <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

formula_age_int <- Surv(tstart, tstop, event_target) ~ exposure + exposure:age10 +
  baseline_social_isolated_2item + exposure_x_proxy + baseline_age + baseline_sex +
  factor(baseline_education_level) + factor(baseline_wealth_tertile) +
  baseline_current_smoking + baseline_current_alcohol + baseline_cancer +
  baseline_lung_disease + cluster(person_id)
formula_sex_int <- Surv(tstart, tstop, event_target) ~ exposure + exposure:baseline_sex +
  baseline_social_isolated_2item + exposure_x_proxy + baseline_age + baseline_sex +
  factor(baseline_education_level) + factor(baseline_wealth_tertile) +
  baseline_current_smoking + baseline_current_alcohol + baseline_cancer +
  baseline_lung_disease + cluster(person_id)

extract_term <- function(fit, term_name) {
  b <- coef(fit)[term_name]
  se <- sqrt(diag(vcov(fit)))[term_name]
  c(beta = unname(b), se = unname(se))
}

progressive_rows <- list()
interaction_rows <- list()
landmark_rows <- list()

for (cohort_name in COHORTS) {
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    tr <- TRANSITION_MAP$transition[i]
    fs <- TRANSITION_MAP$from_state[i]
    x <- prepare_model3_riskset(transitions, cohort_name, tr, fs)
    x$age10 <- x$baseline_age / 10

    for (model_name in c("model1_demographic", "model2_ses", "model3_full")) {
      fml <- switch(model_name, model1_demographic = formula_m1, model2_ses = formula_m2, model3_full = formula_m3)
      z <- fit_capture(function() coxph(fml, data = x, ties = "efron", x = TRUE))
      if (is.null(z$fit)) next
      for (term_spec in list(c("exposure", "mvs"), c("baseline_social_isolated_2item", "proxy"))) {
        est <- extract_term(z$fit, term_spec[1])
        progressive_rows[[length(progressive_rows) + 1]] <- data.frame(
          cohort = cohort_name, transition = tr, model = model_name, term = term_spec[2],
          beta = est["beta"], robust_se = est["se"],
          hazard_ratio = exp(est["beta"]),
          ci_low = exp(est["beta"] - 1.96 * est["se"]),
          ci_high = exp(est["beta"] + 1.96 * est["se"]),
          events = sum(x$event_target), warning_or_error = z$warning, stringsAsFactors = FALSE
        )
      }
    }

    for (int_spec in list(c("age", "exposure:age10"), c("sex", "exposure:baseline_sex"))) {
      fml <- if (int_spec[1] == "age") formula_age_int else formula_sex_int
      z <- fit_capture(function() coxph(fml, data = x, ties = "efron", x = TRUE))
      if (is.null(z$fit)) next
      nm <- names(coef(z$fit))
      int_name <- if (int_spec[2] %in% nm) int_spec[2] else {
        alt <- grep(paste0("^", sub("exposure:", "", int_spec[2]), ":exposure$"), nm, value = TRUE)
        if (length(alt)) alt[1] else NA_character_
      }
      if (is.na(int_name)) next
      est <- extract_term(z$fit, int_name)
      interaction_rows[[length(interaction_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = tr, moderator = int_spec[1],
        interaction_beta = est["beta"], interaction_se = est["se"],
        interaction_hr = exp(est["beta"]),
        ci_low = exp(est["beta"] - 1.96 * est["se"]),
        ci_high = exp(est["beta"] + 1.96 * est["se"]),
        p_value = 2 * pnorm(abs(est["beta"] / est["se"]), lower.tail = FALSE),
        events = sum(x$event_target),
        note = "exploratory; not adjusted for multiple testing",
        warning_or_error = z$warning, stringsAsFactors = FALSE
      )
    }

    # Landmark: drop the first two years (events and person-time before year 2).
    xl <- x[x$tstop > 2, , drop = FALSE]
    xl$tstart <- pmax(xl$tstart, 2)
    xl <- xl[xl$tstop > xl$tstart, , drop = FALSE]
    z <- fit_capture(function() coxph(formula_m3, data = xl, ties = "efron", x = TRUE))
    if (!is.null(z$fit)) {
      for (term_spec in list(c("exposure", "mvs"), c("baseline_social_isolated_2item", "proxy"))) {
        est <- extract_term(z$fit, term_spec[1])
        landmark_rows[[length(landmark_rows) + 1]] <- data.frame(
          cohort = cohort_name, transition = tr, term = term_spec[2],
          landmark_years = 2, beta = est["beta"], robust_se = est["se"],
          hazard_ratio = exp(est["beta"]),
          ci_low = exp(est["beta"] - 1.96 * est["se"]),
          ci_high = exp(est["beta"] + 1.96 * est["se"]),
          events = sum(xl$event_target), warning_or_error = z$warning, stringsAsFactors = FALSE
        )
      }
    }
  }
}

write_csv(do.call(rbind, progressive_rows), "progressive_adjustment_models.csv")
write_csv(do.call(rbind, interaction_rows), "mvs_age_sex_interactions.csv")
write_csv(do.call(rbind, landmark_rows), "landmark_two_year_exclusion.csv")

# Pooled interaction terms (REML + modified HK, same conventions as f04).
reml_pool <- function(yi, sei) {
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[keep]; vi <- sei[keep]^2
  k <- length(yi)
  if (k < 2) return(NULL)
  nll <- function(tau2) {
    w <- 1 / (vi + tau2); mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu)^2))
  }
  tau2 <- optimize(nll, interval = c(0, max(1, var(yi) * 100, max(vi) * 100)), tol = 1e-12)$minimum
  if (tau2 < 1e-10) tau2 <- 0
  w <- 1 / (vi + tau2); mu <- sum(w * yi) / sum(w)
  q <- sum(w * (yi - mu)^2) / (k - 1)
  se_m <- sqrt(max(1, q) / sum(w))
  crit <- qt(.975, df = k - 1)
  data.frame(k = k, pooled_hr = exp(mu),
             modified_hk_ci_low = exp(mu - crit * se_m), modified_hk_ci_high = exp(mu + crit * se_m),
             p_value = 2 * pt(abs(mu / se_m), df = k - 1, lower.tail = FALSE), stringsAsFactors = FALSE)
}
ints <- do.call(rbind, interaction_rows)
pool_rows <- list()
for (mod in unique(ints$moderator)) {
  for (tr in unique(ints$transition)) {
    x <- ints[ints$moderator == mod & ints$transition == tr, , drop = FALSE]
    if (grepl("death", tr)) x <- x[x$cohort %in% c("CHARLS", "HRS"), , drop = FALSE]
    z <- reml_pool(x$interaction_beta, x$interaction_se)
    if (!is.null(z)) {
      z$moderator <- mod; z$transition <- tr; z$cohorts <- paste(x$cohort, collapse = "|")
      pool_rows[[length(pool_rows) + 1]] <- z
    }
  }
}
write_csv(do.call(rbind, pool_rows), "pooled_interaction_terms.csv")

# Fixed-effect vs DL random-effects pooling of the corrected Model 3 (S15).
cox_path <- file.path(OUTPUT_ROOT, "f01_mvs_cox", "mvs_variant_cox_results.csv")
if (file.exists(cox_path)) {
  cox <- read.csv(cox_path, check.names = FALSE)
  cox <- cox[cox$variant == "corrected_per_sd" & cox$term %in% c("mvs", "proxy"), , drop = FALSE]
  fe_dl_rows <- list()
  for (term_name in unique(cox$term)) {
    for (tr in unique(cox$transition)) {
      x <- cox[cox$term == term_name & cox$transition == tr, , drop = FALSE]
      if (grepl("death", tr)) x <- x[x$cohort %in% c("CHARLS", "HRS"), , drop = FALSE]
      yi <- x$beta; vi <- x$robust_se^2
      keep <- is.finite(yi) & is.finite(vi) & vi > 0
      yi <- yi[keep]; vi <- vi[keep]
      k <- length(yi)
      if (k < 2) next
      w_fe <- 1 / vi
      mu_fe <- sum(w_fe * yi) / sum(w_fe); se_fe <- sqrt(1 / sum(w_fe))
      Q <- sum(w_fe * (yi - mu_fe)^2)
      tau2_dl <- max(0, (Q - (k - 1)) / (sum(w_fe) - sum(w_fe^2) / sum(w_fe)))
      w_dl <- 1 / (vi + tau2_dl)
      mu_dl <- sum(w_dl * yi) / sum(w_dl); se_dl <- sqrt(1 / sum(w_dl))
      i2 <- max(0, 100 * (Q - (k - 1)) / Q)
      fe_dl_rows[[length(fe_dl_rows) + 1]] <- data.frame(
        term = term_name, transition = tr, k = k, cohorts = paste(x$cohort[keep], collapse = "|"),
        fixed_hr = exp(mu_fe), fixed_ci_low = exp(mu_fe - 1.96 * se_fe), fixed_ci_high = exp(mu_fe + 1.96 * se_fe),
        dl_hr = exp(mu_dl), dl_ci_low = exp(mu_dl - 1.96 * se_dl), dl_ci_high = exp(mu_dl + 1.96 * se_dl),
        tau2_dl = tau2_dl, i2_percent = i2,
        label = if (k == 2) "descriptive two-cohort summary only" else "descriptive; cohort estimates primary",
        stringsAsFactors = FALSE
      )
    }
  }
  write_csv(do.call(rbind, fe_dl_rows), "fixed_vs_dl_pooling.csv")
}

write_session_log("f10_interactions_progressive_landmark", c(
  "Interactions are exploratory and unadjusted for multiple testing.",
  "Landmark sensitivity: risk intervals truncated to start at year 2; earlier events excluded."
))
message("f10 complete: ", output_dir)
