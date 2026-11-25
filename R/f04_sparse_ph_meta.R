#!/usr/bin/env Rscript
# f04: separation localisation (level x event cells), ordinal-SES / reduced /
# ridge sensitivity, PH diagnostics (global + term), predefined 5-year piecewise
# effects, and REML meta with MODIFIED Hartung-Knapp as the reporting version
# plus prediction intervals; k=2 mortality pooling labelled descriptive.

options(stringsAsFactors = FALSE, width = 170)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f04_sparse_ph_meta")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions

formula_full <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

formula_ordinal <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + baseline_education_level +
  baseline_wealth_tertile + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

formula_reduced_core <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  cluster(person_id)

formula_reduced_minimal <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + cluster(person_id)

formula_set <- list(
  full_factor = formula_full, ordinal_ses = formula_ordinal,
  reduced_core = formula_reduced_core, reduced_minimal = formula_reduced_minimal
)

categorical_vars <- c(
  "baseline_social_isolated_2item", "baseline_sex", "baseline_education_level",
  "baseline_wealth_tertile", "baseline_current_smoking", "baseline_current_alcohol",
  "baseline_cancer", "baseline_lung_disease"
)

split_at_time <- function(x, cut = 5) {
  before <- x[x$tstart < cut, , drop = FALSE]
  if (nrow(before)) {
    original_stop <- before$tstop
    original_event <- before$event_target
    before$tstop <- pmin(before$tstop, cut)
    before$event_target <- ifelse(original_stop <= cut, original_event, 0L)
    before$late_period <- 0L
  }
  after <- x[x$tstop > cut, , drop = FALSE]
  if (nrow(after)) {
    after$tstart <- pmax(after$tstart, cut)
    after$late_period <- 1L
  }
  out <- rbind(before, after)
  out[out$tstop > out$tstart, , drop = FALSE]
}

piecewise_formula <- Surv(tstart, tstop, event_target) ~
  exposure + baseline_social_isolated_2item + exposure_x_proxy +
  exposure:late_period + baseline_social_isolated_2item:late_period + exposure_x_proxy:late_period +
  baseline_age + baseline_sex + factor(baseline_education_level) + factor(baseline_wealth_tertile) +
  baseline_current_smoking + baseline_current_alcohol + baseline_cancer + baseline_lung_disease +
  cluster(person_id)

level_rows <- list(); model_rows <- list(); diagnostic_rows <- list()
coefficient_diag_rows <- list(); ph_rows <- list(); ph_status_rows <- list()
piecewise_rows <- list(); piecewise_diag_rows <- list()
ridge_rows <- list(); ridge_diag_rows <- list()

for (cohort_name in COHORTS) {
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    transition_name <- TRANSITION_MAP$transition[i]
    from_state_name <- TRANSITION_MAP$from_state[i]
    x <- prepare_model3_riskset(transitions, cohort_name, transition_name, from_state_name)

    for (v in categorical_vars) {
      levs <- sort(unique(x[[v]]))
      for (lev in levs) {
        keep <- x[[v]] == lev
        level_rows[[length(level_rows) + 1]] <- data.frame(
          cohort = cohort_name, transition = transition_name, variable = v, level = as.character(lev),
          risk_intervals = sum(keep), event_intervals = sum(x$event_target[keep]),
          unique_persons = length(unique(x$person_id[keep])),
          zero_event_cell = sum(x$event_target[keep]) == 0,
          stringsAsFactors = FALSE
        )
      }
    }

    for (model_name in names(formula_set)) {
      z <- fit_capture(function() coxph(formula_set[[model_name]], data = x, ties = "efron", x = TRUE))
      diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = transition_name, model = model_name,
        risk_intervals = nrow(x), events = sum(x$event_target),
        status = if (is.null(z$fit)) "error" else "fit", warning_or_error = z$warning,
        stringsAsFactors = FALSE
      )
      if (is.null(z$fit)) next
      b_all <- coef(z$fit); se_all <- sqrt(diag(vcov(z$fit)))
      coefficient_diag_rows[[length(coefficient_diag_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = transition_name, model = model_name,
        coefficient = names(b_all), beta = unname(b_all), standard_error = unname(se_all),
        absolute_beta_over_10 = abs(unname(b_all)) > 10, standard_error_over_10 = unname(se_all) > 10,
        stringsAsFactors = FALSE
      )
      term_map <- c(exposure = "mvs", baseline_social_isolated_2item = "proxy", exposure_x_proxy = "mvs_x_proxy")
      b <- b_all[names(term_map)]; se <- se_all[names(term_map)]
      model_rows[[length(model_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = transition_name, model = model_name,
        term = unname(term_map), beta = unname(b), robust_se = unname(se),
        hazard_ratio = exp(unname(b)), ci_low = exp(unname(b) - 1.96 * unname(se)),
        ci_high = exp(unname(b) + 1.96 * unname(se)),
        p_value = 2 * pnorm(abs(unname(b / se)), lower.tail = FALSE),
        risk_intervals = nrow(x), events = sum(x$event_target), stringsAsFactors = FALSE
      )

      if (model_name == "full_factor") {
        zph_warning <- character(); zph_error <- NULL
        zph <- tryCatch(
          withCallingHandlers(
            cox.zph(z$fit, terms = TRUE),
            warning = function(w) { zph_warning <<- c(zph_warning, conditionMessage(w)); invokeRestart("muffleWarning") }
          ),
          error = function(e) { zph_error <<- conditionMessage(e); NULL }
        )
        ph_status_rows[[length(ph_status_rows) + 1]] <- data.frame(
          cohort = cohort_name, transition = transition_name,
          status = if (is.null(zph)) "error" else "fit",
          warning_or_error = paste(unique(c(z$warning, zph_warning, zph_error)), collapse = " | "),
          stringsAsFactors = FALSE
        )
        if (!is.null(zph)) {
          tab <- as.data.frame(zph$table); tab$term <- rownames(tab); rownames(tab) <- NULL
          ph_rows[[length(ph_rows) + 1]] <- data.frame(
            cohort = cohort_name, transition = transition_name, term = tab$term,
            chisq = tab[, "chisq"], df = tab[, "df"], p_value = tab[, "p"],
            risk_intervals = nrow(x), events = sum(x$event_target), stringsAsFactors = FALSE
          )
        }
      }
    }

    xs <- split_at_time(x, cut = 5)
    zp <- fit_capture(function() coxph(piecewise_formula, data = xs, ties = "efron", x = TRUE))
    piecewise_diag_rows[[length(piecewise_diag_rows) + 1]] <- data.frame(
      cohort = cohort_name, transition = transition_name, cut_years = 5,
      early_events = sum(xs$event_target[xs$late_period == 0]),
      late_events = sum(xs$event_target[xs$late_period == 1]),
      status = if (is.null(zp$fit)) "error" else "fit", warning_or_error = zp$warning,
      stringsAsFactors = FALSE
    )
    if (!is.null(zp$fit)) {
      b <- coef(zp$fit); vv <- vcov(zp$fit)
      specs <- data.frame(
        base = c("exposure", "baseline_social_isolated_2item", "exposure_x_proxy"),
        change = c("exposure:late_period", "baseline_social_isolated_2item:late_period", "exposure_x_proxy:late_period"),
        term = c("mvs", "proxy", "mvs_x_proxy"), stringsAsFactors = FALSE
      )
      for (j in seq_len(nrow(specs))) {
        base_name <- specs$base[j]; change_name <- specs$change[j]
        if (!change_name %in% names(b)) {
          alt <- paste("late_period", base_name, sep = ":")
          if (alt %in% names(b)) change_name <- alt
        }
        if (!base_name %in% names(b) || !change_name %in% names(b)) next
        early_beta <- unname(b[base_name]); early_se <- sqrt(vv[base_name, base_name])
        late_beta <- unname(b[base_name] + b[change_name])
        late_var <- vv[base_name, base_name] + vv[change_name, change_name] + 2 * vv[base_name, change_name]
        change_beta <- unname(b[change_name]); change_se <- sqrt(vv[change_name, change_name])
        piecewise_rows[[length(piecewise_rows) + 1]] <- data.frame(
          cohort = cohort_name, transition = transition_name, term = specs$term[j], cut_years = 5,
          early_hr = exp(early_beta), early_ci_low = exp(early_beta - 1.96 * early_se), early_ci_high = exp(early_beta + 1.96 * early_se),
          late_hr = exp(late_beta), late_ci_low = exp(late_beta - 1.96 * sqrt(late_var)), late_ci_high = exp(late_beta + 1.96 * sqrt(late_var)),
          p_time_interaction = 2 * pnorm(abs(change_beta / change_se), lower.tail = FALSE),
          stringsAsFactors = FALSE
        )
      }
    }

    nuisance <- model.matrix(
      ~ baseline_age + baseline_sex + factor(baseline_education_level) + factor(baseline_wealth_tertile) +
        baseline_current_smoking + baseline_current_alcohol + baseline_cancer + baseline_lung_disease,
      data = x
    )[, -1, drop = FALSE]
    x$nuisance_ridge <- nuisance
    for (theta in c(1, 10)) {
      zr <- fit_capture(function() coxph(
        Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item + exposure_x_proxy +
          ridge(nuisance_ridge, theta = theta),
        data = x, ties = "efron", x = TRUE
      ))
      ridge_diag_rows[[length(ridge_diag_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = transition_name, theta = theta,
        status = if (is.null(zr$fit)) "error" else "fit", warning_or_error = zr$warning,
        stringsAsFactors = FALSE
      )
      if (!is.null(zr$fit)) {
        b <- coef(zr$fit)[c("exposure", "baseline_social_isolated_2item", "exposure_x_proxy")]
        ridge_rows[[length(ridge_rows) + 1]] <- data.frame(
          cohort = cohort_name, transition = transition_name, theta = theta,
          term = c("mvs", "proxy", "mvs_x_proxy"), beta = unname(b), hazard_ratio = exp(unname(b)),
          inference = "point-estimate diagnostic; no cluster-robust CI", stringsAsFactors = FALSE
        )
      }
    }
  }
}

write_csv(do.call(rbind, level_rows), "categorical_level_event_cells.csv")
write_csv(do.call(rbind, model_rows), "sparse_model_sensitivity_results.csv")
write_csv(do.call(rbind, diagnostic_rows), "sparse_model_sensitivity_diagnostics.csv")
write_csv(do.call(rbind, coefficient_diag_rows), "all_coefficient_separation_diagnostics.csv")
write_csv(do.call(rbind, ph_rows), "corrected_mvs_ph_tests.csv")
write_csv(do.call(rbind, ph_status_rows), "corrected_mvs_ph_status.csv")
write_csv(do.call(rbind, piecewise_rows), "five_year_piecewise_exposure_effects.csv")
write_csv(do.call(rbind, piecewise_diag_rows), "five_year_piecewise_diagnostics.csv")
write_csv(do.call(rbind, ridge_rows), "nuisance_ridge_point_estimate_sensitivity.csv")
write_csv(do.call(rbind, ridge_diag_rows), "nuisance_ridge_diagnostics.csv")

# ---- REML meta with modified HK as the reporting version ---------------
reml_hk <- function(yi, sei) {
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[keep]; vi <- sei[keep]^2
  k <- length(yi)
  if (k < 2) return(NULL)
  nll <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu)^2))
  }
  upper <- max(1, stats::var(yi) * 100, max(vi) * 100)
  tau2 <- optimize(nll, interval = c(0, upper), tol = 1e-12)$minimum
  if (tau2 < 1e-10) tau2 <- 0
  w <- 1 / (vi + tau2)
  mu <- sum(w * yi) / sum(w)
  q_scale <- sum(w * (yi - mu)^2) / (k - 1)
  se_modified_hk <- sqrt(max(1, q_scale) / sum(w))
  crit <- qt(.975, df = k - 1)
  pi_crit <- if (k >= 3) qt(.975, df = k - 2) else NA_real_
  data.frame(
    k = k, pooled_log_hr = mu, pooled_hr = exp(mu), tau2_reml = tau2, hk_scale_q = q_scale,
    modified_hk_ci_low = exp(mu - crit * se_modified_hk), modified_hk_ci_high = exp(mu + crit * se_modified_hk),
    prediction_interval_low = if (k >= 3) exp(mu - pi_crit * sqrt(tau2 + 1 / sum(w))) else NA_real_,
    prediction_interval_high = if (k >= 3) exp(mu + pi_crit * sqrt(tau2 + 1 / sum(w))) else NA_real_,
    reporting_note = "modified Hartung-Knapp (q* = max(1, q)) is the reporting version",
    inference_label = if (k == 2) "descriptive two-cohort summary only" else "small-k sensitivity; cohort estimates remain primary",
    stringsAsFactors = FALSE
  )
}

cox_result_path <- file.path(OUTPUT_ROOT, "f01_mvs_cox", "mvs_variant_cox_results.csv")
if (file.exists(cox_result_path)) {
  cox <- read.csv(cox_result_path, check.names = FALSE)
  cox <- cox[cox$variant == "corrected_per_sd" & cox$term %in% c("mvs", "proxy", "mvs_x_proxy"), , drop = FALSE]
  meta_rows <- list()
  for (term_name in unique(cox$term)) {
    for (transition_name in unique(cox$transition)) {
      x <- cox[cox$term == term_name & cox$transition == transition_name, , drop = FALSE]
      if (grepl("death", transition_name)) x <- x[x$cohort %in% c("CHARLS", "HRS"), , drop = FALSE]
      x <- x[is.finite(x$beta) & is.finite(x$robust_se), , drop = FALSE]
      z <- reml_hk(x$beta, x$robust_se)
      if (!is.null(z)) {
        z$term <- term_name; z$transition <- transition_name; z$cohorts <- paste(x$cohort, collapse = "|")
        meta_rows[[length(meta_rows) + 1]] <- z
      }
    }
  }
  meta <- do.call(rbind, meta_rows)
  meta <- meta[, c("term", "transition", "cohorts", setdiff(names(meta), c("term", "transition", "cohorts"))), drop = FALSE]
  write_csv(meta, "reml_modified_hk_meta.csv")
}

write_session_log("f04_sparse_ph_meta", c(
  "Meta reporting version: modified Hartung-Knapp with prediction interval; conventional HK not reported.",
  "Ridge: point-estimate diagnostic only (no cluster-robust variance for penalised coxph)."
))
message("f04 complete: ", output_dir)
