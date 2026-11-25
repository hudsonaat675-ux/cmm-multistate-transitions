#!/usr/bin/env Rscript
# f01: primary corrected per-SD Cox models (Model 3) plus MVS construction
# sensitivity (complete-4, two-item scores), correlation/PCA structure, and a
# supplementary main-effects model without the MVS x proxy product term.

options(stringsAsFactors = FALSE, width = 160)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f01_mvs_cox")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
eligible <- inputs$eligible

components <- c("bmi", "waist_cm", "sbp", "dbp")
score_parts <- list()
cor_rows <- list()
pca_load_rows <- list()
pca_var_rows <- list()

for (cohort_name in COHORTS) {
  x <- eligible[eligible$cohort == cohort_name, , drop = FALSE]
  zmat <- sapply(components, function(v) zscore(as.numeric(x[[v]])))
  colnames(zmat) <- components
  ncomp <- rowSums(is.finite(zmat))

  complete4_raw <- rowMeans(zmat)
  complete4_raw[ncomp < 4] <- NA_real_

  make_two <- function(a, b) {
    zz <- rowMeans(zmat[, c(a, b), drop = FALSE])
    zz[!is.finite(zmat[, a]) | !is.finite(zmat[, b])] <- NA_real_
    zscore(zz)
  }

  score_parts[[length(score_parts) + 1]] <- data.frame(
    cohort = cohort_name,
    person_id = x$person_id,
    current_raw_unit = x$mvs_anthro_bp_z,
    corrected_per_sd = zscore(as.numeric(x$mvs_anthro_bp_z)),
    complete4_per_sd = zscore(complete4_raw),
    waist_sbp_per_sd = make_two("waist_cm", "sbp"),
    bmi_sbp_per_sd = make_two("bmi", "sbp"),
    waist_dbp_per_sd = make_two("waist_cm", "dbp"),
    bmi_dbp_per_sd = make_two("bmi", "dbp"),
    component_count = ncomp,
    stringsAsFactors = FALSE
  )

  cm <- cor(x[, components, drop = FALSE], use = "pairwise.complete.obs")
  ij <- which(upper.tri(cm, diag = TRUE), arr.ind = TRUE)
  cor_rows[[length(cor_rows) + 1]] <- data.frame(
    cohort = cohort_name, component_1 = rownames(cm)[ij[, 1]],
    component_2 = colnames(cm)[ij[, 2]], correlation = cm[ij], stringsAsFactors = FALSE
  )

  complete_x <- x[complete.cases(x[, components, drop = FALSE]), components, drop = FALSE]
  pc <- prcomp(complete_x, center = TRUE, scale. = TRUE)
  load <- pc$rotation[, 1]
  if (sum(load) < 0) load <- -load
  pca_load_rows[[length(pca_load_rows) + 1]] <- data.frame(
    cohort = cohort_name, component = names(load), pc1_loading = unname(load),
    complete4_n = nrow(complete_x), stringsAsFactors = FALSE
  )
  ve <- pc$sdev^2 / sum(pc$sdev^2)
  pca_var_rows[[length(pca_var_rows) + 1]] <- data.frame(
    cohort = cohort_name, component = paste0("PC", seq_along(ve)),
    variance_explained = ve, cumulative_variance = cumsum(ve),
    complete4_n = nrow(complete_x), stringsAsFactors = FALSE
  )
}

scores <- do.call(rbind, score_parts)
write_csv(do.call(rbind, cor_rows), "mvs_component_correlations.csv")
write_csv(do.call(rbind, pca_load_rows), "mvs_pca_pc1_loadings.csv")
write_csv(do.call(rbind, pca_var_rows), "mvs_pca_variance_explained.csv")

variant_names <- c("current_raw_unit", "corrected_per_sd", "complete4_per_sd",
                   "waist_sbp_per_sd", "bmi_sbp_per_sd", "waist_dbp_per_sd", "bmi_dbp_per_sd")

score_key <- paste(scores$cohort, scores$person_id, sep = "__")
transition_key <- paste(transitions$cohort, transitions$person_id, sep = "__")
score_index <- match(transition_key, score_key)
for (v in variant_names) transitions[[v]] <- scores[[v]][score_index]

formula_interaction <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

formula_main_effects <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

extract_terms <- function(fit, keep_map, meta) {
  keep <- intersect(names(keep_map), names(coef(fit)))
  b <- coef(fit)[keep]
  se <- sqrt(diag(vcov(fit)))[keep]
  cbind(meta[rep(1, length(keep)), , drop = FALSE],
        data.frame(term = unname(keep_map[keep]), beta = unname(b), robust_se = unname(se),
                   hazard_ratio = exp(unname(b)),
                   ci_low = exp(unname(b) - 1.96 * unname(se)),
                   ci_high = exp(unname(b) + 1.96 * unname(se)),
                   p_value = 2 * pnorm(abs(unname(b / se)), lower.tail = FALSE),
                   stringsAsFactors = FALSE))
}

fit_rows <- list()
diag_rows <- list()
main_effect_rows <- list()

for (cohort_name in COHORTS) {
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    tr <- TRANSITION_MAP$transition[i]
    fs <- TRANSITION_MAP$from_state[i]
    for (variant_name in variant_names) {
      x <- prepare_model3_riskset(transitions, cohort_name, tr, fs, exposure_var = variant_name)
      z <- fit_capture(function() coxph(formula_interaction, data = x, ties = "efron", x = TRUE))
      diag_rows[[length(diag_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = tr, variant = variant_name, model = "interaction",
        risk_intervals = nrow(x), events = sum(x$event_target),
        unique_persons = length(unique(x$person_id)),
        status = if (is.null(z$fit)) "error" else "fit", warning_or_error = z$warning,
        stringsAsFactors = FALSE
      )
      if (!is.null(z$fit)) {
        meta <- data.frame(cohort = cohort_name, transition = tr, from_state = fs,
                           variant = variant_name, risk_intervals = nrow(x),
                           events = sum(x$event_target),
                           unique_persons = length(unique(x$person_id)), stringsAsFactors = FALSE)
        fit_rows[[length(fit_rows) + 1]] <- extract_terms(
          z$fit,
          c(exposure = "mvs", baseline_social_isolated_2item = "proxy", exposure_x_proxy = "mvs_x_proxy"),
          meta
        )
      }

      if (variant_name == "corrected_per_sd") {
        zm <- fit_capture(function() coxph(formula_main_effects, data = x, ties = "efron", x = TRUE))
        diag_rows[[length(diag_rows) + 1]] <- data.frame(
          cohort = cohort_name, transition = tr, variant = variant_name, model = "main_effects",
          risk_intervals = nrow(x), events = sum(x$event_target),
          unique_persons = length(unique(x$person_id)),
          status = if (is.null(zm$fit)) "error" else "fit", warning_or_error = zm$warning,
          stringsAsFactors = FALSE
        )
        if (!is.null(zm$fit)) {
          meta <- data.frame(cohort = cohort_name, transition = tr, from_state = fs,
                             variant = variant_name, risk_intervals = nrow(x),
                             events = sum(x$event_target),
                             unique_persons = length(unique(x$person_id)), stringsAsFactors = FALSE)
          main_effect_rows[[length(main_effect_rows) + 1]] <- extract_terms(
            zm$fit,
            c(exposure = "mvs", baseline_social_isolated_2item = "proxy"),
            meta
          )
        }
      }
    }
  }
}

cox_results <- do.call(rbind, fit_rows)
write_csv(cox_results, "mvs_variant_cox_results.csv")
write_csv(do.call(rbind, diag_rows), "mvs_variant_cox_diagnostics.csv")
write_csv(do.call(rbind, main_effect_rows), "mvs_main_effects_no_interaction_results.csv")

headline <- cox_results[
  cox_results$term == "mvs" &
    cox_results$variant %in% c("current_raw_unit", "corrected_per_sd", "complete4_per_sd", "waist_sbp_per_sd") &
    cox_results$transition %in% c("no_cmd_to_single_cmd", "no_cmd_to_cmm", "single_cmd_to_cmm"),
  , drop = FALSE
]
write_csv(headline, "headline_mvs_sensitivity_comparison.csv")

write_session_log("f01_mvs_cox_primary", c(
  "Primary exposure: corrected per-SD MVS (z-score of composite within cohort eligible sample).",
  "Model 3 retains the MVS x proxy product term; a no-interaction main-effects model is exported as supplementary sensitivity."
))
message("f01 complete: ", output_dir)
