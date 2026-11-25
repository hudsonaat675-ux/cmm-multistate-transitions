#!/usr/bin/env Rscript
# f13: CMD-type heterogeneity.
# (A) no CMD -> first single CMD, split by which disease came first (diabetes,
#     heart disease, stroke): three cause-specific Cox models on the primary
#     Model-3 risk set; the other two first-disease types are censored.
# (B) single CMD -> CMM, split by the index (prevalent) disease.
# Both use the frozen conventions in _common.R (half-month rule, corrected
# per-SD MVS, Model-3 complete cases) and the primary interaction formula.
# Heterogeneity across disease types is tested with a stratified joint model
# (strata = disease type, every coefficient type-specific, cluster(person_id)),
# whose point estimates equal the separate fits; the robust joint covariance
# gives a 2-df Wald test that the exposure coefficient is equal across types.

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f13_cmd_type_heterogeneity")
write_csv <- make_write_csv(output_dir)

PERSON_WAVE_DIR <- env_path("CMM_PERSON_WAVE_DIR")
DISEASES <- c("diabetes", "heart", "stroke")

inputs <- load_frozen_inputs()
transitions <- inputs$transitions

# ---- disease type per person-wave ------------------------------------------
person_wave <- do.call(rbind, lapply(COHORTS, function(cohort_name) {
  p <- read.csv(file.path(PERSON_WAVE_DIR, paste0(tolower(cohort_name), "_person_wave.csv")), check.names = FALSE)
  pattern <- paste0(ifelse(p$diabetes %in% 1, "D", "-"), ifelse(p$heart %in% 1, "H", "-"), ifelse(p$stroke %in% 1, "S", "-"))
  data.frame(key = paste(cohort_name, p$person_id, p$wave, sep = "__"),
             single_type = c("D--" = "diabetes", "-H-" = "heart", "--S" = "stroke")[pattern],
             stringsAsFactors = FALSE)
}))
type_at <- function(cohort, person_id, wave) {
  person_wave$single_type[match(paste(cohort, person_id, wave, sep = "__"), person_wave$key)]
}

formula_terms <- list(
  full_factor = c("exposure", "baseline_social_isolated_2item", "exposure_x_proxy", "baseline_age", "baseline_sex",
                  "factor(baseline_education_level)", "factor(baseline_wealth_tertile)", "baseline_current_smoking",
                  "baseline_current_alcohol", "baseline_cancer", "baseline_lung_disease"),
  ordinal_ses = c("exposure", "baseline_social_isolated_2item", "exposure_x_proxy", "baseline_age", "baseline_sex",
                  "baseline_education_level", "baseline_wealth_tertile", "baseline_current_smoking",
                  "baseline_current_alcohol", "baseline_cancer", "baseline_lung_disease")
)
focal <- c(exposure = "mvs", baseline_social_isolated_2item = "proxy", exposure_x_proxy = "mvs_x_proxy")

term_rows <- function(fit, coef_names, meta) {
  b <- coef(fit)[coef_names]; se <- sqrt(diag(vcov(fit)))[coef_names]
  cbind(meta[rep(1, length(coef_names)), , drop = FALSE],
        data.frame(term = unname(focal[sub("__.*$", "", coef_names)]), beta = unname(b), robust_se = unname(se),
                   hazard_ratio = exp(unname(b)), ci_low = exp(unname(b) - 1.96 * unname(se)),
                   ci_high = exp(unname(b) + 1.96 * unname(se)),
                   p_value = 2 * pnorm(abs(unname(b / se)), lower.tail = FALSE), stringsAsFactors = FALSE))
}

# Joint stratified model: one design block per disease type, so every
# coefficient is type-specific and the partial likelihood factorises by type.
fit_joint <- function(x, rhs_terms) {
  rhs <- paste(rhs_terms, collapse = " + ")
  mm <- model.matrix(as.formula(paste("~", rhs)), data = x)[, -1, drop = FALSE]
  base_names <- make.names(colnames(mm))
  colnames(mm) <- base_names
  blocks <- lapply(DISEASES, function(dz) {
    m <- mm * as.numeric(x$group == dz)
    colnames(m) <- paste0(base_names, "__", dz)
    m
  })
  design <- do.call(cbind, blocks)
  design <- design[, colSums(abs(design)) > 0, drop = FALSE]
  dat <- data.frame(tstart = x$tstart, tstop = x$tstop, event_target = x$event_target,
                    group = x$group, person_id = x$person_id, design, check.names = FALSE)
  f <- as.formula(paste("Surv(tstart, tstop, event_target) ~",
                        paste(sprintf("`%s`", colnames(design)), collapse = " + "),
                        "+ strata(group) + cluster(person_id)"))
  fit_capture(function() coxph(f, data = dat, ties = "efron"))
}

wald_equal <- function(fit, base, groups) {
  nm <- paste0(base, "__", groups)
  nm <- nm[nm %in% names(coef(fit)) & is.finite(coef(fit)[nm])]
  if (length(nm) < 2) return(c(df = NA, chisq = NA, p = NA))
  b <- coef(fit)[nm]; V <- vcov(fit)[nm, nm, drop = FALSE]
  C <- cbind(1, -diag(length(nm) - 1))
  if (length(nm) == 2) C <- matrix(c(1, -1), nrow = 1)
  d <- C %*% b
  stat <- as.numeric(t(d) %*% solve(C %*% V %*% t(C), d))
  c(df = nrow(C), chisq = stat, p = pchisq(stat, nrow(C), lower.tail = FALSE))
}

# ---- build the two analysis data sets per cohort ---------------------------
analysis_sets <- list()
for (cohort_name in COHORTS) {
  a <- prepare_model3_riskset(transitions, cohort_name, "no_cmd_to_single_cmd", "no_cmd")
  first_type <- ifelse(a$event_target == 1, type_at(cohort_name, a$person_id, a$to_wave), NA_character_)
  if (any(a$event_target == 1 & is.na(first_type))) stop("Unlinked first-CMD event in ", cohort_name)
  stacked <- do.call(rbind, lapply(DISEASES, function(dz) {
    y <- a; y$group <- dz
    y$event_target <- as.integer(a$event_target == 1 & first_type == dz)
    y
  }))
  analysis_sets[[length(analysis_sets) + 1]] <- list(cohort = cohort_name, part = "A_first_cmd_type",
    transition = "no_cmd_to_single_cmd", data = stacked, overall_events = sum(a$event_target),
    overall_intervals = nrow(a))

  b <- prepare_model3_riskset(transitions, cohort_name, "single_cmd_to_cmm", "single_cmd")
  b$group <- type_at(cohort_name, b$person_id, b$from_wave)
  if (any(is.na(b$group))) stop("Unlinked index disease in ", cohort_name)
  analysis_sets[[length(analysis_sets) + 1]] <- list(cohort = cohort_name, part = "B_index_cmd_type",
    transition = "single_cmd_to_cmm", data = b, overall_events = sum(b$event_target),
    overall_intervals = nrow(b))
}

# ---- fit ---------------------------------------------------------------------
result_rows <- list(); diag_rows <- list(); het_rows <- list(); count_rows <- list(); check_rows <- list()
for (s in analysis_sets) {
  x <- s$data
  for (dz in DISEASES) {
    g <- x[x$group == dz, , drop = FALSE]
    count_rows[[length(count_rows) + 1]] <- data.frame(
      cohort = s$cohort, part = s$part, transition = s$transition, cmd_type = dz,
      risk_intervals = nrow(g), events = sum(g$event_target), unique_persons = length(unique(g$person_id)),
      overall_intervals = s$overall_intervals, overall_events = s$overall_events, stringsAsFactors = FALSE)
  }
  for (spec in names(formula_terms)) {
    rhs <- formula_terms[[spec]]
    # separate per-type fits (reporting estimates)
    for (dz in DISEASES) {
      g <- x[x$group == dz, , drop = FALSE]
      f <- as.formula(paste("Surv(tstart, tstop, event_target) ~", paste(rhs, collapse = " + "), "+ cluster(person_id)"))
      z <- fit_capture(function() coxph(f, data = g, ties = "efron"))
      meta <- data.frame(cohort = s$cohort, part = s$part, transition = s$transition, cmd_type = dz,
                         specification = spec, risk_intervals = nrow(g), events = sum(g$event_target),
                         unique_persons = length(unique(g$person_id)), stringsAsFactors = FALSE)
      diag_rows[[length(diag_rows) + 1]] <- cbind(meta, status = if (is.null(z$fit)) "error" else "fit",
                                                  warning_or_error = z$warning)
      if (!is.null(z$fit)) result_rows[[length(result_rows) + 1]] <- term_rows(z$fit, names(focal), meta)
      if (!is.null(z$fit) && spec == "full_factor") {
        separate_b <- coef(z$fit)[names(focal)]
        check_rows[[length(check_rows) + 1]] <- data.frame(cohort = s$cohort, part = s$part, cmd_type = dz,
          term = unname(focal), separate_beta = unname(separate_b), stringsAsFactors = FALSE)
      }
    }
    # joint stratified model for the heterogeneity test
    zj <- fit_joint(x, rhs)
    for (base in c("exposure", "baseline_social_isolated_2item")) {
      w <- if (is.null(zj$fit)) c(df = NA, chisq = NA, p = NA) else wald_equal(zj$fit, base, DISEASES)
      het_rows[[length(het_rows) + 1]] <- data.frame(
        cohort = s$cohort, part = s$part, transition = s$transition, specification = spec,
        term = unname(focal[base]), wald_df = w[["df"]], wald_chisq = w[["chisq"]], p_heterogeneity = w[["p"]],
        joint_status = if (is.null(zj$fit)) "error" else "fit", joint_warning = zj$warning, stringsAsFactors = FALSE)
    }
    if (!is.null(zj$fit) && spec == "full_factor") {
      for (dz in DISEASES) for (base in names(focal)) {
        nm <- paste0(base, "__", dz)
        idx <- which(sapply(check_rows, function(r) r$cohort[1] == s$cohort && r$part[1] == s$part && r$cmd_type[1] == dz))
        if (length(idx)) {
          r <- check_rows[[idx]]
          r$joint_beta[r$term == focal[[base]]] <- if (nm %in% names(coef(zj$fit))) coef(zj$fit)[[nm]] else NA_real_
          check_rows[[idx]] <- r
        }
      }
    }
  }
}

results <- do.call(rbind, result_rows)
write_csv(do.call(rbind, count_rows), "cmd_type_riskset_counts.csv")
write_csv(results, "cmd_type_cox_results.csv")
write_csv(do.call(rbind, diag_rows), "cmd_type_cox_diagnostics.csv")
write_csv(do.call(rbind, het_rows), "cmd_type_heterogeneity_tests.csv")
checks <- do.call(rbind, check_rows)
checks$abs_diff <- abs(checks$separate_beta - checks$joint_beta)
write_csv(checks, "joint_vs_separate_fit_check.csv")

# ---- pooled across cohorts (REML + modified Hartung-Knapp, as in f04) -------
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
  data.frame(k = k, pooled_log_hr = mu, pooled_hr = exp(mu), tau2_reml = tau2,
             modified_hk_ci_low = exp(mu - crit * se_modified_hk), modified_hk_ci_high = exp(mu + crit * se_modified_hk),
             inference_label = if (k == 2) "descriptive two-cohort summary only" else "small-k summary; cohort estimates remain primary",
             stringsAsFactors = FALSE)
}
pool_rows <- list()
primary <- results[results$specification == "full_factor" & results$term %in% c("mvs", "proxy"), , drop = FALSE]
for (key in unique(paste(primary$part, primary$cmd_type, primary$term))) {
  x <- primary[paste(primary$part, primary$cmd_type, primary$term) == key, , drop = FALSE]
  z <- reml_hk(x$beta, x$robust_se)
  if (!is.null(z)) pool_rows[[length(pool_rows) + 1]] <- cbind(
    data.frame(part = x$part[1], transition = x$transition[1], cmd_type = x$cmd_type[1], term = x$term[1],
               cohorts = paste(x$cohort, collapse = "|"), total_events = sum(x$events), stringsAsFactors = FALSE), z)
}
write_csv(do.call(rbind, pool_rows), "cmd_type_pooled_reml_hk.csv")

write_session_log("f13_cmd_type_heterogeneity", c(
  "CMD-type-specific transitions on the primary Model-3 risk set.",
  "Disease type linked from cohort person-wave files at to_wave (part A) or from_wave (part B).",
  paste("Person-wave directory:", PERSON_WAVE_DIR)
))
message("f13 complete: ", output_dir)
