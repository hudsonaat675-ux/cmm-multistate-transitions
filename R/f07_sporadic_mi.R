#!/usr/bin/env Rscript
# f07: multiple imputation restricted to sporadically missing covariates among
# exposure-observed persons (structural module missingness is NOT imputed).
# Formal-run change: m = 20 imputations.

options(stringsAsFactors = FALSE, width = 170)
for (pkg in c("survival", "mice")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Package ", pkg, " required")
suppressPackageStartupMessages(library(survival))
source("_common.R")

output_dir <- make_output_dir("f07_sporadic_mi")
write_csv <- make_write_csv(output_dir)

M_IMP <- as.integer(Sys.getenv("FORMAL_MI_M", "20"))
MAXIT <- as.integer(Sys.getenv("FORMAL_MI_MAXIT", "10"))

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
eligible <- inputs$eligible

valid <- is.finite(transitions$tstart) & is.finite(transitions$tstop) & transitions$tstop > transitions$tstart
valid_transitions <- transitions[valid, , drop = FALSE]
valid_keys <- unique(paste(valid_transitions$cohort, valid_transitions$person_id, sep = "__"))
eligible <- eligible[eligible$key %in% valid_keys, , drop = FALSE]

aux_parts <- list()
for (cohort_name in COHORTS) {
  x <- valid_transitions[valid_transitions$cohort == cohort_name, , drop = FALSE]
  ids <- unique(x$person_id)
  aux <- data.frame(cohort = cohort_name, person_id = ids, stringsAsFactors = FALSE)
  aux$max_followup <- vapply(ids, function(id) max(x$tstop[x$person_id == id], na.rm = TRUE), numeric(1))
  for (tr in TRANSITION_MAP$transition) {
    aux[[paste0("event_", tr)]] <- vapply(ids, function(id) {
      as.integer(any(x$person_id == id & x$event == 1 & x$transition == tr))
    }, integer(1))
  }
  aux_parts[[length(aux_parts) + 1]] <- aux
}
aux <- do.call(rbind, aux_parts)
aux_key <- paste(aux$cohort, aux$person_id, sep = "__")

mi_base <- data.frame(
  cohort = eligible$cohort, person_id = eligible$person_id, baseline_state_code = eligible$state_code,
  age = eligible$age, sex = eligible$sex, corrected_per_sd = eligible$corrected_per_sd,
  social_proxy = eligible$social_isolated_2item,
  education = eligible$education_level, wealth = eligible$wealth_tertile,
  smoking = eligible$current_smoking, alcohol = eligible$current_alcohol,
  cancer = eligible$baseline_cancer, lung = eligible$baseline_lung_disease,
  bmi_aux = eligible$bmi, cesd_aux = eligible$cesd_score,
  stringsAsFactors = FALSE
)
mi_key <- paste(mi_base$cohort, mi_base$person_id, sep = "__")
aux_index <- match(mi_key, aux_key)
for (v in setdiff(names(aux), c("cohort", "person_id"))) mi_base[[v]] <- aux[[v]][aux_index]

exposure_complete <- is.finite(mi_base$corrected_per_sd) & is.finite(mi_base$social_proxy)
mi_base <- mi_base[exposure_complete, , drop = FALSE]
sporadic_covars <- c("education", "wealth", "smoking", "alcohol", "cancer", "lung")

scope_rows <- do.call(rbind, lapply(split(mi_base, mi_base$cohort), function(x) {
  cc <- complete.cases(x[, sporadic_covars, drop = FALSE])
  data.frame(
    cohort = unique(x$cohort),
    valid_followup_eligible_n = sum(eligible$cohort == unique(x$cohort)),
    exposure_observed_n = nrow(x),
    sporadic_covariate_complete_n = sum(cc),
    recoverable_by_covariate_mi_n = sum(!cc),
    stringsAsFactors = FALSE
  )
}))
write_csv(scope_rows, "mi_scope_structural_vs_sporadic_missingness.csv")

mi_diagnostic_rows <- list()
imputed_parts <- vector("list", M_IMP)
for (j in seq_len(M_IMP)) imputed_parts[[j]] <- list()

for (cohort_name in COHORTS) {
  x <- mi_base[mi_base$cohort == cohort_name, , drop = FALSE]
  x$education <- ordered(x$education)
  x$wealth <- ordered(x$wealth)
  for (v in c("smoking", "alcohol", "cancer", "lung")) x[[v]] <- factor(x[[v]], levels = c(0, 1))

  # The maxit = 0 call only builds the method/predictor templates, but it emits
  # its own logged events; capture them separately so the diagnostics table
  # cannot be read as "the formal imputation logged nothing" when the log did.
  init_warning <- character()
  ini <- withCallingHandlers(
    mice::mice(x, m = 1, maxit = 0, printFlag = FALSE),
    warning = function(w) { init_warning <<- c(init_warning, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  init_logged <- ini$loggedEvents
  meth <- ini$method
  meth[] <- ""
  for (v in intersect(c("education", "wealth"), names(x))) if (anyNA(x[[v]])) meth[v] <- "polr"
  for (v in intersect(c("smoking", "alcohol", "cancer", "lung"), names(x))) if (anyNA(x[[v]])) meth[v] <- "logreg"
  for (v in intersect(c("bmi_aux", "cesd_aux"), names(x))) if (anyNA(x[[v]])) meth[v] <- "pmm"
  pred <- ini$predictorMatrix
  pred[, c("cohort", "person_id")] <- 0
  pred[c("cohort", "person_id"), ] <- 0
  warning_text <- character()
  imp <- withCallingHandlers(
    mice::mice(x, m = M_IMP, maxit = MAXIT, method = meth, predictorMatrix = pred,
               seed = GLOBAL_SEED + match(cohort_name, COHORTS), printFlag = FALSE),
    warning = function(w) { warning_text <<- c(warning_text, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  logged <- imp$loggedEvents
  mi_diagnostic_rows[[length(mi_diagnostic_rows) + 1]] <- data.frame(
    cohort = cohort_name, exposure_observed_n = nrow(x), imputations = M_IMP, iterations = MAXIT,
    variables_imputed = paste(names(meth)[nzchar(meth)], collapse = "|"),
    logged_event_n = if (is.null(logged)) 0 else nrow(logged),
    initialisation_logged_event_n = if (is.null(init_logged)) 0 else nrow(init_logged),
    initialisation_logged_detail = if (is.null(init_logged)) "" else
      paste(unique(paste(init_logged$meth, init_logged$out, sep = ":")), collapse = " | "),
    initialisation_warning = paste(unique(init_warning), collapse = " | "),
    warning = paste(unique(warning_text), collapse = " | "), stringsAsFactors = FALSE
  )
  for (j in seq_len(M_IMP)) {
    completed <- mice::complete(imp, action = j)
    for (v in sporadic_covars) completed[[v]] <- as.numeric(as.character(completed[[v]]))
    imputed_parts[[j]][[length(imputed_parts[[j]]) + 1]] <- completed
  }
}
write_csv(do.call(rbind, mi_diagnostic_rows), "mice_run_diagnostics.csv")

fit_formula <- Surv(tstart, tstop, event_target) ~ corrected_per_sd + social_proxy +
  mvs_x_proxy + age + sex + education + wealth + smoking + alcohol + cancer + lung + cluster(person_id)

estimate_rows <- list(); fit_diagnostic_rows <- list()
for (j in seq_len(M_IMP)) {
  completed <- do.call(rbind, imputed_parts[[j]])
  completed_key <- paste(completed$cohort, completed$person_id, sep = "__")
  tkey <- paste(valid_transitions$cohort, valid_transitions$person_id, sep = "__")
  index <- match(tkey, completed_key)
  td <- valid_transitions[is.finite(index), , drop = FALSE]
  index <- index[is.finite(index)]
  copy_vars <- c("corrected_per_sd", "social_proxy", "age", "sex", "education", "wealth", "smoking", "alcohol", "cancer", "lung")
  for (v in copy_vars) td[[v]] <- completed[[v]][index]
  td$mvs_x_proxy <- td$corrected_per_sd * td$social_proxy

  for (cohort_name in COHORTS) {
    for (i in seq_len(nrow(TRANSITION_MAP))) {
      tr <- TRANSITION_MAP$transition[i]
      from <- TRANSITION_MAP$from_state[i]
      x <- td[td$cohort == cohort_name & td$from_state == from, , drop = FALSE]
      x$event_target <- as.integer(x$event == 1 & x$transition == tr)
      z <- fit_capture(function() coxph(fit_formula, data = x, ties = "efron", x = TRUE))
      fit_diagnostic_rows[[length(fit_diagnostic_rows) + 1]] <- data.frame(
        imputation = j, cohort = cohort_name, transition = tr,
        risk_intervals = nrow(x), events = sum(x$event_target),
        status = if (is.null(z$fit)) "error" else "fit", warning_or_error = z$warning,
        stringsAsFactors = FALSE
      )
      if (is.null(z$fit)) next
      term_map <- c(corrected_per_sd = "mvs", social_proxy = "proxy", mvs_x_proxy = "mvs_x_proxy")
      b <- coef(z$fit)[names(term_map)]
      vv <- vcov(z$fit)
      estimate_rows[[length(estimate_rows) + 1]] <- data.frame(
        imputation = j, cohort = cohort_name, transition = tr,
        term = unname(term_map), beta = unname(b), within_variance = diag(vv)[names(term_map)],
        stringsAsFactors = FALSE
      )
    }
  }
}

estimates <- do.call(rbind, estimate_rows)
write_csv(do.call(rbind, fit_diagnostic_rows), "mi_cox_fit_diagnostics.csv")
estimates <- estimates[is.finite(estimates$beta) & is.finite(estimates$within_variance) & estimates$within_variance >= 0, , drop = FALSE]

pool_rows <- list()
groups <- split(estimates, paste(estimates$cohort, estimates$transition, estimates$term, sep = "__"))
for (x in groups) {
  m_eff <- nrow(x)
  qbar <- mean(x$beta)
  ubar <- mean(x$within_variance)
  bvar <- if (m_eff > 1) var(x$beta) else 0
  total <- ubar + (1 + 1 / m_eff) * bvar
  r <- if (ubar > 0) ((1 + 1 / m_eff) * bvar) / ubar else Inf
  df <- if (is.finite(r) && r > 0) (m_eff - 1) * (1 + 1 / r)^2 else 1e6
  crit <- qt(.975, df = df)
  pool_rows[[length(pool_rows) + 1]] <- data.frame(
    cohort = x$cohort[1], transition = x$transition[1], term = x$term[1],
    imputations_successful = m_eff, beta = qbar, total_se = sqrt(total), df = df,
    hazard_ratio = exp(qbar), ci_low = exp(qbar - crit * sqrt(total)), ci_high = exp(qbar + crit * sqrt(total)),
    p_value = 2 * pt(abs(qbar / sqrt(total)), df = df, lower.tail = FALSE),
    fraction_missing_information_approx = if (total > 0) ((1 + 1 / m_eff) * bvar) / total else NA_real_,
    stringsAsFactors = FALSE
  )
}
pooled <- do.call(rbind, pool_rows)

cox_result_path <- file.path(OUTPUT_ROOT, "f01_mvs_cox", "mvs_variant_cox_results.csv")
if (file.exists(cox_result_path)) {
  cc <- read.csv(cox_result_path, check.names = FALSE)
  cc <- cc[cc$variant == "corrected_per_sd" & cc$term %in% c("mvs", "proxy", "mvs_x_proxy"),
           c("cohort", "transition", "term", "hazard_ratio", "ci_low", "ci_high", "p_value", "unique_persons"), drop = FALSE]
  names(cc)[4:8] <- paste0("complete_case_", names(cc)[4:8])
  pooled <- merge(pooled, cc, by = c("cohort", "transition", "term"), all.x = TRUE)
}
write_csv(pooled, "sporadic_covariate_mi_pooled_results.csv")

write_session_log("f07_sporadic_mi", c(
  paste("Imputations:", M_IMP, "iterations:", MAXIT),
  "Scope: sporadically missing covariates among exposure-observed persons only; structural module absence not imputed."
))
message("f07 complete: ", output_dir)
