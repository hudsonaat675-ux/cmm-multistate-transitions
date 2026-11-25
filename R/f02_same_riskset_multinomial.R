#!/usr/bin/env Rscript
# f02: discrete-time competing-transition analysis on exactly the Model 3 risk
# sets (pooled multinomial per starting state, person-cluster sandwich variance).
# Formal-run change vs pre-analysis: tightened convergence (reltol 1e-10).

options(stringsAsFactors = FALSE, width = 160)
for (pkg in c("nnet", "sandwich")) if (!requireNamespace(pkg, quietly = TRUE)) stop("Package ", pkg, " required")
source("_common.R")

output_dir <- make_output_dir("f02_multinomial")
write_csv <- make_write_csv(output_dir)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions

covars <- MODEL3_COVARS

model_formula <- outcome ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex +
  factor(baseline_education_level) + factor(baseline_wealth_tertile) +
  baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + factor(interval_index)

outcome_spec <- list(
  no_cmd = c("stay", "no_cmd_to_single_cmd", "no_cmd_to_cmm", "no_cmd_to_death"),
  single_cmd = c("stay", "single_cmd_to_cmm", "single_cmd_to_death")
)

prepare_riskset <- function(cohort_name, from_state_name) {
  x <- transitions[
    transitions$cohort == cohort_name & transitions$from_state == from_state_name &
      is.finite(transitions$tstart) & is.finite(transitions$tstop) &
      transitions$tstop > transitions$tstart,
    , drop = FALSE
  ]
  x$exposure <- as.numeric(x$corrected_per_sd)
  x$exposure_x_proxy <- x$exposure * as.numeric(x$baseline_social_isolated_2item)
  x$outcome <- ifelse(x$event == 1, x$transition, "stay")
  allowed <- outcome_spec[[from_state_name]]
  if (any(!x$outcome %in% allowed)) stop("Unexpected outcome in ", cohort_name, " / ", from_state_name)
  x$outcome <- factor(x$outcome, levels = allowed)
  required <- c("person_id", "outcome", "interval_index", "exposure", "exposure_x_proxy", covars)
  x <- x[complete.cases(x[, required, drop = FALSE]), , drop = FALSE]
  droplevels(x)
}

multinom_cluster_vcov <- function(fit, data) {
  bread <- stats::vcov(fit)
  xmat <- stats::model.matrix(fit)
  prob <- stats::fitted(fit)
  if (is.null(dim(prob))) stop("Multinomial fitted probabilities did not form a matrix")
  outcome_levels <- levels(data$outcome)
  nonbaseline <- outcome_levels[-1]
  prob <- prob[, nonbaseline, drop = FALSE]
  y <- stats::model.matrix(~ outcome - 1, data = data)
  colnames(y) <- sub("^outcome", "", colnames(y))
  y <- y[, nonbaseline, drop = FALSE]
  score_blocks <- lapply(nonbaseline, function(cl) xmat * (y[, cl] - prob[, cl]))
  scores_i <- do.call(cbind, score_blocks)
  score_names <- unlist(lapply(nonbaseline, function(cl) paste(cl, colnames(xmat), sep = ":")), use.names = FALSE)
  colnames(scores_i) <- score_names
  vnames <- colnames(bread)
  if (!setequal(vnames, score_names)) stop("Coefficient/score names do not agree")
  scores_i <- scores_i[, vnames, drop = FALSE]
  cluster_scores <- rowsum(scores_i, group = data$person_id, reorder = FALSE)
  n <- nrow(data); p <- ncol(scores_i); g <- nrow(cluster_scores)
  correction <- if (g > 1 && n > p) (g / (g - 1)) * ((n - 1) / (n - p)) else 1
  meat <- crossprod(cluster_scores) * correction
  robust <- bread %*% meat %*% bread
  dimnames(robust) <- dimnames(bread)
  list(vcov = robust, max_abs_total_score = max(abs(colSums(scores_i))),
       clusters = g, observations = n, parameters = p, correction = correction)
}

extract_multinom <- function(fit, robust, cohort_name, from_state_name, data) {
  bmat <- coef(fit)
  if (is.null(dim(bmat))) bmat <- matrix(bmat, nrow = 1, dimnames = list(levels(data$outcome)[-1], names(bmat)))
  term_map <- c(exposure = "mvs", baseline_social_isolated_2item = "proxy", exposure_x_proxy = "mvs_x_proxy")
  rows <- list()
  for (outcome_name in rownames(bmat)) {
    for (term_name in names(term_map)) {
      coefficient_name <- paste(outcome_name, term_name, sep = ":")
      beta <- unname(bmat[outcome_name, term_name])
      se <- sqrt(robust$vcov[coefficient_name, coefficient_name])
      rows[[length(rows) + 1]] <- data.frame(
        cohort = cohort_name, from_state = from_state_name, transition = outcome_name,
        term = unname(term_map[term_name]), beta = beta, cluster_robust_se = se,
        relative_risk_ratio = exp(beta), ci_low = exp(beta - 1.96 * se), ci_high = exp(beta + 1.96 * se),
        p_value = 2 * pnorm(abs(beta / se), lower.tail = FALSE),
        risk_intervals = nrow(data), events = sum(data$outcome == outcome_name),
        unique_persons = length(unique(data$person_id)), stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

multinom_rows <- list()
diagnostic_rows <- list()
riskset_rows <- list()

for (cohort_name in COHORTS) {
  for (from_state_name in names(outcome_spec)) {
    x <- prepare_riskset(cohort_name, from_state_name)
    level_counts <- table(x$outcome)
    for (outcome_name in names(level_counts)) {
      riskset_rows[[length(riskset_rows) + 1]] <- data.frame(
        cohort = cohort_name, from_state = from_state_name, outcome = outcome_name,
        risk_intervals = nrow(x), events_or_stays = unname(level_counts[outcome_name]),
        unique_persons = length(unique(x$person_id)), stringsAsFactors = FALSE
      )
    }
    z <- fit_capture(function() nnet::multinom(
      model_formula, data = x, Hess = TRUE, trace = FALSE,
      maxit = 2000, MaxNWts = 20000, reltol = 1e-10, abstol = 1e-12
    ))
    robust <- NULL
    robust_error <- NULL
    if (!is.null(z$fit)) {
      robust <- tryCatch(multinom_cluster_vcov(z$fit, x), error = function(e) { robust_error <<- conditionMessage(e); NULL })
    }
    status <- if (is.null(z$fit)) "model_error" else if (is.null(robust)) "robust_vcov_error" else "fit"
    diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
      cohort = cohort_name, from_state = from_state_name, model = "multinomial", status = status,
      convergence_code = if (is.null(z$fit) || !length(z$fit$convergence)) NA_integer_ else z$fit$convergence,
      risk_intervals = nrow(x), unique_persons = length(unique(x$person_id)),
      parameters = if (is.null(robust)) NA_integer_ else robust$parameters,
      max_abs_total_score = if (is.null(robust)) NA_real_ else robust$max_abs_total_score,
      warning_or_error = paste(unique(c(z$warning, robust_error)), collapse = " | "),
      stringsAsFactors = FALSE
    )
    if (!is.null(robust)) {
      multinom_rows[[length(multinom_rows) + 1]] <- extract_multinom(z$fit, robust, cohort_name, from_state_name, x)
    }
  }
}

multinom_results <- do.call(rbind, multinom_rows)
write_csv(multinom_results, "same_riskset_multinomial_results.csv")
write_csv(do.call(rbind, diagnostic_rows), "same_riskset_model_diagnostics.csv")
write_csv(do.call(rbind, riskset_rows), "same_riskset_outcome_counts.csv")

# Side-by-side with the formal corrected Cox once f01 outputs exist.
cox_path <- file.path(OUTPUT_ROOT, "f01_mvs_cox", "mvs_variant_cox_results.csv")
if (file.exists(cox_path)) {
  cox <- read.csv(cox_path, check.names = FALSE)
  cox <- cox[cox$variant == "corrected_per_sd" & cox$term == "mvs",
             c("cohort", "transition", "hazard_ratio", "ci_low", "ci_high", "p_value", "events"), drop = FALSE]
  names(cox)[3:7] <- paste0("cox_", names(cox)[3:7])
  multi <- multinom_results[multinom_results$term == "mvs",
                            c("cohort", "transition", "relative_risk_ratio", "ci_low", "ci_high", "p_value", "events"), drop = FALSE]
  names(multi)[3:7] <- paste0("multinomial_", names(multi)[3:7])
  comparison <- merge(cox, multi, by = c("cohort", "transition"), all = TRUE)
  write_csv(comparison, "mvs_cox_multinomial_comparison.csv")
}

write_session_log("f02_same_riskset_multinomial", c(
  "Convergence: nnet::multinom with reltol 1e-10, abstol 1e-12, maxit 2000."
))
message("f02 complete: ", output_dir)
