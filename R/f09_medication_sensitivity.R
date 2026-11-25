# f09: antihypertensive-medication sensitivity analysis.
#
# Baseline blood-pressure medication is taken from the official release files:
#   CHARLS 2011  health_status_and_functioning.dta
#                da007_1_  doctor-diagnosed hypertension (1 = yes)
#                da011s1   taking Chinese traditional medicine for hypertension
#                da011s2   taking Western modern medicine for hypertension
#                (da011 is asked only of participants with diagnosed hypertension,
#                 so participants without a diagnosis are coded 0)
#   ELSA wave 4  Gateway Harmonized ELSA (h_elsa_g3.dta): r4rxhibp, r4hibpe
#   HRS wave 10  Gateway Harmonized HRS (H_HRS_d.dta):    r10rxhibp, r10hibpe
#
# Definition: bp_med = 1 if antihypertensive medication use was reported at the
# baseline wave; 0 otherwise (no report, including participants not asked).
#
# Two sensitivity models on the same risk sets as the primary model (Model 3):
#   A. medication-covariate model: Model 3 + bp_med
#   B. treatment-corrected score: +15 mmHg systolic / +10 mmHg diastolic for
#      treated participants (Tobin et al.), components re-standardised, composite
#      rebuilt (at least 3 of 4 components) and re-scaled.
#
# Run from the R/ folder:  Rscript f09_medication_sensitivity.R

options(stringsAsFactors = FALSE, width = 170)
library(survival)
library(haven)      # read_dta(), to read the official Stata files
source("_common.R")

output_dir <- make_output_dir("f09_medication")
write_csv <- make_write_csv(output_dir)

# Extra input files for this script (see README.md, "Input files")
CHARLS_HEALTH_DTA <- env_path("CMM_CHARLS_2011_HEALTH_DTA")
H_ELSA_DTA <- env_path("CMM_H_ELSA_DTA")
H_HRS_DTA <- env_path("CMM_H_HRS_DTA")

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
eligible <- inputs$eligible

# Person identifiers are compared as whole numbers written without leading zeros,
# which is how they appear after read.csv() of the analysis dataset.
id_text <- function(x) {
  sprintf("%.0f", suppressWarnings(as.numeric(as.character(x))))
}

# A multiple-choice item is "selected" when it holds a non-missing, non-zero code.
is_selected <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  !is.na(x) & x != 0
}

# ---- 1. Build the baseline medication flag per cohort -----------------------

# ---- 1a. CHARLS 2011 ----
ch <- read_dta(CHARLS_HEALTH_DTA)
ch <- as.data.frame(ch)
names(ch) <- tolower(names(ch))

# The 2011 identifier has 11 characters (9-character household + 2-character
# person). Later releases use 12 characters, with "0" inserted after the
# household part. The analysis dataset uses the 12-character form.
id_2011 <- as.character(ch$id)
id_release <- paste0(substr(id_2011, 1, 9), "0", substr(id_2011, 10, 11))
ch$person_id <- id_text(id_release)

ch_eligible_ids <- id_text(eligible$person_id[eligible$cohort == "CHARLS"])
ch <- ch[ch$person_id %in% ch_eligible_ids, ]

ch_diagnosed <- suppressWarnings(as.numeric(ch$da007_1_)) %in% 1
ch_med_chinese <- is_selected(ch$da011s1)
ch_med_western <- is_selected(ch$da011s2)
ch_med_reported <- ch_med_chinese | ch_med_western
ch_item_missing <- ch_diagnosed & is.na(ch$da011s1) & is.na(ch$da011s2)

med_charls <- data.frame(cohort = "CHARLS", person_id = ch$person_id,
                         bp_med = as.integer(ch_med_reported))
quality_charls <- data.frame(
  cohort = "CHARLS",
  eligible_linked_n = nrow(ch),
  diagnosed_hypertension_n = sum(ch_diagnosed),
  bp_med_reported_n = sum(ch_med_reported),
  diagnosed_item_missing_n = sum(ch_item_missing),
  source = "CHARLS 2011 health_status_and_functioning.dta: da011s1, da011s2"
)

# ---- 1b. ELSA wave 4 ----
el <- read_dta(H_ELSA_DTA, col_select = c("idauniq", "r4rxhibp", "r4hibpe"))
el <- as.data.frame(el)
el$person_id <- id_text(el$idauniq)
el_eligible_ids <- id_text(eligible$person_id[eligible$cohort == "ELSA"])
el <- el[el$person_id %in% el_eligible_ids, ]
el_med_reported <- suppressWarnings(as.numeric(el$r4rxhibp)) %in% 1
el_diagnosed <- suppressWarnings(as.numeric(el$r4hibpe)) %in% 1

med_elsa <- data.frame(cohort = "ELSA", person_id = el$person_id,
                       bp_med = as.integer(el_med_reported))
quality_elsa <- data.frame(
  cohort = "ELSA",
  eligible_linked_n = nrow(el),
  diagnosed_hypertension_n = sum(el_diagnosed),
  bp_med_reported_n = sum(el_med_reported),
  diagnosed_item_missing_n = sum(el_diagnosed & is.na(el$r4rxhibp)),
  source = "Harmonized ELSA: r4rxhibp"
)

# ---- 1c. HRS wave 10 ----
hr <- read_dta(H_HRS_DTA, col_select = c("hhidpn", "r10rxhibp", "r10hibpe"))
hr <- as.data.frame(hr)
hr$person_id <- id_text(hr$hhidpn)
hr_eligible_ids <- id_text(eligible$person_id[eligible$cohort == "HRS"])
hr <- hr[hr$person_id %in% hr_eligible_ids, ]
hr_med_reported <- suppressWarnings(as.numeric(hr$r10rxhibp)) %in% 1
hr_diagnosed <- suppressWarnings(as.numeric(hr$r10hibpe)) %in% 1

med_hrs <- data.frame(cohort = "HRS", person_id = hr$person_id,
                      bp_med = as.integer(hr_med_reported))
quality_hrs <- data.frame(
  cohort = "HRS",
  eligible_linked_n = nrow(hr),
  diagnosed_hypertension_n = sum(hr_diagnosed),
  bp_med_reported_n = sum(hr_med_reported),
  diagnosed_item_missing_n = sum(hr_diagnosed & is.na(hr$r10rxhibp)),
  source = "Harmonized HRS: r10rxhibp"
)

# ---- 1d. Combine the three cohorts ----
med <- rbind(med_charls, med_elsa, med_hrs)
quality <- rbind(quality_charls, quality_elsa, quality_hrs)
quality$eligible_n <- NA
for (k in seq_len(nrow(quality))) {
  quality$eligible_n[k] <- sum(eligible$cohort == quality$cohort[k])
}
quality$linked_percent_of_eligible <- 100 * quality$eligible_linked_n / quality$eligible_n
quality$bp_med_percent_of_linked <- 100 * quality$bp_med_reported_n / quality$eligible_linked_n
write_csv(quality, "bp_medication_source_quality.csv")

# Participants who are eligible but absent from the medication source: bp_med = 0.
med_key <- paste(med$cohort, med$person_id, sep = "__")
eligible_key <- paste(eligible$cohort, id_text(eligible$person_id), sep = "__")
eligible$bp_med <- med$bp_med[match(eligible_key, med_key)]
eligible$bp_med[is.na(eligible$bp_med)] <- 0L

# ---- 2. Treatment-corrected MVS composite (Tobin +15/+10) -------------------

components <- c("bmi", "waist_cm", "sbp", "dbp")
tobin_parts <- list()
for (cohort_name in COHORTS) {
  x <- eligible[eligible$cohort == cohort_name, , drop = FALSE]
  x$sbp_adj <- as.numeric(x$sbp) + ifelse(x$bp_med == 1, 15, 0)
  x$dbp_adj <- as.numeric(x$dbp) + ifelse(x$bp_med == 1, 10, 0)
  zmat <- cbind(
    bmi = zscore(as.numeric(x$bmi)),
    waist_cm = zscore(as.numeric(x$waist_cm)),
    sbp = zscore(x$sbp_adj),
    dbp = zscore(x$dbp_adj)
  )
  ncomp <- rowSums(is.finite(zmat))
  raw <- rowMeans(zmat, na.rm = TRUE)
  raw[ncomp < 3] <- NA_real_
  tobin_parts[[cohort_name]] <- data.frame(
    cohort = cohort_name, person_id = x$person_id,
    tobin_per_sd = zscore(raw), bp_med = x$bp_med, stringsAsFactors = FALSE
  )
}
tobin <- do.call(rbind, tobin_parts)
tobin_key <- paste(tobin$cohort, tobin$person_id, sep = "__")
transition_key <- paste(transitions$cohort, transitions$person_id, sep = "__")
transitions$tobin_per_sd <- tobin$tobin_per_sd[match(transition_key, tobin_key)]
transitions$bp_med <- tobin$bp_med[match(transition_key, tobin_key)]
transitions$bp_med[is.na(transitions$bp_med)] <- 0L

# ---- 3. Sensitivity Cox fits on exact Model 3 risk sets ---------------------

formula_medcov <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + bp_med + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

formula_plain <- Surv(tstart, tstop, event_target) ~ exposure + baseline_social_isolated_2item +
  exposure_x_proxy + baseline_age + baseline_sex + factor(baseline_education_level) +
  factor(baseline_wealth_tertile) + baseline_current_smoking + baseline_current_alcohol +
  baseline_cancer + baseline_lung_disease + cluster(person_id)

result_rows <- list()
diagnostic_rows <- list()
for (cohort_name in COHORTS) {
  for (i in seq_len(nrow(TRANSITION_MAP))) {
    tr <- TRANSITION_MAP$transition[i]
    fs <- TRANSITION_MAP$from_state[i]
    for (model_name in c("medication_covariate", "tobin_corrected_score")) {
      exposure_var <- if (model_name == "tobin_corrected_score") "tobin_per_sd" else "corrected_per_sd"
      x <- prepare_model3_riskset(transitions, cohort_name, tr, fs, exposure_var = exposure_var)
      fml <- if (model_name == "medication_covariate") formula_medcov else formula_plain
      z <- fit_capture(function() coxph(fml, data = x, ties = "efron", x = TRUE))
      diagnostic_rows[[length(diagnostic_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = tr, model = model_name,
        risk_intervals = nrow(x), events = sum(x$event_target),
        unique_persons = length(unique(x$person_id)),
        status = if (is.null(z$fit)) "error" else "fit", warning_or_error = z$warning,
        stringsAsFactors = FALSE
      )
      if (is.null(z$fit)) next
      keep_map <- c(exposure = "mvs", baseline_social_isolated_2item = "proxy", exposure_x_proxy = "mvs_x_proxy")
      if (model_name == "medication_covariate") keep_map <- c(keep_map, bp_med = "bp_med")
      keep <- intersect(names(keep_map), names(coef(z$fit)))
      b <- coef(z$fit)[keep]
      se <- sqrt(diag(vcov(z$fit)))[keep]
      result_rows[[length(result_rows) + 1]] <- data.frame(
        cohort = cohort_name, transition = tr, model = model_name,
        term = unname(keep_map[keep]), beta = unname(b), robust_se = unname(se),
        hazard_ratio = exp(unname(b)), ci_low = exp(unname(b) - 1.96 * unname(se)),
        ci_high = exp(unname(b) + 1.96 * unname(se)),
        p_value = 2 * pnorm(abs(unname(b / se)), lower.tail = FALSE),
        risk_intervals = nrow(x), events = sum(x$event_target), stringsAsFactors = FALSE
      )
    }
  }
}
results <- do.call(rbind, result_rows)
write_csv(results, "medication_sensitivity_cox_results.csv")
write_csv(do.call(rbind, diagnostic_rows), "medication_sensitivity_cox_diagnostics.csv")

# Compact comparison against the primary corrected model for the MVS term.
cox_path <- file.path(OUTPUT_ROOT, "f01_mvs_cox", "mvs_variant_cox_results.csv")
if (file.exists(cox_path)) {
  cox <- read.csv(cox_path, check.names = FALSE)
  cox <- cox[cox$variant == "corrected_per_sd" & cox$term == "mvs",
             c("cohort", "transition", "hazard_ratio", "ci_low", "ci_high"), drop = FALSE]
  names(cox)[3:5] <- paste0("primary_", names(cox)[3:5])
  sens <- results[results$term == "mvs", c("cohort", "transition", "model", "hazard_ratio", "ci_low", "ci_high"), drop = FALSE]
  comparison <- merge(sens, cox, by = c("cohort", "transition"), all.x = TRUE)
  write_csv(comparison, "medication_vs_primary_mvs_comparison.csv")
}

write_session_log("f09_medication_sensitivity", c(
  "bp_med definition: reported antihypertensive medication use at baseline wave; no report coded 0 (structural zeros included).",
  "Tobin correction: +15 mmHg SBP / +10 mmHg DBP for treated before component z-scoring.",
  "Sources: CHARLS 2011 health_status_and_functioning.dta (da011s1, da011s2); Harmonized ELSA r4rxhibp; Harmonized HRS r10rxhibp."
))
message("f09 complete: ", output_dir)
