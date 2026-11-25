#!/usr/bin/env Rscript
# f16: included vs excluded participants within the analysis-eligible no-CMM
# pool, on health and socioeconomic characteristics.
# Excluded participants are split into stage A (exposure module
# unavailable) and stage B (exposure observed, covariate or interval missing).
# Standardised mean differences (SMD) are relative to the included group.

options(stringsAsFactors = FALSE, width = 180)
source("_common.R")

output_dir <- make_output_dir("f16_included_excluded_health")
write_csv <- make_write_csv(output_dir)

PERSON_WAVE_DIR <- env_path("CMM_PERSON_WAVE_DIR")
BASE_WAVE <- c(CHARLS = 1, ELSA = 4, HRS = 10)

inputs <- load_frozen_inputs()
transitions <- inputs$transitions
eligible <- inputs$eligible

model3_persons <- function(cohort_name) {
  unique(unlist(lapply(seq_len(nrow(TRANSITION_MAP)), function(i) {
    prepare_model3_riskset(transitions, cohort_name, TRANSITION_MAP$transition[i], TRANSITION_MAP$from_state[i])$person_id
  })))
}

smd <- function(a, b, binary) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  if (binary) {
    p1 <- mean(a); p0 <- mean(b); den <- sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
  } else {
    p1 <- mean(a); p0 <- mean(b); den <- sqrt((stats::var(a) + stats::var(b)) / 2)
  }
  if (!is.finite(den) || den == 0) return(NA_real_)
  (p1 - p0) / den
}

rows <- list()
for (cohort_name in COHORTS) {
  x <- eligible[eligible$cohort == cohort_name, , drop = FALSE]
  p <- read.csv(file.path(PERSON_WAVE_DIR, paste0(tolower(cohort_name), "_person_wave.csv")), check.names = FALSE)
  p <- p[p$wave == BASE_WAVE[[cohort_name]], c("person_id", "diabetes", "heart", "stroke")]
  x <- merge(x, p, by = "person_id", all.x = TRUE, sort = FALSE)
  cc <- model3_persons(cohort_name)
  exposure_ok <- is.finite(x$corrected_per_sd) & !is.na(x$social_isolated_2item)
  x$group <- ifelse(x$person_id %in% cc, "included",
             ifelse(!exposure_ok, "excluded_stage_a_exposure_unavailable", "excluded_stage_b_covariate_or_interval"))
  x$excluded_any <- x$group != "included"
  x$male <- as.integer(x$sex == 1)
  x$single_cmd <- as.integer(x$state == "single_cmd")
  x$died_during_followup <- as.integer(x$has_future_death_record %in% 1)
  for (lv in sort(unique(na.omit(x$education_level)))) x[[paste0("education_", lv)]] <- as.integer(x$education_level == lv)
  for (lv in 1:3) x[[paste0("wealth_tertile_", lv)]] <- as.integer(x$wealth_tertile == lv)

  vars <- c(age = FALSE, male = TRUE, single_cmd = TRUE, diabetes = TRUE, heart = TRUE, stroke = TRUE,
            baseline_cancer = TRUE, baseline_lung_disease = TRUE, cesd_score = FALSE,
            current_smoking = TRUE, current_alcohol = TRUE, living_alone = TRUE,
            died_during_followup = TRUE)
  for (nm in grep("^(education|wealth_tertile)_[0-9]+$", names(x), value = TRUE)) vars[nm] <- TRUE

  inc <- x[x$group == "included", , drop = FALSE]
  comparisons <- list(all_excluded = x[x$excluded_any, , drop = FALSE],
                      stage_a = x[x$group == "excluded_stage_a_exposure_unavailable", , drop = FALSE],
                      stage_b = x[x$group == "excluded_stage_b_covariate_or_interval", , drop = FALSE])
  for (v in names(vars)) {
    binary <- vars[[v]]
    summ <- function(z) if (binary) 100 * mean(z, na.rm = TRUE) else mean(z, na.rm = TRUE)
    sdv <- function(z) if (binary) NA_real_ else stats::sd(z, na.rm = TRUE)
    rows[[length(rows) + 1]] <- data.frame(
      cohort = cohort_name, variable = v, type = if (binary) "percent" else "mean",
      included_n = nrow(inc), included_value = summ(inc[[v]]), included_sd = sdv(inc[[v]]),
      included_missing = sum(is.na(inc[[v]])),
      excluded_n = nrow(comparisons$all_excluded), excluded_value = summ(comparisons$all_excluded[[v]]),
      excluded_sd = sdv(comparisons$all_excluded[[v]]), excluded_missing = sum(is.na(comparisons$all_excluded[[v]])),
      smd_excluded_vs_included = smd(comparisons$all_excluded[[v]], inc[[v]], binary),
      stage_a_n = nrow(comparisons$stage_a), stage_a_value = summ(comparisons$stage_a[[v]]),
      smd_stage_a_vs_included = smd(comparisons$stage_a[[v]], inc[[v]], binary),
      stage_b_n = nrow(comparisons$stage_b), stage_b_value = summ(comparisons$stage_b[[v]]),
      smd_stage_b_vs_included = smd(comparisons$stage_b[[v]], inc[[v]], binary))
  }
}
out <- do.call(rbind, rows)
write_csv(out, "included_vs_excluded_health_ses.csv")
write_csv(out[abs(out$smd_excluded_vs_included) >= 0.1 & !is.na(out$smd_excluded_vs_included), ],
          "included_vs_excluded_smd_ge_0.1.csv")
write_session_log("f16_included_excluded_health", c(
  "Groups within the analysis-eligible no-CMM pool; included = Model-3 complete case.",
  "CES-D is on cohort-specific scales (CHARLS 0-30; ELSA/HRS 0-8): compare within cohort only.",
  "died_during_followup = any later death record (from baseline-candidate file)."))
message("f16 complete: ", output_dir)
