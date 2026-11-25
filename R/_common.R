# Shared settings and helper functions. Every f*.R script starts with
# source("_common.R"), so the scripts must be run from this folder.
#
# What this file does:
#   - reads the two analysis datasets (transition intervals and baseline records);
#   - applies the half-month rule: a death recorded in the same month as the
#     preceding interview is given a duration of half a month;
#   - defines the analysis-eligible population (baseline without CMM);
#   - computes the exposure used in all models: the metabolic vulnerability
#     score standardised within each cohort's eligible population (per 1 SD);
#   - builds the risk set of one cohort and one transition for the primary
#     model (Model 3);
#   - writes a session log for each script.

ANALYSIS_ID <- "REV1-FORMAL-20260828"
GLOBAL_SEED <- 20260828
HALF_MONTH_YEARS <- 0.5 / 12

# File paths are read from environment variables, so no script contains a
# path. README.md lists the variables; config.example.sh sets all of them.
env_path <- function(var, default = NULL, must_exist = TRUE) {
  value <- Sys.getenv(var, unset = NA_character_)
  if (is.na(value) || !nzchar(value)) {
    if (is.null(default)) {
      stop(sprintf("Set %s (see README.md and config.example.sh).", var), call. = FALSE)
    }
    value <- default
  }
  if (must_exist && !file.exists(value)) {
    stop(sprintf("%s points at a path that does not exist: %s", var, value), call. = FALSE)
  }
  value
}

PROJECT_ROOT <- Sys.getenv("CMM_PROJECT_ROOT", unset = getwd())
OUTPUT_ROOT <- Sys.getenv("CMM_OUTPUT_ROOT", unset = file.path(PROJECT_ROOT, "outputs"))
TRANSITION_PATH <- env_path("CMM_TRANSITION_CSV")
BASELINE_PATH <- env_path("CMM_BASELINE_CSV")

COHORTS <- c("CHARLS", "ELSA", "HRS")

TRANSITION_MAP <- data.frame(
  transition = c("no_cmd_to_single_cmd", "no_cmd_to_cmm", "no_cmd_to_death",
                 "single_cmd_to_cmm", "single_cmd_to_death", "cmm_to_death"),
  from_state = c("no_cmd", "no_cmd", "no_cmd", "single_cmd", "single_cmd", "cmm"),
  stringsAsFactors = FALSE
)

MODEL3_COVARS <- c(
  "baseline_age", "baseline_sex", "baseline_education_level",
  "baseline_wealth_tertile", "baseline_current_smoking",
  "baseline_current_alcohol", "baseline_cancer", "baseline_lung_disease",
  "baseline_social_isolated_2item"
)

zscore <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  if (sum(ok) > 1) {
    s <- stats::sd(x[ok])
    if (is.finite(s) && s > 0) out[ok] <- (x[ok] - mean(x[ok])) / s
  }
  out
}

make_output_dir <- function(step_name) {
  out <- file.path(OUTPUT_ROOT, step_name)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

make_write_csv <- function(output_dir) {
  function(x, name) write.csv(x, file.path(output_dir, name), row.names = FALSE, na = "")
}

fit_capture <- function(expr_fun) {
  warning_text <- character()
  error_text <- NULL
  fit <- tryCatch(
    withCallingHandlers(
      expr_fun(),
      warning = function(w) {
        warning_text <<- c(warning_text, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      error_text <<- conditionMessage(e)
      NULL
    }
  )
  list(fit = fit, warning = paste(unique(c(warning_text, error_text)), collapse = " | "))
}

# Load transitions and baseline; apply the half-month rule BEFORE any
# validity filtering so models and descriptive counts share one convention.
load_frozen_inputs <- function() {
  transitions <- read.csv(TRANSITION_PATH, check.names = FALSE)
  baseline <- read.csv(BASELINE_PATH, check.names = FALSE)

  nonpositive <- is.finite(transitions$tstart) & is.finite(transitions$tstop) &
    transitions$tstop <= transitions$tstart
  death_row <- transitions$event == 1 & grepl("_death$", transitions$transition)
  fix_rows <- nonpositive & death_row
  drop_rows <- nonpositive & !death_row
  transitions$half_month_rule_applied <- as.integer(fix_rows)
  transitions$tstop[fix_rows] <- transitions$tstart[fix_rows] + HALF_MONTH_YEARS
  transitions$duration_years[fix_rows] <- transitions$tstop[fix_rows] - transitions$tstart[fix_rows]
  if (any(drop_rows)) transitions <- transitions[!drop_rows, , drop = FALSE]

  eligible <- baseline[baseline$preferred_baseline_wave == 1 & baseline$analysis_eligible_no_cmm == 1, , drop = FALSE]
  eligible <- eligible[!duplicated(paste(eligible$cohort, eligible$person_id, sep = "__")), , drop = FALSE]
  eligible$key <- paste(eligible$cohort, eligible$person_id, sep = "__")

  scores <- do.call(rbind, lapply(COHORTS, function(cohort_name) {
    x <- eligible[eligible$cohort == cohort_name, , drop = FALSE]
    data.frame(cohort = cohort_name, person_id = x$person_id,
               corrected_per_sd = zscore(as.numeric(x$mvs_anthro_bp_z)), stringsAsFactors = FALSE)
  }))
  score_key <- paste(scores$cohort, scores$person_id, sep = "__")
  transition_key <- paste(transitions$cohort, transitions$person_id, sep = "__")
  transitions$corrected_per_sd <- scores$corrected_per_sd[match(transition_key, score_key)]
  eligible$corrected_per_sd <- scores$corrected_per_sd[match(eligible$key, score_key)]

  list(
    transitions = transitions,
    baseline = baseline,
    eligible = eligible,
    scores = scores,
    n_half_month_fixed = sum(fix_rows),
    n_nonpositive_dropped = sum(drop_rows),
    input_md5 = c(transitions = unname(tools::md5sum(TRANSITION_PATH)),
                  baseline = unname(tools::md5sum(BASELINE_PATH)))
  )
}

# Risk set of the primary model (Model 3) for one cohort and one transition,
# using the per-SD exposure, after the half-month rule.
prepare_model3_riskset <- function(transitions, cohort_name, transition_name, from_state_name,
                                   exposure_var = "corrected_per_sd") {
  x <- transitions[
    transitions$cohort == cohort_name & transitions$from_state == from_state_name &
      is.finite(transitions$tstart) & is.finite(transitions$tstop) &
      transitions$tstop > transitions$tstart,
    , drop = FALSE
  ]
  x$event_target <- as.integer(x$event == 1 & x$transition == transition_name)
  x$exposure <- as.numeric(x[[exposure_var]])
  x$exposure_x_proxy <- x$exposure * as.numeric(x$baseline_social_isolated_2item)
  required <- c("person_id", "tstart", "tstop", "event_target", "exposure", "exposure_x_proxy", MODEL3_COVARS)
  x[complete.cases(x[, required, drop = FALSE]), , drop = FALSE]
}

write_session_log <- function(step_name, extra_lines = character()) {
  log_dir <- file.path(OUTPUT_ROOT, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_text <- capture.output({
    cat("Analysis ID:", ANALYSIS_ID, "\n")
    cat("Run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
    cat("Half-month rule: same-time death intervals set to tstart + 0.5/12 years.\n")
    cat("Transition input md5:", unname(tools::md5sum(TRANSITION_PATH)), "\n")
    cat("Baseline input md5:", unname(tools::md5sum(BASELINE_PATH)), "\n")
    for (line in extra_lines) cat(line, "\n")
    print(sessionInfo())
  })
  writeLines(log_text, file.path(log_dir, paste0(step_name, "_session.txt")))
}
