# Analysis code: metabolic vulnerability, social connection and cardiometabolic multimorbidity transitions

Analysis code for the manuscript

> Metabolic vulnerability, a harmonised social-connection proxy and cardiometabolic
> multimorbidity progression across ageing cohorts: a multistate cohort study.

The study uses three cohorts: the China Health and Retirement Longitudinal Study
(CHARLS), the English Longitudinal Study of Ageing (ELSA) and the Health and
Retirement Study (HRS).

## What is in this repository

```
R/                   R scripts for all statistical analyses (f00 to f18)
config.example.sh    example of the environment variables the scripts read
```

This repository contains analysis code only. It contains no data and no result
files: running the scripts on the input files writes the results to a folder
of your choice.

**No participant-level data are included.** CHARLS, ELSA and HRS data are
distributed by the cohort teams under data-use agreements that do not allow
redistribution.

## What is not in this repository

The scripts start from two analysis datasets (described under "Input files").
The code that builds these two datasets from the original cohort files is not
included here. Supplementary Table S26 of the manuscript lists, for every
variable, the release file and variable it comes from in each cohort and the
derivation rule.

## Software

The analyses were run in R 4.5.3 and checked again in R 4.6.1.

| R package | version | used by |
|---|---|---|
| survival | 3.8.6 | Cox models and multistate analyses |
| nnet | 7.3.20 | f02 |
| mice | 3.19.0 | f07 |
| msm | 1.8.2 | f08 |
| haven | 2.5 or later | f09 and f12, to read the official Stata files |
| data.table | 1.18 | f18, to read large CSV files |
| splines, parallel | base R | f03, f05 |

## Input files

The scripts read file locations from environment variables, so no script
contains a path. Set them in the shell (see `config.example.sh`) or in R with
`Sys.setenv()`.

| variable | needed by | what it points to |
|---|---|---|
| `CMM_OUTPUT_ROOT` | f00-f17 | folder for results (created if absent) |
| `CMM_TRANSITION_CSV` | f00-f17 | transition dataset: one row per person and follow-up interval |
| `CMM_BASELINE_CSV` | f00-f17 | baseline dataset: one row per person and wave |
| `CMM_PERSON_WAVE_DIR` | f13, f14, f16 | folder with `charls_person_wave.csv`, `elsa_person_wave.csv`, `hrs_person_wave.csv` |
| `CMM_CHARLS_2011_HEALTH_DTA` | f09 | CHARLS 2011 `health_status_and_functioning.dta` |
| `CMM_H_ELSA_DTA` | f09 | Gateway Harmonized ELSA, `h_elsa_g3.dta` |
| `CMM_H_HRS_DTA` | f09 | Gateway Harmonized HRS, `H_HRS_d.dta` |
| `CMM_CHARLS_2011_WEIGHT_DTA` | f12 | CHARLS 2011 `weight.dta` |

Optional: `FORMAL_BOOTSTRAP_B` (default 500) and `FORMAL_BOOTSTRAP_CORES`
(default 6) for f05; `FORMAL_MI_M` (default 20) for f07.

**Transition dataset.** Required columns: `cohort`, `person_id`, `tstart`,
`tstop` (years since baseline), `from_state`, `to_state`, `transition`, `event`,
`baseline_age`, `baseline_sex`, `baseline_education_level`,
`baseline_wealth_tertile`, `baseline_current_smoking`, `baseline_current_alcohol`,
`baseline_cancer`, `baseline_lung_disease`, `baseline_social_isolated_2item`,
`baseline_mvs_anthro_bp_z`. The six transitions are `no_cmd_to_single_cmd`,
`no_cmd_to_cmm`, `no_cmd_to_death`, `single_cmd_to_cmm`, `single_cmd_to_death`
and `cmm_to_death`.

**Baseline dataset.** Required columns: `cohort`, `person_id`,
`preferred_baseline_wave`, `analysis_eligible_no_cmm`, `age`, `sex`, `state`,
`living_alone`, `social_participation`, `social_isolated_2item`,
`education_level`, `wealth_tertile`, `current_smoking`, `current_alcohol`,
`baseline_cancer`, `baseline_lung_disease`, `cesd_score`, `bmi`, `waist_cm`,
`sbp`, `dbp`, `mvs_anthro_bp_z`.

## How to run

Run the scripts from the `R/` folder, because each one starts with
`source("_common.R")`.

```sh
source config.sh        # your copy of config.example.sh
cd R
Rscript f00_riskset_and_counts.R
Rscript f01_mvs_cox_primary.R
# ... and so on
```

The scripts are independent of each other. Each script writes its tables to a
sub-folder of `CMM_OUTPUT_ROOT` and a session log to `CMM_OUTPUT_ROOT/logs`.

f05 is slow (500 bootstrap samples). On Linux and macOS it uses several cores;
on Windows it runs on one core and takes a few hours. f07 and f08 also take
longer than the others.

f18 uses only public data and takes its folders as arguments:

```sh
Rscript f18_immune_crp_sensitivity.R <folder with public CSV files> <output folder>
```

The public files are those released with Gong et al., *Nature* 2025 (Sound Life
and Immunobiology of Aging cohorts, Allen Institute for Immunology).

## What each script does

| script | analysis | where reported |
|---|---|---|
| `_common.R` | shared settings and helper functions | |
| `f00_riskset_and_counts.R` | risk sets and event counts for every cohort and transition | Table 1, Table 2 |
| `f01_mvs_cox_primary.R` | primary transition-specific Cox models; alternative constructions of the score | Table 2, Fig. 2 |
| `f02_same_riskset_multinomial.R` | discrete-time multinomial models on the same risk sets | sensitivity analysis |
| `f03_two_stage_selection_ipw.R` | selection analysis and inverse-probability weighting | sensitivity analysis |
| `f04_sparse_ph_meta.R` | sparse-event sensitivity, proportional-hazards tests, random-effects pooling | sensitivity analysis |
| `f05_multistate_gcomp.R` | five-year multistate probabilities by g-computation, with bootstrap | Table 3, Fig. 4 |
| `f06_clock_reset.R` | clock-forward versus clock-reset time scales | sensitivity analysis |
| `f07_sporadic_mi.R` | multiple imputation of sporadically missing covariates | sensitivity analysis |
| `f08_msm_feasibility.R` | continuous-time multistate model for panel data | sensitivity analysis |
| `f09_medication_sensitivity.R` | antihypertensive-medication sensitivity analysis | sensitivity analysis |
| `f10_interactions_progressive_landmark.R` | progressive adjustment, age and sex interactions, two-year landmark | sensitivity analysis |
| `f11_occupancy_cmdfree.R` | state occupancy and years free of cardiometabolic disease | Fig. 3 |
| `f12_weighted_descriptive_and_followup.R` | weighted descriptive estimates; follow-up summary | sensitivity analysis |
| `f13_cmd_type_heterogeneity.R` | condition-specific models (diabetes, heart disease, stroke) | Supplementary Table S7 |
| `f14_participant_flow.R` | participant flow from all baseline respondents | Fig. 1a, Supplementary Table S1 |
| `f15_mvs_clinical_threshold.R` | clinical-threshold and BMI-free versions of the score | Supplementary Table S25 |
| `f16_included_excluded_health.R` | included versus excluded participants | Supplementary Table S3 |
| `f17_incremental_discrimination.R` | incremental discrimination of the score | Supplementary Table S21 |
| `f18_immune_crp_sensitivity.R` | external immune-ageing analysis with and without CRP | Supplementary Table S19 |

## Reproducibility notes

- Random seeds are fixed (f05 bootstrap, f07 imputation). In f05 every
  bootstrap sample sets its own seed, so the result does not depend on the
  number of cores.
- Results can differ in the last digits between R versions and operating
  systems, because of floating-point arithmetic and iterative optimisation
  (most visibly in the multinomial models of f02).
- f09 and f12 read official release files directly. See the header of each
  script for the files and variables used.

## Licence

See `LICENSE`.
