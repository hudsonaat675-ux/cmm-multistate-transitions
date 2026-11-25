#!/bin/bash
# Example settings. Copy this file to config.sh, edit the paths, then run
#   source config.sh
# before running the R scripts. See README.md for what each file must contain.

# --- needed by f00 to f17 ------------------------------------------------------
export CMM_OUTPUT_ROOT="/path/to/outputs"                      # results are written here
export CMM_TRANSITION_CSV="/path/to/transitions_selected_survival_ready.csv"
export CMM_BASELINE_CSV="/path/to/baseline_candidates.csv"

# --- needed by f13, f14 and f16 --------------------------------------------------
export CMM_PERSON_WAVE_DIR="/path/to/folder_with_person_wave_csv_files"

# --- needed only by f09 (antihypertensive-medication sensitivity analysis) -------
export CMM_CHARLS_2011_HEALTH_DTA="/path/to/CHARLS/2011/health_status_and_functioning.dta"
export CMM_H_ELSA_DTA="/path/to/h_elsa_g3.dta"
export CMM_H_HRS_DTA="/path/to/H_HRS_d.dta"

# --- needed only by f12 (CHARLS weighted descriptive analysis) -------------------
export CMM_CHARLS_2011_WEIGHT_DTA="/path/to/CHARLS/2011/weight.dta"

# --- optional ---------------------------------------------------------------------
export FORMAL_BOOTSTRAP_B=500        # f05: number of bootstrap samples
export FORMAL_BOOTSTRAP_CORES=6      # f05: number of parallel workers
export FORMAL_MI_M=20                # f07: number of imputations
