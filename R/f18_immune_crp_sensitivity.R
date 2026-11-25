# f18_immune_crp_sensitivity.R
#
# Question: in the external immune-ageing analysis (Sound Life, healthy adults),
# the metabolic-inflammatory proxy contains high-sensitivity CRP. Do the
# metabolic components carry information about IL6 / TNF beyond CRP itself?
#
# This script fits, for each of the 15 prespecified immune features:
#   Model A  outcome ~ original proxy (BMI, glucose, triglycerides, 1/HDL, log hs-CRP)
#   Model B  outcome ~ CRP-free metabolic score (BMI, glucose, triglycerides, 1/HDL)
#   Model C  outcome ~ log hs-CRP alone
#   Model D  outcome ~ CRP-free metabolic score + log hs-CRP (both in one model)
# All models adjust for age, sex and CMV status, as in the original analysis.
#
# Conventions are the same as in the original analysis (Supplementary Table S19):
#   - one baseline sample per participant (smallest days since first visit);
#   - Olink values: olink.NPX_raw, median over duplicate panels on the baseline day;
#   - cell frequencies: AIFI level-3 CLR values, mean per participant;
#   - exposure, outcome and age are z-standardised, so betas are in SD units;
#   - each score is the mean of the available standardised components;
#   - Benjamini-Hochberg FDR within each block of 15 features.
#
# Input: public CSV files of Gong et al., Nature 2025 (Allen Institute).
# No participant-level CHARLS, ELSA or HRS data are used.
#
# How to run:
#   Rscript f18_immune_crp_sensitivity.R <folder with public CSV files> <output folder>

library(data.table)   # used only for fread(), to read the large CSV files quickly

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript f18_immune_crp_sensitivity.R <public_csv_dir> <output_dir>")
}
raw_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

proteins <- c("IL6", "TNF", "IL1B", "IL11", "CXCL17")
cells <- c(
  "CM CD4 T cell", "GZMB- CD27- EM CD4 T cell", "ISG+ memory CD4 T cell",
  "CD27+ effector B cell", "ISG+ memory CD8 T cell", "Core naive CD4 T cell",
  "GZMK+ CD27+ EM CD8 T cell", "CM CD8 T cell", "CD27- effector B cell",
  "Core naive CD8 T cell"
)

# z-standardise a numeric vector
z <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# 1. Olink proteins: one value per participant and protein at baseline
# ---------------------------------------------------------------------------
olink <- fread(
  file.path(raw_dir, "sound_life_all_olink.csv"),
  select = c("olink.assay", "olink.NPX_raw", "olink.panel", "subject.subjectGuid",
             "subject.biologicalSex", "subject.ageAtFirstDraw", "subject.cmv",
             "sample.daysSinceFirstVisit"),
  showProgress = FALSE
)
olink <- as.data.frame(olink)
olink <- olink[olink$olink.assay %in% proteins, ]
olink$day <- suppressWarnings(as.numeric(olink$sample.daysSinceFirstVisit))
olink$value <- suppressWarnings(as.numeric(olink$olink.NPX_raw))

# covariates: first record of each participant
meta <- olink[!duplicated(olink$subject.subjectGuid),
              c("subject.subjectGuid", "subject.biologicalSex",
                "subject.cmv", "subject.ageAtFirstDraw")]

# keep the baseline day of each participant-protein pair
olink$key <- paste(olink$subject.subjectGuid, olink$olink.assay, sep = "||")
first_day <- tapply(olink$day, olink$key, min, na.rm = TRUE)
olink$first_day <- first_day[olink$key]
olink_base <- olink[!is.na(olink$day) & olink$day == olink$first_day, ]

# median over duplicate panels
protein_long <- aggregate(value ~ subject.subjectGuid + olink.assay,
                          data = olink_base, FUN = median, na.action = na.pass)

# long to wide: one column per protein
dat <- meta
for (p in proteins) {
  one <- protein_long[protein_long$olink.assay == p, c("subject.subjectGuid", "value")]
  names(one)[2] <- paste0("protein__", p)
  dat <- merge(dat, one, by = "subject.subjectGuid", all = TRUE)
}

# ---------------------------------------------------------------------------
# 2. Immune-cell frequencies: mean CLR value per participant and cell type
# ---------------------------------------------------------------------------
cell <- fread(
  file.path(raw_dir, "sound_life_AIFI_L3_frequencies.csv"),
  select = c("subject.subjectGuid", "AIFI_L3", "AIFI_L3_clr"),
  showProgress = FALSE
)
cell <- as.data.frame(cell)
cell <- cell[cell$AIFI_L3 %in% cells, ]
cell$clr <- suppressWarnings(as.numeric(cell$AIFI_L3_clr))
cell_long <- aggregate(clr ~ subject.subjectGuid + AIFI_L3, data = cell,
                       FUN = function(v) mean(v, na.rm = TRUE), na.action = na.pass)
for (cc in cells) {
  one <- cell_long[cell_long$AIFI_L3 == cc, c("subject.subjectGuid", "clr")]
  names(one)[2] <- paste0("cell__", cc)
  dat <- merge(dat, one, by = "subject.subjectGuid", all = TRUE)
}

# ---------------------------------------------------------------------------
# 3. Clinical laboratory data: baseline record per participant, then scores
# ---------------------------------------------------------------------------
labs <- fread(
  file.path(raw_dir, "sound_life_labs_metadata.csv"),
  select = c("subject.subjectGuid", "sample.daysSinceFirstVisit", "am.bmi",
             "chem.glucose", "lip.triglycerides", "lip.cholesterol_hdl", "infl.hs_crp"),
  showProgress = FALSE
)
labs <- as.data.frame(labs)
labs$day <- suppressWarnings(as.numeric(labs$sample.daysSinceFirstVisit))
labs <- labs[order(labs$day, seq_len(nrow(labs)), na.last = TRUE), ]
labs <- labs[!duplicated(labs$subject.subjectGuid), ]

z_bmi <- z(labs$am.bmi)
z_glucose <- z(labs$chem.glucose)
z_tg <- z(labs$lip.triglycerides)
z_inv_hdl <- z(1 / suppressWarnings(as.numeric(labs$lip.cholesterol_hdl)))
log_crp <- log(suppressWarnings(as.numeric(labs$infl.hs_crp)))
z_log_crp <- z(log_crp)

scores <- data.frame(
  subject.subjectGuid = labs$subject.subjectGuid,
  # original proxy: five components, including CRP
  metabolic_inflammatory_score = rowMeans(cbind(z_bmi, z_glucose, z_tg, z_inv_hdl, z_log_crp),
                                          na.rm = TRUE),
  # CRP-free score: the four metabolic components only
  metabolic_score_crp_free = rowMeans(cbind(z_bmi, z_glucose, z_tg, z_inv_hdl), na.rm = TRUE),
  log_hs_crp = log_crp,
  n_metabolic_components = rowSums(!is.na(cbind(z_bmi, z_glucose, z_tg, z_inv_hdl))),
  stringsAsFactors = FALSE
)
dat <- merge(dat, scores, by = "subject.subjectGuid", all = TRUE)

# ---------------------------------------------------------------------------
# 4. Covariates and standardised exposures
# ---------------------------------------------------------------------------
dat$sex_male <- as.numeric(dat$subject.biologicalSex == "Male")
dat$cmv_positive <- NA
dat$cmv_positive[dat$subject.cmv == "Positive"] <- 1
dat$cmv_positive[dat$subject.cmv == "Negative"] <- 0
dat$age_z <- z(dat$subject.ageAtFirstDraw)

dat$x_original <- z(dat$metabolic_inflammatory_score)
dat$x_crp_free <- z(dat$metabolic_score_crp_free)
dat$x_crp <- z(dat$log_hs_crp)

outcomes <- c(paste0("protein__", proteins), paste0("cell__", cells))

# ---------------------------------------------------------------------------
# 5. Models
# ---------------------------------------------------------------------------
# Fit one linear model and return the row of the requested term.
fit_one <- function(data, outcome, exposure_terms, report_term, model_label) {
  data$y <- z(data[[outcome]])
  f <- as.formula(paste("y ~", paste(c(exposure_terms, "sex_male", "cmv_positive", "age_z"),
                                    collapse = " + ")))
  fit <- lm(f, data = data)
  co <- summary(fit)$coefficients
  ci <- confint(fit)
  data.frame(
    model = model_label,
    term = report_term,
    outcome = outcome,
    n = nobs(fit),
    beta = co[report_term, "Estimate"],
    se = co[report_term, "Std. Error"],
    ci_low = ci[report_term, 1],
    ci_high = ci[report_term, 2],
    p = co[report_term, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

# Run one block of 15 outcomes and add the FDR within the block.
run_block <- function(data, exposure_terms, report_term, model_label, sample_label) {
  rows <- list()
  for (o in outcomes) {
    rows[[length(rows) + 1]] <- fit_one(data, o, exposure_terms, report_term, model_label)
  }
  out <- do.call(rbind, rows)
  out$fdr <- p.adjust(out$p, method = "BH")
  out$sample <- sample_label
  out
}

# Sample 1: all participants available to each model (as in the original analysis)
res1 <- rbind(
  run_block(dat, "x_original", "x_original", "A. original proxy (with CRP)", "available cases"),
  run_block(dat, "x_crp_free", "x_crp_free", "B. CRP-free metabolic score", "available cases"),
  run_block(dat, "x_crp", "x_crp", "C. log hs-CRP alone", "available cases"),
  run_block(dat, c("x_crp_free", "x_crp"), "x_crp_free",
            "D. joint model: CRP-free metabolic score", "available cases"),
  run_block(dat, c("x_crp_free", "x_crp"), "x_crp",
            "D. joint model: log hs-CRP", "available cases")
)

# Sample 2: the same participants in every model
# (hs-CRP and all four metabolic components observed)
common <- dat[!is.na(dat$log_hs_crp) & !is.na(dat$n_metabolic_components) &
                dat$n_metabolic_components == 4, ]
res2 <- rbind(
  run_block(common, "x_original", "x_original", "A. original proxy (with CRP)", "common sample"),
  run_block(common, "x_crp_free", "x_crp_free", "B. CRP-free metabolic score", "common sample"),
  run_block(common, "x_crp", "x_crp", "C. log hs-CRP alone", "common sample"),
  run_block(common, c("x_crp_free", "x_crp"), "x_crp_free",
            "D. joint model: CRP-free metabolic score", "common sample"),
  run_block(common, c("x_crp_free", "x_crp"), "x_crp",
            "D. joint model: log hs-CRP", "common sample")
)

results <- rbind(res1, res2)
write.csv(results, file.path(out_dir, "f18_immune_crp_sensitivity_results.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Sample description and checks
# ---------------------------------------------------------------------------
has_il6 <- !is.na(dat$protein__IL6)
summary_rows <- data.frame(
  item = c(
    "participants in merged table",
    "participants with IL6 measured",
    "with IL6 and original proxy",
    "with IL6 and CRP-free metabolic score",
    "with IL6 and hs-CRP",
    "with IL6, hs-CRP and all four metabolic components",
    "correlation: CRP-free metabolic score vs log hs-CRP (participants with IL6)",
    "correlation: original proxy vs CRP-free metabolic score (participants with IL6)"
  ),
  value = c(
    nrow(dat),
    sum(has_il6),
    sum(has_il6 & !is.na(dat$x_original)),
    sum(has_il6 & !is.na(dat$x_crp_free)),
    sum(has_il6 & !is.na(dat$x_crp)),
    sum(has_il6 & !is.na(dat$log_hs_crp) & !is.na(dat$n_metabolic_components) &
          dat$n_metabolic_components == 4),
    cor(dat$metabolic_score_crp_free[has_il6], dat$log_hs_crp[has_il6], use = "complete.obs"),
    cor(dat$metabolic_inflammatory_score[has_il6], dat$metabolic_score_crp_free[has_il6],
        use = "complete.obs")
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_rows, file.path(out_dir, "f18_sample_summary.csv"), row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "f18_session.txt"))

cat("\nSample summary:\n")
print(summary_rows, row.names = FALSE)
cat("\nIL6 and TNF, all models:\n")
print(results[results$outcome %in% c("protein__IL6", "protein__TNF"),
              c("sample", "model", "outcome", "n", "beta", "ci_low", "ci_high", "p", "fdr")],
      row.names = FALSE, digits = 3)
