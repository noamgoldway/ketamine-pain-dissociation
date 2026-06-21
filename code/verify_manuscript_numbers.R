#!/usr/bin/env Rscript
# Spot-check committed outputs against statistics reported in Final_files docx.
# Run from repository root: Rscript code/verify_manuscript_numbers.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

repo_root <- Sys.getenv("ROOT_DIR", unset = "")
if (!nzchar(repo_root)) {
  repo_root <- getwd()
}
repo_root <- normalizePath(repo_root, mustWork = FALSE)

data_dir <- file.path(repo_root, "data")
tab_dir <- file.path(repo_root, "output", "revision", "tables")
main_tab <- file.path(repo_root, "output", "tables")

checks <- list()

fail <- function(name, expected, observed, tol = 0.001) {
  if (is.logical(expected)) {
    ok <- isTRUE(expected) == isTRUE(observed)
  } else {
    ok <- (is.na(expected) && is.na(observed)) || abs(expected - observed) <= tol
  }
  checks[[name]] <<- list(ok = ok, expected = expected, observed = observed)
  if (!ok) {
    message("FAIL ", name, ": expected ", expected, ", got ", observed)
  } else {
    message("OK   ", name)
  }
  invisible(ok)
}

# --- Supp Fig S1: weight vs ketamine post-bolus CADSS (r ~ 0.04, p ~ 0.86) ---
we <- read_csv(file.path(data_dir, "CADSS_Weight_DoseEquivalence_Data.csv"), show_col_types = FALSE)
ct <- cor.test(we$weight_kg, we$CADSS_total_post_session2)
fail("S1_n", 28, nrow(we), tol = 0)
fail("S1_pearson_r", 0.04, unname(ct$estimate), tol = 0.01)
fail("S1_pearson_p", 0.86, ct$p.value, tol = 0.02)

# --- Demographics: mean weight 61.6 ± 9.4 kg (n may be 36 enrolled vs 28 in some analyses) ---
pw <- read_csv(file.path(data_dir, "participants.csv"), show_col_types = FALSE)
fail("demog_weight_mean", 61.6, mean(pw$weight_kg, na.rm = TRUE), tol = 0.15)
fail("demog_weight_sd", 9.4, sd(pw$weight_kg, na.rm = TRUE), tol = 0.15)

# --- Steiger: all ROI p > 0.11 (main text); recomputed from data should match committed CSV ---
steiger <- read_csv(file.path(tab_dir, "steiger_roi_session_corr_diff.csv"), show_col_types = FALSE)
fail("steiger_min_p_above_0.11", 1, as.numeric(min(steiger$p) > 0.11))
fail("steiger_dlPFC_p", 0.1139, steiger$p[steiger$ROI == "Right dlPFC"])

# --- Pain model calib_temp p = 0.070 (supp table) ---
pain_m <- read_csv(file.path(main_tab, "supp_table_tempadj_pain_model.csv"), show_col_types = FALSE)
calib_row <- pain_m %>% filter(grepl("calib_temp", Effect, ignore.case = TRUE))
if (nrow(calib_row)) {
  fail("pain_calib_temp_p", 0.07, as.numeric(calib_row$p[1]), tol = 0.005)
}

# --- Rebuttal: ketamine post-bolus CADSS descriptives (mean 24.94, SD 10.77, range 5–48) ---
cad <- read_csv(file.path(data_dir, "CADSS.csv"), show_col_types = FALSE)
cad$post_bolus_ket <- cad$CADSS_Amnesiasum_2_post +
  cad$CADSS_Depersonalizationsum_2_post +
  cad$CADSS_Derealisationsum_2_post
cad_vals <- cad$post_bolus_ket[!is.na(cad$post_bolus_ket)]
fail("cadss_post_bolus_mean", 24.94, mean(cad_vals), tol = 0.05)
fail("cadss_post_bolus_sd", 10.77, sd(cad_vals), tol = 0.05)
fail("cadss_post_bolus_min", 5, min(cad_vals), tol = 0.01)
fail("cadss_post_bolus_max", 48, max(cad_vals), tol = 0.01)

# --- Main text: NPS vs pain (high intensity) Spearman rho .64 / .60 ---
nps_corr <- read_csv(file.path(main_tab, "nps_pain_corr_results.csv"), show_col_types = FALSE)
fail("nps_high_pain_placebo_rho", 0.64, nps_corr$r[nps_corr$session == "placebo"], tol = 0.02)
fail("nps_high_pain_ketamine_rho", 0.60, nps_corr$r[nps_corr$session == "ketamine"], tol = 0.02)

# --- Main / rebuttal: ROI activation + calib_temp model comparison and 3-way interaction ---
roi_cmp <- read_csv(file.path(main_tab, "model_comparison_roi_calib_temp.csv"), show_col_types = FALSE)
plus_row <- roi_cmp %>% filter(model == "plus_calib_temp")
fail("roi_calib_temp_chisq", 9.11, plus_row$Chisq[1], tol = 0.05)
fail("roi_calib_temp_chisq_p", 0.003, plus_row$`Pr(>Chisq)`[1], tol = 0.002)

roi_anova <- read_csv(file.path(main_tab, "roi_m2_nice_anova.csv"), show_col_types = FALSE)
int_row <- roi_anova %>% filter(Effect == "intensity:session:calib_temp")
fail("roi_intensity_session_calib_F", 7.16, parse_number(int_row$F[1]), tol = 0.05)
fail("roi_intensity_session_calib_p", 0.008, parse_number(int_row$p.value[1]), tol = 0.002)

# --- Rebuttal: pain ratings + calib_temp LRT (not significant) ---
pain_cmp <- read_csv(file.path(main_tab, "model_comparison_pain_calib_temp.csv"), show_col_types = FALSE)
pain_plus <- pain_cmp %>% filter(model == "plus_calib_temp")
fail("pain_calib_temp_chisq", 2.66, pain_plus$Chisq[1], tol = 0.05)
fail("pain_calib_temp_chisq_p", 0.103, pain_plus$`Pr(>Chisq)`[1], tol = 0.01)

# --- Summary ---
n_fail <- sum(!vapply(checks, function(x) isTRUE(x$ok), logical(1)))
if (n_fail > 0) {
  message("\n", n_fail, " check(s) failed.")
  quit(status = 1)
}
message("\nAll manuscript spot-checks passed.")
