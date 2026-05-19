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

# --- Supp Fig S1: weight vs CADSS (manuscript: r = -0.056, p = 0.777) ---
we <- read_csv(file.path(data_dir, "CADSS_Weight_DoseEquivalence_Data.csv"), show_col_types = FALSE)
ct <- cor.test(we$weight_kg, we$CADSS_total_post_avg)
fail("S1_pearson_r", -0.056, unname(ct$estimate), tol = 0.01)
fail("S1_pearson_p", 0.777, ct$p.value, tol = 0.01)

# Ketamine post-bolus only (session2): should NOT match manuscript if caption stats are from post_avg
ct2 <- cor.test(we$weight_kg, we$CADSS_total_post_session2)
message(
  "NOTE S1 ketamine post-bolus only: r = ",
  round(ct2$estimate, 3),
  ", p = ",
  round(ct2$p.value, 3),
  " (differs from manuscript caption)"
)

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

# --- Summary ---
n_fail <- sum(!vapply(checks, function(x) isTRUE(x$ok), logical(1)))
if (n_fail > 0) {
  message("\n", n_fail, " check(s) failed.")
  quit(status = 1)
}
message("\nAll manuscript spot-checks passed.")
