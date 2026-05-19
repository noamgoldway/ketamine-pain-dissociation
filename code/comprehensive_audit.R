#!/usr/bin/env Rscript
# Extract numeric claims from Final_files docx and compare to repo outputs.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

repo_root <- Sys.getenv("ROOT_DIR", unset = "")
if (!nzchar(repo_root)) repo_root <- normalizePath("..", mustWork = FALSE)
repo_root <- normalizePath(repo_root, mustWork = TRUE)

final_dir <- file.path(dirname(dirname(repo_root)), "revision", "Final_files")
if (!dir.exists(final_dir)) {
  final_dir <- file.path(repo_root, "..", "..", "revision", "Final_files")
}
final_dir <- normalizePath(final_dir, mustWork = FALSE)

docx_files <- c(
  Main = "Main_revision_clean.docx",
  Supplementary = "Supplementary_Information_revision_clean.docx",
  Rebuttal = "rebuttal_letter.docx"
)

extract_text <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  tmp <- tempfile(fileext = ".xml")
  cmd <- sprintf(
    "unzip -p %s word/document.xml 2>/dev/null | sed 's/<[^>]*>/ /g' | tr -s ' '",
    shQuote(path)
  )
  system(cmd, intern = TRUE)
}

patterns <- list(
  s1_r = "Pearson r\\s*=\\s*0\\.04",
  s1_p = "p\\s*=\\s*0\\.86",
  steiger_dlPFC = "0\\.1139",
  calib_pain_p = "p\\s*=\\s*0\\.070",
  nps_64 = "0\\.64",
  nps_60 = "0\\.60",
  cadss_mean = "24\\.94",
  chi2_911 = "9\\.11",
  F_716 = "7\\.16"
)

cat("=== Docx presence (revision/Final_files) ===\n")
for (nm in names(docx_files)) {
  p <- file.path(final_dir, docx_files[[nm]])
  cat(nm, ":", if (file.exists(p)) "found" else "MISSING", "\n")
}

cat("\n=== Repo spot-checks (independent of verify_manuscript_numbers.R) ===\n")

we <- read_csv(file.path(repo_root, "data", "CADSS_Weight_DoseEquivalence_Data.csv"), show_col_types = FALSE)
ct <- cor.test(we$weight_kg, we$CADSS_total_post_session2)
cat(sprintf("S1: n=%d r=%.3f p=%.3f\n", nrow(we), ct$estimate, ct$p.value))

pw <- read_csv(file.path(repo_root, "data", "participants.csv"), show_col_types = FALSE)
cat(sprintf("Weight: %.1f +/- %.1f (n=%d)\n", mean(pw$weight_kg), sd(pw$weight_kg), nrow(pw)))

cadss <- read_csv(file.path(repo_root, "data", "CADSS.csv"), show_col_types = FALSE)
post2 <- rowSums(cadss[, grep("_2_post$", names(cadss))], na.rm = TRUE)
cat(sprintf("CADSS post-bolus ketamine: mean=%.2f sd=%.2f min=%d max=%d (n=%d)\n",
            mean(post2), sd(post2), min(post2), max(post2), length(post2)))

nps <- read_csv(file.path(repo_root, "output", "tables", "nps_pain_corr_results.csv"), show_col_types = FALSE)
cat("NPS high pain correlations:\n")
print(nps %>% filter(grepl("high", contrast, ignore.case = TRUE)) %>% select(contrast, rho, p))

steiger <- read_csv(file.path(repo_root, "output", "revision", "tables", "steiger_roi_session_corr_diff.csv"), show_col_types = FALSE)
cat(sprintf("Steiger min p=%.4f; dlPFC p=%.4f\n", min(steiger$p), steiger$p[steiger$ROI == "Right dlPFC"]))

pain_m <- read_csv(file.path(repo_root, "output", "tables", "supp_table_tempadj_pain_model.csv"), show_col_types = FALSE)
calib <- pain_m %>% filter(grepl("calib_temp", Effect, ignore.case = TRUE))
cat(sprintf("Pain calib_temp p=%.3f\n", calib$p[1]))

roi_lrt <- read_csv(file.path(repo_root, "output", "tables", "model_comparison_roi_calib_temp.csv"), show_col_types = FALSE)
cat("ROI calib_temp LRT:\n")
print(roi_lrt)

anova_m2 <- read_csv(file.path(repo_root, "output", "tables", "roi_m2_nice_anova.csv"), show_col_types = FALSE)
row <- anova_m2 %>% filter(grepl("intensity.*session.*calib", `Effect`, ignore.case = TRUE))
if (nrow(row)) cat(sprintf("ROI 3-way: F=%.2f p=%.3f\n", row$`F`[1], row$`p`[1]))

cat("\n=== Table file counts ===\n")
cat("output/tables:", length(list.files(file.path(repo_root, "output", "tables"), pattern = "\\.csv$")), "csv\n")
cat("output/figures:", length(list.files(file.path(repo_root, "output", "figures"))), "files\n")
cat("output/revision/tables:", length(list.files(file.path(repo_root, "output", "revision", "tables"), pattern = "\\.csv$")), "csv\n")
cat("output/revision/figures:", length(list.files(file.path(repo_root, "output", "revision", "figures"))), "files\n")
