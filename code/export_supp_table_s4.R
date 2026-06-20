#!/usr/bin/env Rscript
# Export Supplementary Table S4 (CADSS timepoint contrasts) for Word paste-in.
# Source: output/tables/cadss_posthoc_timepoint_contrasts.csv

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

repo_root <- Sys.getenv("ROOT_DIR", unset = "")
if (!nzchar(repo_root)) repo_root <- getwd()
repo_root <- normalizePath(repo_root, mustWork = FALSE)

src <- file.path(repo_root, "output", "tables", "cadss_posthoc_timepoint_contrasts.csv")
out <- file.path(repo_root, "output", "revision", "tables", "supp_table_s4_manuscript_format.csv")

if (!file.exists(src)) {
  stop("Run make main first; missing ", src)
}

raw <- read_csv(src, show_col_types = FALSE)

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<.001")
  s <- sprintf("%.3f", p)
  sub("^0\\.", ".", s)
}

# Order and labels matching supplementary docx layout (ketamine/placebo blocks per subscale).
manuscript <- raw %>%
  mutate(
    Contrast = case_when(
      str_detect(contrast, fixed("(Post-bolus) - End")) ~ "Post-bolus\u2013End of infusion",
      str_detect(contrast, fixed("(Pre-infusion) - (Post-bolus)")) ~ "Pre-infusion\u2013Post-bolus",
      TRUE ~ contrast
    ),
    Session = str_to_title(session),
    Subscale = sub_scale,
    Estimate = sprintf("%.3f", estimate),
    SE = sprintf("%.3f", SE),
    df = sprintf("%.0f", df),
    `t-ratio` = sprintf("%.3f", t.ratio),
    `p-value` = vapply(p.value, format_p, character(1))
  ) %>%
  select(Contrast, Session, Subscale, Estimate, SE, df, `t-ratio`, `p-value`) %>%
  arrange(Subscale, Contrast, Session)

dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
write_csv(manuscript, out)
message("Wrote ", out, " (", nrow(manuscript), " rows)")
