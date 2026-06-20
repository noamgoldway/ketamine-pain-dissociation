#!/usr/bin/env Rscript
# Export the 12-row S4 subset shown in the supplementary docx (2 contrasts × 6 blocks).
# Use this to paste-replace Table S4 in Word.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

repo_root <- Sys.getenv("ROOT_DIR", unset = "")
if (!nzchar(repo_root)) repo_root <- getwd()
repo_root <- normalizePath(repo_root, mustWork = FALSE)

src <- file.path(repo_root, "output", "revision", "tables", "supp_table_s4_manuscript_format.csv")
out <- file.path(repo_root, "output", "revision", "tables", "supp_table_s4_word_paste.csv")

format_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<.001")
  s <- sprintf("%.3f", p)
  sub("^0\\.", ".", s)
}

raw <- read_csv(src, show_col_types = FALSE)

word_order <- tibble(
  Contrast = c(
    rep("Post-bolus\u2013End of infusion", 6),
    rep("Pre-infusion\u2013Post-bolus", 6)
  ),
  Session = rep(rep(c("Placebo", "Ketamine"), each = 3), 2),
  Subscale = rep(c("Amnesia", "Depersonalization", "Derealisation"), 4)
)

out_tbl <- word_order %>%
  left_join(raw, by = c("Contrast", "Session", "Subscale")) %>%
  select(Contrast, Session, Subscale, Estimate, SE, df, `t-ratio`, `p-value`)

write_csv(out_tbl, out)
message("Wrote ", out, " (", nrow(out_tbl), " rows — matches Word S4 layout)")
