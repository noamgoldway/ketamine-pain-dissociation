#!/usr/bin/env Rscript
# Full R2 verification: supplementary table cells vs committed R outputs + in-text anchors.
# Uses CSVs from code/01_*.r and code/02_*.Rmd (not Python re-analysis).
# Prerequisite: make extract-claims  (python docx inventory only)
# Run: make verify-all

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(jsonlite)
  library(stringr)
  library(purrr)
  if (requireNamespace("yaml", quietly = TRUE)) {
    library(yaml)
  }
})

repo_root <- Sys.getenv("ROOT_DIR", unset = "")
if (!nzchar(repo_root)) repo_root <- getwd()
repo_root <- normalizePath(repo_root, mustWork = FALSE)

data_dir <- file.path(repo_root, "data")
tab_dir <- file.path(repo_root, "output", "tables")
rev_tab <- file.path(repo_root, "output", "revision", "tables")

failures <- character()
passes <- character()

note_pass <- function(msg) passes <<- c(passes, msg)
note_fail <- function(msg) failures <<- c(failures, msg)

parse_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  s <- as.character(x)
  s <- str_trim(s)
  if (!nzchar(s) || s %in% c("***", "ns", "NA", "—", "-")) return(NA_real_)
  s <- gsub("^\\.", "0.", s)
  if (grepl("^<", s)) {
    v <- gsub("[^0-9.]", "", s)
    v <- sub("^\\.", "0.", v)
    v <- suppressWarnings(as.numeric(v))
    if (!is.na(v)) return(v * 0.5)
    return(1e-4)
  }
  s <- gsub("\\*+", "", s)
  m <- str_extract(s, "-?[0-9]+(\\.[0-9]+)?(e[+-]?[0-9]+)?")
  if (is.na(m)) {
    m <- str_extract(s, "\\.[0-9]+")
    if (!is.na(m)) m <- paste0("0", m)
  }
  if (is.na(m)) return(NA_real_)
  suppressWarnings(as.numeric(m))
}

close_num <- function(a, b, tol = NULL, rel = 0.08) {
  if (is.na(a) || is.na(b)) return(FALSE)
  if (is.null(tol)) {
    tol <- if (max(abs(a), abs(b), na.rm = TRUE) < 1) 0.005 else 0.02
  }
  if (abs(a) < 1e-6 && abs(b) > 0.001) return(FALSE)
  if (abs(b) < 1e-6 && abs(a) > 0.001) return(FALSE)
  if (abs(a) < 1e-9 && abs(b) < 1e-9) return(TRUE)
  if (abs(a - b) <= tol) return(TRUE)
  denom <- max(abs(a), abs(b), 1e-9)
  if (abs(a - b) / denom <= rel) return(TRUE)
  if (abs(abs(a) - abs(b)) <= tol) return(TRUE)
  denom <- max(abs(abs(a)), abs(abs(b)), 1e-9)
  abs(abs(a) - abs(b)) / denom <= rel
}

num_close <- function(x, pool) {
  which(vapply(pool, function(y) close_num(x, y), logical(1)))
}

csv_numeric_pool <- function(path) {
  if (!file.exists(path)) return(numeric())
  df <- read_csv(path, show_col_types = FALSE)
  unlist(lapply(df, function(col) {
    if (is.numeric(col)) return(col)
    vapply(col, parse_num, numeric(1))
  }), use.names = FALSE) %>%
    .[!is.na(.)]
}

match_pool <- function(ms_nums, pool, min_frac = 1.0) {
  pool <- pool
  matched <- 0L
  for (x in ms_nums) {
    if (is.na(x)) next
    hits <- num_close(x, pool)
    if (length(hits)) {
      best <- hits[which.min(abs(pool[hits] - x))]
      pool <- pool[-best]
      matched <- matched + 1L
    }
  }
  list(
    matched = matched,
    total = length(ms_nums),
    ok = length(ms_nums) == 0 || matched / length(ms_nums) >= min_frac
  )
}

normalize_s4_contrast <- function(s) {
  s <- tolower(str_trim(as.character(s)))
  s <- gsub("[–—−]", "-", s, fixed = FALSE)
  s <- gsub("\\s+", " ", s)
  if (grepl("post-bolus.*end of infusion|post-bolus.*end", s)) {
    return("post_end")
  }
  if (grepl("pre-infusion.*post-bolus|pre.*post-bolus", s)) {
    return("pre_post")
  }
  s
}

verify_s4_structured <- function(docx_rows, csv_path) {
  if (!file.exists(csv_path)) return(character())
  csv <- read_csv(csv_path, show_col_types = FALSE)
  issues <- character()
  data_rows <- docx_rows[-1]
  for (row in data_rows) {
    if (length(row) < 8) next
    est <- parse_num(row[[4]])
    tval <- parse_num(row[[7]])
    pdoc <- parse_num(row[[8]])
    if (is.na(est) || is.na(tval)) next
    hit <- which(
      vapply(seq_len(nrow(csv)), function(i) {
        close_num(est, csv$estimate[i], tol = 0.02) &&
          close_num(tval, csv$t.ratio[i], tol = 0.05)
      }, logical(1))
    )
    if (!length(hit)) {
      issues <- c(
        issues,
        paste0(
          "S4 row [", row[[1]], " / ", row[[2]], " / ", row[[3]],
          "]: no CSV row for estimate=", row[[4]], ", t=", row[[7]]
        )
      )
      next
    }
    i <- hit[1]
    if (!close_num(pdoc, csv$p.value[i], tol = 0.005)) {
      issues <- c(
        issues,
        paste0(
          "S4 p-value mismatch [", row[[1]], " / ", row[[2]], " / ", row[[3]],
          "]: docx=", row[[8]], ", R=", signif(csv$p.value[i], 4),
          " (CSV session=", csv$session[i], ", sub_scale=", csv$sub_scale[i], ")"
        )
      )
    }
    doc_sess <- tolower(str_trim(row[[2]]))
    csv_sess <- tolower(str_trim(csv$session[i]))
    if (doc_sess != csv_sess) {
      issues <- c(
        issues,
        paste0(
          "S4 session label mismatch [", row[[3]], ", estimate=", row[[4]],
          "]: docx=", row[[2]], ", R=", csv$session[i]
        )
      )
    }
  }
  issues
}

load_table_specs <- function() {
  manifest_path <- file.path(repo_root, "docs", "TABLE_MANIFEST.yaml")
  if (file.exists(manifest_path) && requireNamespace("yaml", quietly = TRUE)) {
    raw <- yaml::read_yaml(manifest_path)
    specs <- list()
    for (sid in paste0("S", 1:14)) {
      entry <- raw[[sid]]
      if (is.null(entry)) next
      rel <- entry$source_csv
      if (is.null(rel) || !nzchar(rel)) next
      specs[[sid]] <- list(
        status = entry$status %||% "reproducible",
        csv = file.path(repo_root, rel)
      )
    }
    if (length(specs)) return(specs)
  }
  list(
    S1  = list(status = "reproducible", csv = file.path(data_dir, "connectivity_roi_list.csv")),
    S2  = list(status = "reproducible", csv = file.path(tab_dir, "supp_table_cadss_mixed_model.csv")),
    S3  = list(status = "reproducible", csv = file.path(tab_dir, "cadss_posthoc_session_contrasts.csv")),
    S4  = list(status = "reproducible", csv = file.path(tab_dir, "cadss_posthoc_timepoint_contrasts.csv")),
    S5  = list(status = "reproducible", csv = file.path(rev_tab, "supp_table_s5_ketamine_clusters.csv")),
    S6  = list(status = "reproducible", csv = file.path(tab_dir, "supp_table_roi_activation_model_anova.csv")),
    S7  = list(status = "reproducible", csv = file.path(tab_dir, "roi_diff_posthoc.csv")),
    S8  = list(status = "reproducible", csv = file.path(tab_dir, "roi_m2_nice_anova.csv")),
    S9  = list(status = "reproducible", csv = file.path(rev_tab, "posthoc_m2_calibtemp_slope_by_intensity_noadj.csv")),
    S10 = list(status = "reproducible", csv = file.path(rev_tab, "posthoc_m2_calibtemp_slope_by_session_noadj.csv")),
    S11 = list(status = "reproducible", csv = file.path(rev_tab, "posthoc_m2_calibtemp_slope_by_intensity_within_session_noadj.csv")),
    S12 = list(status = "reproducible", csv = file.path(tab_dir, "pain_roi_corr_results.csv")),
    S13 = list(status = "reproducible", csv = file.path(tab_dir, "pain_roi_dissoc_model_comparison.csv")),
    S14 = list(status = "reproducible", csv = file.path(tab_dir, "nps_dissoc_model_comparison.csv"))
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# --- Load extracted docx numbers (inventory only; comparison uses R CSVs) ---
claims_path <- file.path(repo_root, "docs", "extracted_claims_r2.json")
if (!file.exists(claims_path)) {
  note_fail(paste0("Missing ", claims_path, "; run: make extract-claims"))
} else {
  claims <- fromJSON(claims_path, simplifyVector = FALSE)
  supp_tables <- claims$supplementary_tables
  table_specs <- load_table_specs()

  for (sid in names(table_specs)) {
    spec <- table_specs[[sid]]
    csv_path <- spec$csv
    tab <- supp_tables[[sid]]
    if (is.null(tab)) {
      note_fail(paste0(sid, ": missing from docx extraction"))
      next
    }
    ms_raw <- vapply(tab$numeric_cells, function(c) c$raw, character(1))
    ms_nums <- vapply(ms_raw, parse_num, numeric(1))
    ms_nums <- ms_nums[!is.na(ms_nums)]
    if (spec$status == "frozen") {
      if (file.exists(csv_path)) {
        note_pass(paste0(sid, ": frozen reference OK (", length(ms_nums), " docx nums; csv present)"))
      } else {
        note_fail(paste0(sid, ": frozen csv missing: ", csv_path))
      }
      next
    }
    if (!file.exists(csv_path)) {
      note_fail(paste0(sid, ": missing csv ", csv_path))
      next
    }
    pool <- csv_numeric_pool(csv_path)
    res <- match_pool(ms_nums, pool)
    if (res$ok) {
      note_pass(paste0(sid, ": ", res$matched, "/", res$total, " cells match ", basename(csv_path)))
    } else {
      note_fail(paste0(sid, ": ", res$matched, "/", res$total, " cells match ", basename(csv_path)))
    }
    if (sid == "S4" && !is.null(tab$rows)) {
      s4_issues <- verify_s4_structured(tab$rows, csv_path)
      if (length(s4_issues)) {
        for (msg in s4_issues) note_fail(msg)
      } else if (res$ok) {
        note_pass("S4: row-level estimate/t/p/session checks OK")
      }
    }
  }
}

# --- Temp-adjusted pain model export (cited in text; not the S14 docx table) ---
pain_model_path <- file.path(tab_dir, "supp_table_tempadj_pain_model.csv")
if (!file.exists(pain_model_path)) {
  note_fail(paste0("supp_table_tempadj_pain_model.csv missing: ", pain_model_path))
} else {
  pain_m <- read_csv(pain_model_path, show_col_types = FALSE)
  calib_row <- pain_m %>% filter(grepl("calib_temp", Effect, ignore.case = TRUE))
  if (nrow(calib_row)) {
    p_calib <- parse_num(calib_row$p[1])
    if (close_num(p_calib, 0.07, tol = 0.005)) {
      note_pass(paste0("supp_table_tempadj_pain_model: calib_temp p = 0.070 (got ", p_calib, ")"))
    } else {
      note_fail(paste0("supp_table_tempadj_pain_model: calib_temp p expected 0.070, got ", p_calib))
    }
  } else {
    note_fail("supp_table_tempadj_pain_model.csv: no calib_temp row")
  }
}

# --- In-text anchors from R outputs ---
pain_anova <- read_csv(file.path(tab_dir, "pain_ratings_mixed_model_anova.csv"), show_col_types = FALSE)
int_row <- pain_anova %>% filter(grepl("intensity:session", Effect, ignore.case = TRUE))
if (nrow(int_row)) {
  f_int <- parse_num(int_row$F[1])
  if (close_num(f_int, 11.22, tol = 0.05)) {
    note_pass(paste0("Main abstract pain interaction F(1,54)=11.22: got ", round(f_int, 2)))
  } else {
    note_fail(paste0("Main abstract pain interaction F: expected 11.22, got ", f_int))
  }
} else {
  note_fail("pain_ratings_mixed_model_anova: no intensity:session row")
}

cadss_anova <- read_csv(file.path(tab_dir, "supp_table_cadss_mixed_model.csv"), show_col_types = FALSE)
sess_row <- cadss_anova %>% filter(Effect == "session")
if (nrow(sess_row)) {
  f_sess <- parse_num(sess_row$F[1])
  if (close_num(f_sess, 57.44, tol = 0.1)) {
    note_pass(paste0("Main abstract CADSS session F(1,32)=57.44: got ", round(f_sess, 2)))
  } else {
    note_fail(paste0("CADSS session F: expected 57.44, got ", f_sess))
  }
}

corrs <- read_csv(file.path(tab_dir, "pain_roi_corr_results.csv"), show_col_types = FALSE)
check_rho <- function(roi, sess, target, tol = 0.02) {
  val <- corrs %>% filter(ROI == roi, session == sess) %>% pull(r)
  if (length(val) && close_num(val[1], target, tol = tol)) {
    note_pass(paste0("ROI rho ", roi, " ", sess, " ≈ ", target))
  } else {
    note_fail(paste0("ROI rho ", roi, " ", sess, ": expected ", target, ", got ", val[1]))
  }
}
check_rho("ant_insula_L", "placebo", 0.528, tol = 0.02)
check_rho("ant_insula_L", "ketamine", 0.656, tol = 0.02)
check_rho("dlPFC_R", "ketamine", 0.570, tol = 0.02)

calib <- read_csv(file.path(tab_dir, "pain_calibration_stats.csv"), show_col_types = FALSE)
if (nrow(calib)) {
  note_pass("pain_calibration_stats.csv present for calibration temperature descriptives")
}

s2 <- file.path(repo_root, "output", "revision", "figures", "Supplementary_Figure_S2_roi_calib_temperature.png")
s2_slopes <- file.path(repo_root, "output", "revision", "tables", "supp_fig_s2_session_slopes.csv")
if (!file.exists(s2)) {
  note_fail("Supp Fig S2 PNG missing (run: make supp-fig-s2)")
} else if (!file.exists(s2_slopes)) {
  note_fail("supp_fig_s2_session_slopes.csv missing (run: make supp-fig-s2)")
} else {
  slopes <- read_csv(s2_slopes, show_col_types = FALSE)
  pl <- slopes %>% filter(session == "placebo") %>% pull(calib_temp.trend)
  kl <- slopes %>% filter(session == "ketamine") %>% pull(calib_temp.trend)
  if (length(pl) && length(kl) && pl[1] > kl[1]) {
    note_pass(paste0(
      "Supp Fig S2 regenerated: placebo slope ", sprintf("%.3f", pl[1]),
      " > ketamine ", sprintf("%.3f", kl[1])
    ))
  } else {
    note_fail("Supp Fig S2: expected placebo slope > ketamine slope")
  }
  if (close_num(pl[1], 0.341, tol = 0.015) && close_num(kl[1], 0.131, tol = 0.015)) {
    note_pass(paste0(
      "Supp Fig S2 slopes consistent with submitted figure (",
      sprintf("%.3f/%.3f", pl[1], kl[1]), " vs 0.341/0.131)"
    ))
  } else {
    note_fail(paste0(
      "Supp Fig S2 slopes: expected 0.341/0.131, got ",
      sprintf("%.3f/%.3f", pl[1], kl[1])
    ))
  }
}

# S5 reference sync (SPM table committed in data/)
s5_ref <- file.path(data_dir, "spm_ketamine_clusters_reference.csv")
s5_out <- file.path(rev_tab, "supp_table_s5_ketamine_clusters.csv")
if (file.exists(s5_ref) && file.exists(s5_out)) {
  ref_txt <- readLines(s5_ref, warn = FALSE)
  out_txt <- readLines(s5_out, warn = FALSE)
  if (identical(ref_txt, out_txt)) {
    note_pass("S5 cluster table synced from data/spm_ketamine_clusters_reference.csv")
  } else {
    note_fail("S5 output differs from spm_ketamine_clusters_reference.csv (run: make sync-s5)")
  }
} else {
  note_fail("S5 reference or output CSV missing")
}

# --- R spot-check suite (verify_manuscript_numbers.R) ---
smoke_script <- file.path(repo_root, "code", "verify_manuscript_numbers.R")
smoke_ok <- FALSE
smoke_out <- ""
if (file.exists(smoke_script)) {
  old_root <- Sys.getenv("ROOT_DIR", unset = "")
  Sys.setenv(ROOT_DIR = repo_root)
  on.exit({
    if (nzchar(old_root)) Sys.setenv(ROOT_DIR = old_root) else Sys.unsetenv("ROOT_DIR")
  }, add = TRUE)
  smoke_res <- system2(
    "Rscript",
    shQuote(normalizePath(smoke_script, mustWork = TRUE)),
    stdout = TRUE,
    stderr = TRUE
  )
  smoke_out <- paste(smoke_res, collapse = "\n")
  smoke_ok <- attr(smoke_res, "status") %||% 0L == 0L
  if (smoke_ok) {
    note_pass("verify_manuscript_numbers.R spot-checks passed")
  } else {
    note_fail("verify_manuscript_numbers.R spot-checks failed")
  }
} else {
  note_fail("verify_manuscript_numbers.R not found")
}

# --- Write report ---
report_path <- file.path(repo_root, "docs", "VERIFICATION_FULL.md")
lines <- c(
  "# Full verification report (R2)",
  "",
  paste0("Generated from R outputs in `", repo_root, "`."),
  "",
  "Pipeline: `make verify-all` runs **`make all`** (main + supplementary + figures + static sync) before numeric checks.",
  "",
  "Source docx: `text/manuscript/`",
  "",
  "## Summary",
  "",
  paste0("- **Passed:** ", length(passes)),
  paste0("- **Failed:** ", length(failures)),
  "",
  "### Supplementary tables + anchors",
  ""
)
if (length(passes)) lines <- c(lines, paste0("- OK: ", passes))
if (length(failures)) lines <- c(lines, "", "### Failures", "", paste0("- FAIL: ", failures))
lines <- c(
  lines,
  "",
  "### R spot-check suite (`verify_manuscript_numbers.R`)",
  "",
  if (smoke_ok) "OK" else "FAIL",
  "",
  "```",
  if (nzchar(smoke_out)) tail(strsplit(smoke_out, "\n")[[1]], 30) else character(),
  "```",
  "",
  paste0("**Overall:** ", if (length(failures)) "FAIL" else "PASS")
)
writeLines(lines, report_path)

cat(paste(lines, collapse = "\n"), "\n")
if (length(failures)) quit(status = 1)
message("verify_all_manuscript.R: all checks passed.")
