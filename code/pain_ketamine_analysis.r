#!/usr/bin/env Rscript
# =============================================================================
#  Pain × Ketamine Study – Publication‑Ready Analysis Pipeline
#  Author: Noam Goldway
#  Date: Jul 24 2025
#  Description: End‑to‑end, reproducible script for demographics, behavioral,
#               CADSS, ROI/fMRI, connectivity, and correlation analyses.
#               Outputs: clean tables (CSV/RTF) and figures (PNG/PDF/SVG).
# =============================================================================

# ------------------------------ 0. Housekeeping ------------------------------
# Clear environment (optional when sourcing)
rm(list = ls())

# Reproducibility
set.seed(20240101)

# Project root handling --------------------------------------------------------
# Allow the user to set the root directory when running the script.
# Priority order:
#   1. Command line argument  --root_dir /path/to/project
#   2. Environment variable   ROOT_DIR
#   3. getOption("root_dir") if set in .Rprofile
#   4. here::here() fallback

# Small parser for --key value style args (no extra deps like optparse)
parse_cli_arg <- function(flag = "--root_dir") {
  args <- commandArgs(trailingOnly = TRUE)
  hit  <- which(args == flag)
  if (length(hit) == 1 && length(args) >= hit + 1) return(normalizePath(args[hit + 1], mustWork = FALSE))
  NULL
}

root_dir <- parse_cli_arg()
if (is.null(root_dir)) root_dir <- Sys.getenv("ROOT_DIR", unset = NA_character_)
if (is.na(root_dir) || root_dir == "") root_dir <- getOption("root_dir", default = NA)
if (is.na(root_dir) || root_dir == "") {
  if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
  root_dir <- here::here()
}

root_dir <- normalizePath(root_dir, mustWork = FALSE)

# ---- Paths ---------------------------------------------------------------
# Option A: use root_dir logic above (default)
# Option B: hard-set everything to an OSF folder. Set OSF_ROOT to "" to disable.
OSF_ROOT <- ""

if (nzchar(OSF_ROOT)) {
  PATHS <- list(
    data    = file.path(OSF_ROOT, "data"),
    tables  = file.path(OSF_ROOT, "output", "tables"),
    figures = file.path(OSF_ROOT, "output", "figures")
  )
} else {
  PATHS <- list(
    data    = file.path(root_dir, "data"),
    tables  = file.path(root_dir, "output", "tables"),
    figures = file.path(root_dir, "output", "figures")
  )
}

# Create output dirs if missing
invisible(lapply(PATHS[c("tables", "figures")], dir.create, recursive = TRUE, showWarnings = FALSE))

# ------------------------------ 1. Libraries --------------------------------- ---------------------------------
# Only load what is used. Install on the fly if missing.
needed_pkgs <- c(
  "tidyverse", "readxl", "broom", "broom.mixed", "lme4", "emmeans",
  "afex", "Hmisc", "patchwork", "ggtext", "janitor", "corrplot",
  "ComplexHeatmap", "circlize", "grid", "gridExtra"
)

installed <- rownames(installed.packages())
if (any(!needed_pkgs %in% installed)) {
  install.packages(setdiff(needed_pkgs, installed))
}

# Load
lapply(needed_pkgs, library, character.only = TRUE)

# afex default options for Kenward-Roger etc.
afex::set_sum_contrasts()  # sum-to-zero contrasts
options(afex.type = 3)     # Type III SS

# ------------------------------ 2. Helpers -----------------------------------
# 2.1 Simple summary function for demographics
get_demog_summary <- function(df, age_var = age, gender_var = gender) {
  df <- dplyr::rename(df, age = {{ age_var }}, gender = {{ gender_var }})
  tibble(
    N        = nrow(df),
    mean_age = mean(df$age, na.rm = TRUE),
    sd_age   = sd(df$age, na.rm = TRUE)
  ) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
    mutate(gender_1 = sum(df$gender == 1, na.rm = TRUE),
           gender_2 = sum(df$gender == 2, na.rm = TRUE))
}

# 2.2 Convenience wrappers -----------------------------------------------------
# Error bar function for ggplot
mean_se <- function(x, na.rm = TRUE) {
  m  <- mean(x, na.rm = na.rm)
  se <- sd(x, na.rm = na.rm) / sqrt(sum(!is.na(x)))
  c(y = m, ymin = m - se, ymax = m + se)
}

save_plot <- function(p, filename, width = 7, height = 5, dpi = 300, formats = c("png", "pdf")) {
  purrr::walk(formats, ~ ggsave(filename = file.path(PATHS$figures, paste0(filename, ".", .x)),
                                plot = p, width = width, height = height, dpi = dpi, device = .x))
}

save_table <- function(df, filename) {
  readr::write_csv(df, file.path(PATHS$tables, paste0(filename, ".csv")))
}

# Significance stars helper
p_stars <- function(p) dplyr::case_when(
  p < 0.001 ~ "***",
  p < 0.01  ~ "**",
  p < 0.05  ~ "*",
  TRUE      ~ ""
)

# ------------------------------ 3. Data In -----------------------------------
# Subject ID vectors -----------------------------------------------------------

both_sessions <- c("0435", "0614", "0642", "0839", "0863", "0991", "1008", "1137", "1175", "1405",
                   "1842", "2060", "2494", "2711", "3440", "3571", "3911", "5005", "5071", "5407",
                   "7826", "7856", "7932", "8295", "8298", "8446", "8754", "9018", "9235", "9364",
                   "9501", "9616", "9653")

# Demographics -----------------------------------------------------------------
demog <- readr::read_csv(file.path(PATHS$data, "demog.csv"))

demog_both <- demog %>% filter(subji %in% as.numeric(both_sessions))

summary_all  <- get_demog_summary(demog)
summary_both <- get_demog_summary(demog_both)

save_table(summary_all,  "demographics_all")
save_table(summary_both, "demographics_both")

# ------------------------------ 4. Pain Calibration ---------------------------
Pain_calibration <- readxl::read_excel(file.path(PATHS$data, "pain_calibration.xlsx"))

# Paired t-tests
vas2_test <- t.test(Pain_calibration$vas2_p, Pain_calibration$vas2_k, paired = TRUE)
vas8_test <- t.test(Pain_calibration$vas8_p, Pain_calibration$vas8_k, paired = TRUE)

pain_calib_tbl <- tibble(
  level  = c("VAS2", "VAS8"),
  mean_p = c(mean(Pain_calibration$vas2_p), mean(Pain_calibration$vas8_p)),
  sd_p   = c(sd(Pain_calibration$vas2_p),   sd(Pain_calibration$vas8_p)),
  mean_k = c(mean(Pain_calibration$vas2_k), mean(Pain_calibration$vas8_k)),
  sd_k   = c(sd(Pain_calibration$vas2_k),   sd(Pain_calibration$vas8_k)),
  t      = c(vas2_test$statistic, vas8_test$statistic),
  df     = c(vas2_test$parameter,  vas8_test$parameter),
  p      = c(vas2_test$p.value,    vas8_test$p.value)
) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)), stars = p_stars(p))

save_table(pain_calib_tbl, "pain_calibration_stats")

# ------------------------------ 5. Pain Ratings -------------------------------
Pain_ratings <- readr::read_csv(file.path(PATHS$data, "pain_ratings.csv"))

Pain_ratings_wide <- Pain_ratings %>%
  dplyr::select(subji, rating.response_heat_high_1, rating.response_heat_high_2,
         rating.response_heat_low_1,  rating.response_heat_low_2) %>%
  rename_with(~ gsub("rating.response_heat_", "", .x))

Pain_data <- Pain_ratings_wide %>%
  tidyr::pivot_longer(
    cols = starts_with("high") | starts_with("low"),
    names_to = c("intensity", "session"),
    names_pattern = "(\\w+)_(\\d+)",
    values_to = "rating"
  ) %>%
  mutate(
    subji    = factor(subji),
    session  = factor(if_else(session == "1", "placebo", "ketamine"), levels = c("placebo", "ketamine")),
    intensity = dplyr::recode(intensity, high = "Pain", low = "No pain")
  )

RatingBySessionIntensity <- Pain_data %>%
  dplyr::group_by(session, intensity) %>%
  dplyr::summarise(N = n(), meanRating = mean(rating, na.rm = TRUE),
            seRating = sd(rating, na.rm = TRUE)/sqrt(N), .groups = "drop") %>%
  mutate(intensity = factor(intensity, levels = c("Pain", "No pain")))

# Plot
pain_rating_plot <- ggplot(RatingBySessionIntensity, aes(x = intensity, y = meanRating, fill = session)) +
  geom_col(position = position_dodge(width = 0.9), color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = meanRating - seRating, ymax = meanRating + seRating),
                position = position_dodge(width = 0.9), width = 0.2, size = 0.6) +
  geom_point(data = Pain_data, aes(y = rating), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9),
             size = 1.8, alpha = 0.5, shape = 21, stroke = 0.2) +
  scale_fill_manual(values = c(placebo = "#ADD8E6", ketamine = "#FFB6C1")) +
  labs(y = "Subjective pain rating", x = NULL) +
  theme_classic(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

save_plot(pain_rating_plot, "pain_rating_plot", width = 5.5, height = 4)

# Mixed model
Pain_ratings_m <- afex::mixed(
  rating ~ intensity * session + (intensity + 1 | subji),
  data   = Pain_data,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
  method = "KR"
)

posthoc_pain <- emmeans(Pain_ratings_m, consec ~ session | intensity)

posthoc_contr <- emmeans::contrast(
  posthoc_pain,
  method = "pairwise",
  by = "intensity",
) %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(
    p.value = sprintf("%.4f", p.value),  
    sig     = p_stars(as.numeric(p.value)) # keep stars if you want
  )


save_table(posthoc_contr,  "pain_ratings_posthoc_contrasts")

# ------------------------------ 6. CADSS -------------------------------------
CADSS_wide <- readr::read_csv(file.path(PATHS$data, "CADSS.csv"))

CADSS <- CADSS_wide %>%
  pivot_longer(!subji, names_to = c("scale", "sub_scale", "session", "time_point"), names_sep = "_",
               values_to = "rating") %>%
  mutate(
    sub_scale = dplyr::recode(sub_scale,
                              Amnesiasum = "Amnesia",
                              Depersonalizationsum = "Depersonalization",
                              Derealisationsum = "Derealisation"),
    session   = factor(if_else(session == "2", "Placebo", "Ketamine"), levels = c("Placebo", "Ketamine")),
    time_point = dplyr::recode(time_point,
                               baseline = "Pre-infusion",
                               post     = "Post-bolus",
                               end      = "End of infusion")
  ) %>%
  mutate(time_point = factor(time_point, levels = c("Pre-infusion", "Post-bolus", "End of infusion"))) %>%
  drop_na()


RatingBySessionSubScaleTimePoint <- CADSS %>%
  dplyr::group_by(session, sub_scale, time_point) %>%
  dplyr::summarise(N = n(), meanRating = mean(rating), seRating = sd(rating)/sqrt(N), .groups = "drop")

cadss_plot <- ggplot(CADSS, aes(x = time_point, y = rating, group = subji, color = session)) +
  geom_point(size = 1.5, alpha = 0.6) +
  geom_line(alpha = 0.4) +
  geom_line(data = RatingBySessionSubScaleTimePoint,
            aes(x = time_point, y = meanRating, group = session),
            color = "black", size = 0.9, inherit.aes = FALSE) +
  geom_errorbar(data = RatingBySessionSubScaleTimePoint,
                aes(x = time_point, y = meanRating,
                    ymin = meanRating - seRating, ymax = meanRating + seRating),
                width = 0.15, inherit.aes = FALSE) +
  facet_grid(session ~ sub_scale) +
  scale_color_manual(values = c("lightgray", "darkgray")) +
  labs(y = "Rating", x = "Time point") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank())

save_plot(cadss_plot, "cadss_timecourse", width = 9, height = 6)

CADSS_m <- afex::mixed(
  rating ~ time_point * session * sub_scale + (1 + session | subji),
  data = CADSS,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
  method = "KR"
)


cadss_emm_sess_tp <- emmeans(CADSS_m, ~ session | time_point)

cadss_post_session_means <- summary(cadss_emm_sess_tp) %>% as.data.frame()

cadss_post_session <- contrast(cadss_emm_sess_tp, method = "pairwise", by = "time_point") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(
    p.value = sprintf("%.4f", p.value),  
    sig     = p_stars(as.numeric(p.value)) # keep stars if you want
  )
# Time point differences within each session and sub_scale
cadss_emm_tp_sess_scale <- emmeans(CADSS_m, ~ time_point | session | sub_scale)

cadss_post_tp_means <- summary(cadss_emm_tp_sess_scale) %>% as.data.frame()

cadss_post_tp <- contrast(cadss_emm_tp_sess_scale, method = "pairwise", by = c("session", "sub_scale"), adjust = "fdr") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(
    p.value = sprintf("%.4f", p.value),  
    sig     = p_stars(as.numeric(p.value)) # keep stars if you want
  )

save_table(cadss_post_session_means, "cadss_posthoc_session_means")
save_table(cadss_post_session,       "cadss_posthoc_session_contrasts")
save_table(cadss_post_tp_means,      "cadss_posthoc_timepoint_means")
save_table(cadss_post_tp,            "cadss_posthoc_timepoint_contrasts")
# ------------------------------ 7. Pain–CADSS Correlations --------------------
# Prepare matrices
CADSS_wide <- CADSS_wide %>%
  mutate(
    Amnesia_k_min_p          = CADSS_Amnesiasum_2_post         - CADSS_Amnesiasum_1_post,
    Depersonalization_k_min_p= CADSS_Depersonalizationsum_2_post - CADSS_Depersonalizationsum_1_post,
    Derealisation_k_min_p    = CADSS_Derealisationsum_2_post   - CADSS_Derealisationsum_1_post,
    sum                      = CADSS_Amnesiasum_2_post + CADSS_Depersonalizationsum_2_post + CADSS_Derealisationsum_2_post
  )

CADSS_for_corr <- CADSS_wide %>% dplyr::select(subji, CADSS_Amnesiasum_2_post, CADSS_Depersonalizationsum_2_post,
                                        CADSS_Derealisationsum_2_post, sum)

Pain_ratings_wide <- Pain_ratings_wide %>%
  mutate(
    high_k_min_p      = high_2 - high_1,
    high_low_k_min_p  = (high_2 - low_2) - (high_1 - low_1)
  )

Pain_for_corr <- Pain_ratings_wide %>% dplyr::select(subji, high_k_min_p, high_low_k_min_p, high_2)

Pain_CADSS_for_corr <- inner_join(Pain_for_corr, CADSS_for_corr, by = "subji")

names(Pain_CADSS_for_corr) <- c("subji", "Pain high ketamine-placebo", "Pain high-low ketamine-placebo",
                                "Pain high ketamine", "CADSS amnesia", "CADSS depersonalization",
                                "CADSS derealisation", "CADSS sum")

corr_mat <- Hmisc::rcorr(as.matrix(Pain_CADSS_for_corr[,-1]), type = "spearman")

# Heatmap with ComplexHeatmap
col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#6baed6", "white", "#de2d26"))

ht <- ComplexHeatmap::Heatmap(corr_mat$r, name = "Spearman ρ", col = col_fun,
                              cluster_rows = FALSE, cluster_columns = FALSE,
                              show_row_names = TRUE, show_column_names = TRUE,
                              heatmap_legend_param = list(legend_direction = "horizontal"),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                grid::grid.text(sprintf("%.2f", corr_mat$r[i, j]), x, y, gp = grid::gpar(fontsize = 8))
                              })

# Save heatmap as PDF/PNG
pdf(file.path(PATHS$figures, "corr_heatmap.pdf"), width = 6, height = 5)
ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
dev.off()

png(file.path(PATHS$figures, "corr_heatmap.png"), width = 1800, height = 1500, res = 300)
ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
dev.off()

# ------------------------------ 8. fMRI ROI Analyses -------------------------
ROI_wide <- readr::read_csv(file.path(PATHS$data, "roi_beta_values_by_condition.csv"))

ROI_data <- ROI_wide %>%
  pivot_longer(cols = -c(subji, Condition), names_to = "ROI", values_to = "beta") %>%
  mutate(
    subji    = factor(subji),
    ROI      = factor(ROI),
    intensity= if_else(str_detect(Condition, "high"), "pain", "no-pain"),
    session  = if_else(str_detect(Condition, "ketamine"), "ketamine", "placebo"),
    beta_scale = scale(beta)[,1]
  )

ROI_activation_m <- afex::mixed(beta_scale ~ intensity * session * ROI + (1 | subji),
                                data = ROI_data,
                                control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
                                method = "KR")


# Pain minus no-pain difference per ROI ---------------------------------------
ROI_diff <- ROI_data %>%
  dplyr::filter(intensity %in% c("pain", "no-pain")) %>%
  pivot_wider(id_cols = c(subji, session, ROI), names_from = intensity, values_from = beta) %>%
  mutate(diff_beta = pain - `no-pain`)

pain_diff_model <- afex::mixed(diff_beta ~ session * ROI + (1 | subji),
                               data = ROI_diff,
                               control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
                               method = "KR")


emm_roi <- emmeans(pain_diff_model, ~ session | ROI)
roi_posthoc <- broom::tidy(contrast(emm_roi, method = "pairwise"))
save_table(roi_posthoc, "roi_diff_posthoc")

# ---- ROI labels (old_code = new label) ----
roi_lookup <- c(
  ant_insula_L = "Ant. Insula (L)",
  ant_insula_R = "Ant. Insula (R)",
  dlPFC_R      = "DLPFC (R)",
  dACC         = "Dorsal ACC",
  s1_R         = "S1 (R)",
  s2_L         = "S2 (L)",
  s2_R         = "S2 (R)"
)

# Clean ROI names in the wide → diff data
ROI_diff <- ROI_diff %>%
  dplyr::mutate(
    ROI_clean = dplyr::recode(ROI, !!!roi_lookup),
    session   = factor(session, levels = c("placebo", "ketamine"))
  )

# Summaries
ROI_summary <- ROI_diff %>%
  dplyr::group_by(session, ROI_clean) %>%
  dplyr::summarise(
    mean_diff = mean(diff_beta, na.rm = TRUE),
    se_diff   = sd(diff_beta, na.rm = TRUE) / sqrt(n()),
    .groups   = "drop"
  )

asterisk_pos <- ROI_summary %>%
  dplyr::group_by(ROI_clean) %>%
  dplyr::summarise(y = max(mean_diff + se_diff) + 0.1, .groups = "drop")

# Post-hoc contrasts → tidy
roi_posthoc <- emmeans::contrast(emm_roi, method = "pairwise") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(
    sig = p_stars(p.value),
    ROI_clean = dplyr::recode(ROI, !!!roi_lookup)
  ) %>%
  dplyr::select(ROI_clean, sig) %>%
  dplyr::distinct()

# Join sig stars & y positions
ROI_summary <- ROI_summary %>%
  dplyr::left_join(roi_posthoc,  by = "ROI_clean") %>%
  dplyr::left_join(asterisk_pos, by = "ROI_clean")

roi_labels <- c(
  ant_insula_L = "Ant. Insula (L)",
  ant_insula_R = "Ant. Insula (R)",
  dlPFC_R      = "DLPFC (R)",
  dACC         = "Dorsal ACC",
  s1_R         = "S1 (R)",
  s2_L         = "S2 (L)",
  s2_R         = "S2 (R)"
)


# Plot
roi_plot <- ggplot(ROI_summary, aes(x = session, y = mean_diff, fill = session)) +
  geom_col(width = 0.55, color = "black") +
  geom_errorbar(aes(ymin = mean_diff - se_diff, ymax = mean_diff + se_diff), width = 0.18) +
  geom_text(data = ROI_summary %>% dplyr::filter(!is.na(sig) & sig != ""),
            aes(x = 1.5, y = y, label = sig),
            inherit.aes = FALSE, size = 6) +
  scale_fill_manual(values = c(placebo = "#ADD8E6", ketamine = "#FFB6C1")) +
  facet_wrap(~ ROI_clean, ncol = 4) +
  labs(x = "Session", y = "Beta difference (pain – no-pain)") +
  theme_classic(base_size = 14) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


save_plot(roi_plot, "roi_diff_plot", width = 9, height = 6)

# ------------------------------ 9. Pain vs ROI Corrs --------------------------
# Prepare pain difference data
pain_long_diff <- Pain_ratings_wide %>%
  mutate(
    subji              = as.integer(subji),
    pain_diff_placebo  = high_1 - low_1,
    pain_diff_ketamine = high_2 - low_2
  ) %>%
  dplyr::select(subji, pain_diff_placebo, pain_diff_ketamine) %>%
  pivot_longer(
    cols = starts_with("pain_diff"),
    names_to    = "session",
    names_pattern = "pain_diff_(.*)",
    values_to   = "pain_diff_rating"
  ) %>%
  mutate(session = as.character(session))



roi_diff_long <- ROI_diff %>%
  dplyr::rename(pain_diff_activation = diff_beta) %>%
  mutate(
    subji   = str_remove(as.character(subji), "^0+"),
    subji   = as.integer(subji),
    session = as.character(session)
  ) %>%
  dplyr::select(subji, session, ROI, pain_diff_activation)

pain_roi_diff <- left_join(roi_diff_long, pain_long_diff,
                           by = c("subji", "session"))


cor_results <- pain_roi_diff %>%
  dplyr::group_by(ROI, session) %>%
  dplyr::summarise(
    n    = sum(stats::complete.cases(pain_diff_activation, pain_diff_rating)),
    test = list(suppressWarnings(
      cor.test(pain_diff_activation, pain_diff_rating,
               method = "spearman", exact = FALSE)
    )),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    r  = purrr::map_dbl(test, "estimate"),
    p  = purrr::map_dbl(test, "p.value"),
    sig = p_stars(p)
  ) %>%
  dplyr::select(ROI, session, n, r, p, sig)

save_table(cor_results, "pain_roi_corr_results")


session_r_summary <- cor_results %>%
  dplyr::group_by(session) %>%
  dplyr::summarise(
    n_roi  = dplyr::n(),
    mean_r = mean(r, na.rm = TRUE),
    sd_r   = sd(r,   na.rm = TRUE),
    min_r  = min(r,  na.rm = TRUE),
    max_r  = max(r,  na.rm = TRUE)
  )



corr_plot <- ggplot(pain_roi_diff, aes(x = pain_diff_rating, y = pain_diff_activation, color = session)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, size = 1) +
  facet_wrap(~ ROI, scales = "fixed", ncol = 7, labeller = as_labeller(roi_labels)) +
  scale_color_manual(values = c(placebo = "#136580", ketamine = "#b82339")) +
  labs(x = "Pain rating difference (high - low)", y = "ROI activation difference (high - low)", color = "Session") +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(), strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(corr_plot, "pain_roi_corr_plot", width = 11, height = 6)


# ------------------------------ 10. Pain ROI–Dissociation Correspondence ------------------------------

#  Prepare CADSS difference scores

cadss_diff <- CADSS_wide %>%
  mutate(
    # Strip leading zeros then convert to integer
    subji = as.integer(str_remove(as.character(subji), "^0+")),
    
    # Compute difference scores
    amnesia_diff           = CADSS_Amnesiasum_2_post - CADSS_Amnesiasum_1_post,
    depersonalization_diff = CADSS_Depersonalizationsum_2_post - CADSS_Depersonalizationsum_1_post,
    derealisation_diff     = CADSS_Derealisationsum_2_post   - CADSS_Derealisationsum_1_post
  ) %>%
  dplyr::select(subji, amnesia_diff, depersonalization_diff, derealisation_diff)

# Merge pain activation differences with CADSS deltas, keep only ketamine
pain_dissoc <- pain_roi_diff %>%
  mutate(
    # strip any leading zeros and convert to integer
    subji   = as.integer(str_remove(as.character(subji), "^0+")),
    session = as.character(session)
  ) %>%
  dplyr::filter(session == "ketamine") %>%
  left_join(cadss_diff, by = "subji")





pain_roi_dissoc_results <- pain_dissoc %>%
  # ensure ROI is character
  mutate(ROI = as.character(ROI)) %>%
  # split into a list of data‐frames, one per ROI
  group_split(ROI) %>%
  # map over that list and bind rows
  map_dfr(function(df) {
    # fit the two nested models
    m1   <- lm(pain_diff_activation ~ pain_diff_rating, data = df)
    m2   <- lm(
      pain_diff_activation ~ 
        pain_diff_rating + amnesia_diff + depersonalization_diff + derealisation_diff,
      data = df
    )
    comp <- anova(m1, m2)
    coefs <- coef(summary(m2))
    
    # build a one‐row tibble of results
    tibble(
      ROI       = unique(df$ROI),
      pain_p    = coefs['pain_diff_rating',        'Pr(>|t|)'],
      dep_p     = coefs['depersonalization_diff',  'Pr(>|t|)'],
      der_p     = coefs['derealisation_diff',      'Pr(>|t|)'],
      amn_p     = coefs['amnesia_diff',            'Pr(>|t|)'],
      p_anova   = comp$`Pr(>F)`[2],
      f_stat    = comp$F[2],
      df1       = comp$Df[2],
      df2       = comp$Res.Df[2],
      adj_r2_m1 = summary(m1)$adj.r.squared,
      adj_r2_m2 = summary(m2)$adj.r.squared,
      aic_m1    = AIC(m1),
      aic_m2    = AIC(m2)
    )
  }) %>%
  ungroup() %>%
  # add stars
  mutate(
    pain_sig   = case_when(
      pain_p  < 0.001 ~ "***",
      pain_p  < 0.01  ~ "**",
      pain_p  < 0.05  ~ "*",
      TRUE            ~ ""
    ),
    dissoc_sig = case_when(
      p_anova < 0.001 ~ "***",
      p_anova < 0.01  ~ "**",
      p_anova < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

# inspect & save
print(pain_roi_dissoc_results, n = nrow(pain_roi_dissoc_results))
save_table(pain_roi_dissoc_results, "pain_roi_dissoc_model_comparison")


# ------------------------------ 11. NPS Analyses ------------------------------
NPS_resp <- readr::read_csv(file.path(PATHS$data, "NPS.csv"))

NPS_long <- NPS_resp %>% pivot_longer(!subji, names_to = c("signature", "session", "intensity"), names_sep = "_", values_to = "value")

NPS_m <- afex::mixed(value ~ intensity * session + (intensity + session + 1 | subji),
                     data = NPS_long,
                     control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
                     method = "KR")

save_table(broom.mixed::tidy(NPS_m$full_model, effects = "fixed"), "nps_mixed_model")

# EMMs for NPS (avoid broom::tidy on lists)
NPS_emm_sess_int_obj <- emmeans(NPS_m, ~ session | intensity, adjust = "fdr")
NPS_emm_int_sess_obj <- emmeans(NPS_m, ~ intensity | session, adjust = "fdr")

NPS_emm_sess_int <- summary(NPS_emm_sess_int_obj) %>% as.data.frame()
NPS_emm_int_sess <- summary(NPS_emm_int_sess_obj) %>% as.data.frame()

# Pairwise contrasts 
NPS_contr_sess_int <- contrast(NPS_emm_sess_int_obj, method = "pairwise", by = "intensity", adjust = "fdr") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(sig = p_stars(p.value))

NPS_contr_int_sess <- contrast(NPS_emm_int_sess_obj, method = "pairwise", by = "session", adjust = "fdr") %>%
  summary(infer = TRUE) %>%
  as.data.frame() %>%
  dplyr::mutate(sig = p_stars(p.value))

save_table(NPS_emm_sess_int,    "nps_emm_session_intensity")
save_table(NPS_emm_int_sess,    "nps_emm_intensity_session")
save_table(NPS_contr_sess_int,  "nps_contr_session_within_intensity")
save_table(NPS_contr_int_sess,  "nps_contr_intensity_within_session")

NPS_summary <- NPS_long %>% dplyr::group_by(session, intensity) %>% dplyr::summarise(mean = mean(value), se = sd(value)/sqrt(n()), .groups = "drop")


fill_colors <- c("placebo.low" = "#FFB6C1", "placebo.high" = "#ADD8E6",
                 "ketamine.low" = "white", "ketamine.high" = "white")

NPS_plot <- ggplot(NPS_summary, aes(x = session, y = mean, color = intensity, fill = interaction(session, intensity))) +
  geom_col(position = position_dodge(width = 0.9), size = 0.8, width = 0.7) +
  scale_fill_manual(values = fill_colors) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), position = position_dodge(width = 0.9), width = 0.2) +
  geom_hline(yintercept = 0, size = 0.6) +
  geom_point(data = NPS_long, aes(y = value), position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.9), size = 1.5, alpha = 0.5) +
  labs(y = "Pain Signature", x = NULL) +
  scale_x_discrete(limits = c("placebo", "ketamine")) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(NPS_plot, "nps_plot", width = 5.5, height = 4)

# NPS correlations -------------------------------------------------------------
Pain_ratings_wide <- Pain_ratings_wide %>%
  mutate(subji = sprintf("%04d", subji))  # ensures 4-digit character IDs with leading zeros

Pain_NPS_corr <- merge(NPS_resp, Pain_ratings_wide, by = "subji") %>%
  mutate(high_p_min_k = high_1 - high_2)

nps_corr_df <- tibble(
  session = c("placebo", "ketamine"),
  x = c("NPS_placebo_high", "NPS_ketamine_high"),
  y = c("high_1",           "high_2")
) %>%
  mutate(
    test = purrr::map2(x, y, ~ suppressWarnings(
      cor.test(Pain_NPS_corr[[.x]], Pain_NPS_corr[[.y]], method = "spearman", exact = FALSE)
    )),
    r = purrr::map_dbl(test, "estimate"),
    p = purrr::map_dbl(test, "p.value"),
    sig = p_stars(p)
  ) %>%
  dplyr::select(session, r, p, sig)

save_table(nps_corr_df, "nps_pain_corr_results")


# ------------------------------ 12. Pain NPS–Dissociation Correspondence ------------------------------

# Merge NPS, pain ratings & CADSS (you should have CADSS diffs already in cadss_diff)
nps_data <- Pain_NPS_corr %>%
  mutate(
    subji = as.integer(str_remove(as.character(subji), "^0+"))
  ) %>%
  left_join(cadss_diff, by = "subji")

# Fit nested models and extract stats
nps_dissoc_results <- list(nps_data) %>%
  map_dfr(function(df) {
    # Model 1: NPS ~ pain rating only
    m1 <- lm(NPS_ketamine_high ~ high_2, data = df)
    # Model 2: + dissociation diffs
    m2 <- lm(
      NPS_ketamine_high ~ high_2 +
        amnesia_diff + depersonalization_diff + derealisation_diff,
      data = df
    )
    cmp   <- anova(m1, m2)
    coefs <- coef(summary(m2))
    
    tibble(
      outcome      = "NPS",
      pain_p       = coefs["high_2",                       "Pr(>|t|)"],
      dep_p        = coefs["depersonalization_diff",       "Pr(>|t|)"],
      der_p        = coefs["derealisation_diff",           "Pr(>|t|)"],
      amn_p        = coefs["amnesia_diff",                 "Pr(>|t|)"],
      p_anova      = cmp$`Pr(>F)`[2],
      f_stat       = cmp$F[2],
      df1          = cmp$Df[2],
      df2          = cmp$Res.Df[2],
      adj_r2_m1    = summary(m1)$adj.r.squared,
      adj_r2_m2    = summary(m2)$adj.r.squared,
      aic_m1       = AIC(m1),
      aic_m2       = AIC(m2)
    )
  }) %>%
  # Add significance stars for the omnibus test
  mutate(
    dissoc_sig = case_when(
      p_anova < 0.001 ~ "***",
      p_anova < 0.01  ~ "**",
      p_anova < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

# Inspect & save
print(nps_dissoc_results)
save_table(nps_dissoc_results, "nps_dissoc_model_comparison")


# ------------------------------ 11. Connectivity ------------------------------
conn_path <- file.path(PATHS$data, "within_connectivity.xlsx")
all_conn <- readxl::read_excel(conn_path) %>% janitor::clean_names()

# First row holds session tags (PLC/KET)
session_tags <- all_conn[1, -1] %>% unlist() %>% as.character()
all_conn <- all_conn[-1, ]
colnames(all_conn)[1] <- "subji"

# Define networks (adjust if needed)
networks <- c("DefaultMode", "Salience", "Frontoparietal", "SensoryMotor", "DorsalAttention")
# Repeat for PLC/KET
actual_names <- rep(networks, each = 2)
colnames(all_conn)[-1] <- paste(actual_names, session_tags, sep = "_")

all_conn_long <- all_conn %>%
  mutate(subji = as.numeric(subji)) %>%
  pivot_longer(-subji, names_to = c("network", "session"), names_sep = "_", values_to = "connectivity") %>%
  mutate(connectivity = as.numeric(connectivity), session = dplyr::recode(session, PLC = "Placebo", KET = "Ketamine")) %>%
  mutate(session = factor(session, levels = c("Placebo", "Ketamine")))

conn_m <- afex::mixed(connectivity ~ session * network + (1 | subji),
                      data = all_conn_long,
                      control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e6)),
                      method = "KR")

save_table(broom.mixed::tidy(conn_m$full_model, effects = "fixed"), "connectivity_mixed_model")

conn_ttests_paired <- all_conn_long %>%
  mutate(session = as.character(session)) %>%
  group_split(network) %>%
  map_dfr(function(df) {
    wide <- df %>%
      dplyr::select(subji, session, connectivity) %>%
      pivot_wider(names_from = session, values_from = connectivity) %>%
      drop_na(Placebo, Ketamine)    # use exact capitalization
    
    tt <- t.test(wide$Ketamine, wide$Placebo, paired = TRUE)
    
    broom::tidy(tt) %>%
      mutate(network = df$network[1])
  })

conn_ttests_paired

save_table(conn_ttests_paired, "connectivity_ttests")
 
# Plot all networks

all_networks_plot <- ggplot(all_conn_long, aes(x = session, y = connectivity, group = subji)) +
  geom_point(aes(color = session), position = position_dodge(width = 0.2)) +
  geom_line(aes(group = subji, color = session), alpha = 0.5, position = position_dodge(width = 0.2)) +
  stat_summary(fun = mean, geom = "point", size = 3, aes(group = 1), color = "red", shape = 18, position = position_dodge(width = 0.2)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, aes(group = 1), position = position_dodge(width = 0.2)) +
  labs(y = "Within-Network Connectivity", x = "Condition") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_manual(values = c("Placebo" = "lightgray", "Ketamine" = "darkgray")) +
  facet_wrap(~network, scales = "free_y")

# Display
all_networks_plot



save_plot(all_networks_plot, "connectivity_all_networks", width = 8, height = 6)

# DMN focused correlation with CADSS sum --------------------------------------

all_conn <- all_conn %>%
  mutate(
    DefaultMode_diff = as.numeric(DefaultMode_PLC) - as.numeric(DefaultMode_KET)
  )


DMN_diff <- all_conn %>% dplyr::select(subji, DefaultMode_diff) %>% mutate(subji = as.character(subji))

Pain_CADSS_for_corr <- Pain_CADSS_for_corr %>% mutate(subji = as.character(subji))
Pain_NPS_CADSS <- left_join(Pain_CADSS_for_corr, DMN_diff, by = "subji") %>%
  mutate(DMN_diff = as.numeric(DefaultMode_diff), CADSS_sum = as.numeric(`CADSS sum`))

dmn_test <- suppressWarnings(cor.test(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS$CADSS_sum, method = "spearman"))

r_text <- paste0("ρ = ", round(dmn_test$estimate, 2), p_stars(dmn_test$p.value))

DMN_corr_plot <-ggplot(Pain_NPS_CADSS, aes(x = DefaultMode_diff, y = CADSS_sum)) +
  geom_point(color = "#383637", alpha = 0.6,size=3) +
  geom_smooth(method = "lm", se = TRUE, color = "#383637", size = 2) +
  geom_text(aes(x = 0.03, y = min(CADSS_sum, na.rm = TRUE)*0.95, label = r_text),
            inherit.aes = FALSE, size = 10, color = "#383637") +
  labs(x = "DMN within-network connectivity \n (placebo - ketamine)",
       y = "CADSS sum") +
  ylim(0, NA) +
  theme_minimal(base_size = 24) +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


save_plot(DMN_corr_plot, "dmn_cadss_corr", width = 5.5, height = 4)

# Individual subscales correlations (optional reporting)
subscale_tests <- list(
  amnesia          = cor.test(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS$`CADSS amnesia`,          method = "spearman"),
  depersonalization= cor.test(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS$`CADSS depersonalization`, method = "spearman"),
  derealisation    = cor.test(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS$`CADSS derealisation`,    method = "spearman")
)

subscale_tbl <- tibble(
  subscale = names(subscale_tests),
  r = purrr::map_dbl(subscale_tests, ~ .x$estimate),
  p = purrr::map_dbl(subscale_tests, ~ .x$p.value),
  n = purrr::map_int(names(subscale_tests), function(var) {
    sum(complete.cases(Pain_NPS_CADSS$DMN_diff, Pain_NPS_CADSS[[paste("CADSS", var)]]))
  })
) %>%
  mutate(sig = p_stars(p))
save_table(subscale_tbl, "dmn_cadss_subscale_corrs")


# Correlation between change in DMN connectivity and change in pain ratings between sessions
dmn_pain_diff_corr <- cor.test(
  Pain_NPS_CADSS$DMN_diff,
  Pain_NPS_CADSS$`Pain high ketamine-placebo`,
  method = "spearman"
)

# Correlation between DMN connectivity change and high pain rating during the ketamine session
dmn_high_pain_ket_corr <- cor.test(
  Pain_NPS_CADSS$DMN_diff,
  Pain_NPS_CADSS$`Pain high ketamine`,
  method = "spearman"
)

# Show the results
dmn_pain_diff_corr
dmn_high_pain_ket_corr

# ------------------------------ 12. Session Info ------------------------------
cat("\nAnalysis complete. Tables in:", PATHS$tables, "\nFigures in:", PATHS$figures, "\n")

# END OF SCRIPT ----------------------------------------------------------------
