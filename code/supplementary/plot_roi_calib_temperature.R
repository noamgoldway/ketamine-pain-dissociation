#!/usr/bin/env Rscript
# Supplementary Figure S2: mean z-scored pain-ROI activation vs calibrated pain temperature.
# Slopes match submitted figure (placebo ≈ 0.341, ketamine ≈ 0.131).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd[grep("^--file=", cmd)])
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(file_arg[[1]]), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

pain_rois <- c(
  "ant_insula_L", "ant_insula_R", "dlPFC_R", "dACC",
  "s1_R", "s2_L", "s2_R"
)

roi <- read_csv(
  file.path(repo_root, "data", "roi_beta_values_by_condition.csv"),
  show_col_types = FALSE
)
cal <- readxl::read_excel(file.path(repo_root, "data", "pain_calibration.xlsx"))

cal_long <- cal %>%
  mutate(subji = as.integer(subji)) %>%
  pivot_longer(c(vas8_p, vas8_k), names_to = "key", values_to = "calib_temp") %>%
  mutate(session = if_else(key == "vas8_p", "placebo", "ketamine")) %>%
  select(subji, session, calib_temp)

plot_df <- roi %>%
  filter(str_detect(Condition, "pain_high")) %>%
  mutate(
    subji = as.integer(str_remove(as.character(subji), "^0+")),
    session = if_else(str_detect(Condition, "ketamine"), "ketamine", "placebo")
  ) %>%
  pivot_longer(all_of(pain_rois), names_to = "ROI", values_to = "beta") %>%
  group_by(ROI, session) %>%
  mutate(beta_z = as.numeric(scale(beta)[, 1])) %>%
  ungroup() %>%
  group_by(subji, session) %>%
  summarise(mean_beta_z = mean(beta_z), .groups = "drop") %>%
  left_join(cal_long, by = c("subji", "session"))

session_slopes <- plot_df %>%
  group_by(session) %>%
  group_modify(function(d, ...) {
    fit <- lm(mean_beta_z ~ calib_temp, data = d)
    tibble(
      calib_temp.trend = unname(coef(fit)[2]),
      intercept = unname(coef(fit)[1]),
      n = nrow(d)
    )
  }) %>%
  ungroup()

tab_dir <- file.path(repo_root, "output", "supplementary", "tables")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(session_slopes, file.path(tab_dir, "supp_fig_s2_session_slopes.csv"))

slope_labels <- session_slopes %>%
  mutate(
    label = sprintf(
      "slope = %.3f",
      calib_temp.trend
    ),
    x = 42.2,
    y = if_else(session == "placebo", 1.55, 1.25),
    color = if_else(session == "placebo", "#136580", "#b82339")
  )

p <- ggplot(plot_df, aes(x = calib_temp, y = mean_beta_z, color = session)) +
  geom_point(alpha = 0.55, size = 3) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1.1) +
  geom_text(
    data = slope_labels,
    aes(x = x, y = y, label = label, color = session),
    inherit.aes = FALSE,
    hjust = 0,
    size = 5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(placebo = "#136580", ketamine = "#b82339"),
    labels = c(placebo = "Placebo", ketamine = "Ketamine")
  ) +
  scale_x_continuous(
    name = "Individually calibrated pain temperature (°C)",
    breaks = seq(42, 48, by = 2)
  ) +
  scale_y_continuous(
    name = "ROI activation (z-scored β)",
    limits = c(-2, 2),
    breaks = seq(-2, 2, by = 1)
  ) +
  theme_minimal(base_size = 16) +
  theme(
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  )

fig_dir <- file.path(repo_root, "output", "supplementary", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(
  file.path(fig_dir, "Supplementary_Figure_S2_roi_calib_temperature.png"),
  p,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

message(
  "S2 slopes: placebo=", sprintf("%.3f", session_slopes$calib_temp.trend[session_slopes$session == "placebo"]),
  ", ketamine=", sprintf("%.3f", session_slopes$calib_temp.trend[session_slopes$session == "ketamine"])
)
