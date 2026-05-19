#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd[grep("^--file=", cmd)])
repo_root <- if (length(file_arg)) {
  normalizePath(file.path(dirname(file_arg[[1]]), "..", ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
data_path <- file.path(repo_root, "data", "CADSS_Weight_DoseEquivalence_Data.csv")
fig_dir <- file.path(repo_root, "output", "revision", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

df <- read_csv(data_path, show_col_types = FALSE) %>%
  transmute(
    participant_id,
    weight_kg,
    cadss_post_ketamine = CADSS_total_post_session2
  )

fit <- lm(cadss_post_ketamine ~ weight_kg, data = df)
fit_sum <- summary(fit)
r_val <- unname(cor(df$weight_kg, df$cadss_post_ketamine, use = "complete.obs"))
p_val <- coef(fit_sum)[2, 4]
r2_val <- fit_sum$r.squared
n_val <- nrow(df)

p <- ggplot(df, aes(x = weight_kg, y = cadss_post_ketamine)) +
  geom_point(
    shape = 21,
    size = 3.2,
    stroke = 0.9,
    fill = "#2C7F95",
    color = "#1E5D6E",
    alpha = 0.55
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    color = "#B12243",
    fill = "#D88FA1",
    linetype = "solid",
    linewidth = 1.2,
    alpha = 0.18
  ) +
  annotate(
    "text",
    x = max(df$weight_kg),
    y = min(df$cadss_post_ketamine) + 0.8,
    label = paste0("Linear fit  (R^2=", sprintf("%.3f", r2_val), ")"),
    hjust = 1,
    vjust = 0,
    size = 4.2,
    color = "#444444"
  ) +
  labs(
    x = "Body Weight (kg)",
    y = "Post-Bolus CADSS Total Score"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "#DDDDDD", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.text = element_text(color = "#2D2D2D")
  )

ggsave(
  filename = file.path(fig_dir, "Figure_S_calibtemp_pain_slopes.png"),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path(fig_dir, "Figure_S_calibtemp_pain_slopes.pdf"),
  plot = p,
  width = 10,
  height = 7
)
