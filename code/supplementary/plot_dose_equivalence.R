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
fig_dir <- file.path(repo_root, "output", "supplementary", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Ketamine post-bolus CADSS (session 2); matches manuscript caption (r ~ 0.04, p ~ 0.86).
df <- read_csv(data_path, show_col_types = FALSE) %>%
  transmute(
    participant_id,
    weight_kg,
    cadss_post = CADSS_total_post_session2
  )

ct <- cor.test(df$weight_kg, df$cadss_post)
fit <- lm(cadss_post ~ weight_kg, data = df)
fit_sum <- summary(fit)

r_val <- unname(ct$estimate)
p_val <- ct$p.value
r2_val <- fit_sum$r.squared
n_val <- nrow(df)

x_rng <- range(df$weight_kg, na.rm = TRUE)
pred_line <- data.frame(weight_kg = seq(x_rng[1], x_rng[2], length.out = 100))
pred_line$fit <- predict(fit, newdata = pred_line)
pred_ci <- as.data.frame(predict(fit, newdata = pred_line, interval = "confidence"))
pred_line$ymin <- pred_ci$lwr
pred_line$ymax <- pred_ci$upr

p <- ggplot(df, aes(x = weight_kg, y = cadss_post)) +
  geom_ribbon(
    data = pred_line,
    aes(x = weight_kg, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#D88FA1",
    alpha = 0.22
  ) +
  geom_line(
    data = pred_line,
    aes(x = weight_kg, y = fit),
    inherit.aes = FALSE,
    color = "#B12243",
    linewidth = 1.1
  ) +
  geom_point(
    shape = 21,
    size = 3.2,
    stroke = 0.9,
    fill = "#2C7F95",
    color = "#1E5D6E",
    alpha = 0.55
  ) +
  annotate(
    "text",
    x = max(df$weight_kg),
    y = min(df$cadss_post) + 1.5,
    label = paste0("Linear fit (R\u00b2=", sprintf("%.3f", r2_val), ")"),
    hjust = 1,
    vjust = 0,
    size = 4.2,
    color = "#444444"
  ) +
  scale_x_continuous(
    name = "Body Weight (kg)",
    breaks = seq(50, 80, by = 10)
  ) +
  scale_y_continuous(
    name = "Post-Bolus CADSS Total Score",
    limits = c(0, 50),
    breaks = seq(0, 50, by = 10)
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

stem <- "Supplementary_Figure_S1_weight_cadss_submitted"
ggsave(
  filename = file.path(fig_dir, paste0(stem, ".png")),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

message(
  "Dose equivalence (CADSS_total_post_session2): n=", n_val,
  ", r=", sprintf("%.3f", r_val),
  ", p=", sprintf("%.3f", p_val),
  ", R2=", sprintf("%.3f", r2_val)
)
message("Wrote figures under: ", fig_dir)
