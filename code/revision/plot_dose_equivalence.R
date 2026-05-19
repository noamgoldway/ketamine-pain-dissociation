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

# Session-averaged post-infusion CADSS; n = 28, r = -0.056, p = 0.777, R² = 0.003.
df <- read_csv(data_path, show_col_types = FALSE) %>%
  transmute(
    participant_id,
    weight_kg,
    cadss_post = CADSS_total_post_avg
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

stats_label <- paste0(
  "n = ", n_val, "\n",
  "Pearson r = ", sprintf("%.3f", r_val), "\n",
  "p = ", sprintf("%.3f", p_val), "\n",
  "ns"
)

p <- ggplot(df, aes(x = weight_kg, y = cadss_post)) +
  geom_line(
    data = pred_line,
    aes(x = weight_kg, y = fit),
    inherit.aes = FALSE,
    color = "#CC0000",
    linewidth = 0.85,
    linetype = "dashed"
  ) +
  geom_point(
    shape = 21,
    size = 4,
    stroke = 0.55,
    fill = "#A8D4F0",
    color = "#1A1A1A"
  ) +
  annotate(
    "label",
    x = 43.5,
    y = 29.5,
    label = stats_label,
    hjust = 0,
    vjust = 1,
    fill = "#F5F0E1",
    color = "black",
    label.size = 0.35,
    label.padding = unit(0.45, "lines"),
    size = 3.9
  ) +
  annotate(
    "segment",
    x = 70.5,
    xend = 76.5,
    y = 2.8,
    yend = 2.8,
    color = "#CC0000",
    linewidth = 0.85,
    linetype = "dashed"
  ) +
  annotate(
    "label",
    x = 81.5,
    y = 2.8,
    label = paste0("Linear fit (R\u00b2=", sprintf("%.3f", r2_val), ")"),
    hjust = 1,
    vjust = 0.5,
    fill = "white",
    color = "black",
    label.size = 0.35,
    size = 3.7
  ) +
  scale_x_continuous(
    name = "Body Weight (kg)",
    limits = c(42, 83),
    breaks = seq(45, 80, by = 5),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    name = "Post-Infusion CADSS Total Score",
    limits = c(0, 30),
    breaks = seq(0, 30, by = 5),
    expand = expansion(mult = 0.02)
  ) +
  labs(
    title = "Dose Equivalence: Body Weight vs. Dissociative Symptom Response (Ketamine Post-Infusion)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12.5, margin = margin(b = 8)),
    axis.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "#DDDDDD", linewidth = 0.5, linetype = "dashed"),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "#2D2D2D")
  )

out_primary <- list(
  png = file.path(fig_dir, "Figure_S_dose_equivalence_weight_cadss.png"),
  pdf = file.path(fig_dir, "Figure_S_dose_equivalence_weight_cadss.pdf")
)
out_legacy <- list(
  png = file.path(fig_dir, "Figure_S_calibtemp_pain_slopes.png"),
  pdf = file.path(fig_dir, "Figure_S_calibtemp_pain_slopes.pdf")
)

for (paths in list(out_primary, out_legacy)) {
  ggsave(
    filename = paths[["png"]],
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
  ggsave(
    filename = paths[["pdf"]],
    plot = p,
    width = 8,
    height = 6,
    bg = "white"
  )
}

message(
  "Dose equivalence (CADSS_total_post_avg): n=", n_val,
  ", r=", sprintf("%.3f", r_val),
  ", p=", sprintf("%.3f", p_val),
  ", R2=", sprintf("%.3f", r2_val)
)
message("Wrote: ", out_primary$png)
message("Wrote: ", out_primary$pdf)
