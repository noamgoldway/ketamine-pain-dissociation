#!/usr/bin/env Rscript
# Export Supplementary Table S1 (connectivity ROI atlas) as tidy CSV.
# Run after: make extract-claims

suppressPackageStartupMessages(library(jsonlite))

repo_root <- Sys.getenv("ROOT_DIR", unset = "")
if (!nzchar(repo_root)) repo_root <- getwd()
repo_root <- normalizePath(repo_root, mustWork = FALSE)

claims_path <- file.path(repo_root, "docs", "extracted_claims_r2.json")
out_path <- file.path(repo_root, "data", "connectivity_roi_list.csv")

if (!file.exists(claims_path)) {
  stop("Missing ", claims_path, " — run: make extract-claims")
}

claims <- fromJSON(claims_path, simplifyVector = FALSE)
rows <- claims$supplementary_tables$S1$rows
if (is.null(rows) || length(rows) < 3) {
  stop("S1 table missing or too short in extracted claims")
}

# Skip two header rows; data rows: label, X, Y, Z, network
data_rows <- rows[(seq_along(rows) > 2)]
lines <- c("roi_label,x,y,z,network")
for (row in data_rows) {
  if (length(row) < 5) next
  label <- gsub(",", ";", row[[1]])
  x <- row[[2]]
  y <- row[[3]]
  z <- row[[4]]
  net <- row[[5]]
  if (!nzchar(label) || label == "ROI Label") next
  lines <- c(lines, paste(label, x, y, z, net, sep = ","))
}

writeLines(lines, out_path, useBytes = TRUE)
message("Wrote ", out_path, " (", length(lines) - 1, " ROI rows)")
