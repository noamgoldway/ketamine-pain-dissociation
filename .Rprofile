# Activate renv only when the project library is populated (after renv::restore()).
if (dir.exists("renv/library") && length(list.dirs("renv/library", recursive = FALSE)) > 0) {
  source("renv/activate.R")
}
