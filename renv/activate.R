
local({

  # the requested version of renv
  version <- "1.0.7"
  attr(version, "sha") <- NULL

  # the project directory
  project <- getwd()

  # use start-up diagnostics if enabled
  diagnostics <- Sys.getenv("RENV_STARTUP_DIAGNOSTICS", unset = "FALSE")
  diagnostics <- identical(diagnostics, "TRUE")

  # load renv from the project
  if (!requireNamespace("renv", quietly = TRUE)) {
    message("renv not installed; install with install.packages('renv')")
    return(invisible(FALSE))
  }

  renv::load(project)

  invisible(TRUE)

})
