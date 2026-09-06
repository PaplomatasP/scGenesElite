# Install the CRAN dependencies used by scGeneFinder.
# The script is idempotent: packages that can already be loaded are skipped.

InstallCran <- function() {
  cran_packages <- c(
    "shiny", "shinyalert", "shinythemes", "enrichR", "fastshap",
    "SeuratObject", "Seurat", "ggplot2", "gridExtra", "glue",
    "tidyverse", "readr", "shinyWidgets", "shinydashboard", "DT",
    "ggiraph", "visNetwork", "png", "shinyjs", "votesys",
    "shinycustomloader", "igraph", "BiocManager", "remotes", "pkgbuild", "caret",
    "randomForest", "xgboost", "C50", "glmnet", "doParallel", "foreach", "zip"
  )

  missing <- cran_packages[
    !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing) == 0L) {
    message("All CRAN packages are already installed.")
    return(invisible(cran_packages))
  }

  message("Installing missing CRAN packages: ", paste(missing, collapse = ", "))
  install.packages(
    missing,
    repos = "https://cloud.r-project.org",
    dependencies = NA,
    Ncpus = max(1L, parallel::detectCores(logical = FALSE) - 1L)
  )

  failed <- missing[
    !vapply(missing, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(failed) > 0L) {
    stop("CRAN packages failed to install: ", paste(failed, collapse = ", "))
  }

  invisible(cran_packages)
}

InstallCran()
