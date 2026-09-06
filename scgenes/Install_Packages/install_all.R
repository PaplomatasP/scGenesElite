# One-command dependency installer for scGeneFinder.
#
# From the project directory run:
#   Rscript Install_Packages/install_all.R
#
# The component scripts install only missing packages. This file then verifies
# every package needed by the application, including transitive packages used
# directly by the server scripts.

options(timeout = max(1800L, getOption("timeout", 60L)))
options(repos = c(CRAN = "https://cloud.r-project.org"))

installer_dir <- if (dir.exists("Install_Packages")) {
  normalizePath("Install_Packages", mustWork = TRUE)
} else {
  normalizePath(".", mustWork = TRUE)
}

message("Installing scGeneFinder dependencies with R ", getRversion(), " ...")

source(file.path(installer_dir, "cran_pkg.R"), local = TRUE)
source(file.path(installer_dir, "bioc_pkg.R"), local = TRUE)
source(file.path(installer_dir, "git_pkg.R"), local = TRUE)

required_packages <- c(
  "shiny", "shinyalert", "shinythemes", "enrichR", "ggplot2",
  "gridExtra", "glue", "tidyverse", "shinyWidgets", "shinydashboard",
  "twoddpcr", "SCMarker", "scran", "DT", "Seurat", "pathview", "png",
  "ggiraph", "AnnotationDbi", "AnnotationFilter", "Biobase",
  "BiocFileCache", "BiocGenerics", "BiocParallel", "BiocStyle",
  "BiocManager", "fastAdaboost", "votesys", "M3Drop", "ComplexHeatmap",
  "igraph", "visNetwork", "SingleR", "shinyjs", "STRINGdb", "fastshap",
  "xgboost", "randomForest", "SingleCellExperiment", "SummarizedExperiment",
  "MAST", "DESeq2", "BPSC", "scPNMF", "SelfE", "org.Hs.eg.db",
  "caret", "zip", "foreach", "doParallel"
)

failed <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(failed) > 0L) {
  stop(
    "Dependency installation is incomplete. Packages that still cannot load: ",
    paste(failed, collapse = ", ")
  )
}

message(
  "Success: all ", length(required_packages),
  " scGeneFinder dependency packages can be loaded."
)
