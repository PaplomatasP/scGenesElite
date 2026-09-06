# Install the Bioconductor dependencies used by scGeneFinder.
# Run cran_pkg.R first so that BiocManager is available.

InstallBioc <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  bioc_packages <- c(
    "scran", "M3Drop", "ComplexHeatmap", "pathview", "ensembldb",
    "celldex", "org.Mm.eg.db", "org.Hs.eg.db", "AnnotationFilter",
    "AnnotationDbi", "twoddpcr", "EnsDb.Mmusculus.v79",
    "EnsDb.Hsapiens.v79", "Biobase", "BiocFileCache", "BiocGenerics",
    "BiocParallel", "BiocStyle", "SingleCellExperiment",
    "SummarizedExperiment", "SingleR", "STRINGdb", "MAST", "DESeq2"
  )

  missing <- bioc_packages[
    !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing) == 0L) {
    message("All Bioconductor packages are already installed.")
    return(invisible(bioc_packages))
  }

  message(
    "Installing missing Bioconductor packages: ",
    paste(missing, collapse = ", ")
  )
  BiocManager::install(missing, ask = FALSE, update = FALSE)

  failed <- missing[
    !vapply(missing, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(failed) > 0L) {
    stop(
      "Bioconductor packages failed to install: ",
      paste(failed, collapse = ", ")
    )
  }

  invisible(bioc_packages)
}

InstallBioc()
