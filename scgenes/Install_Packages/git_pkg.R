# Install research packages that are distributed from their upstream GitHub
# repositories rather than Bioconductor or the current CRAN repository.

InstallGit <- function() {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }

  github_packages <- c(
    fastAdaboost = "souravc83/fastAdaboost",
    scPNMF = "JSB-UCLA/scPNMF",
    SelfE = "Priyadarshini-Rai/SelfE",
    BPSC = "nghiavtr/BPSC",
    SCMarker = "KChen-lab/SCMarker"
  )

  missing <- names(github_packages)[
    !vapply(names(github_packages), requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing) == 0L) {
    message("All GitHub packages are already installed.")
    return(invisible(github_packages))
  }

  if (.Platform$OS.type == "windows" &&
      requireNamespace("pkgbuild", quietly = TRUE) &&
      !pkgbuild::has_build_tools(debug = FALSE)) {
    stop(
      "Rtools is required for the missing GitHub packages. ",
      "Run Install_Packages/install_all.ps1 to install Rtools and all packages automatically."
    )
  }

  for (package in missing) {
    repository <- unname(github_packages[[package]])
    message("Installing ", package, " from GitHub repository ", repository)
    remotes::install_github(
      repository,
      dependencies = NA,
      upgrade = "never",
      build_vignettes = FALSE,
      force = TRUE
    )
  }

  failed <- missing[
    !vapply(missing, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(failed) > 0L) {
    stop("GitHub packages failed to install: ", paste(failed, collapse = ", "))
  }

  invisible(github_packages)
}

InstallGit()
