# Render a KEGG pathway in an isolated temporary directory.
#
# pathview downloads the KEGG XML/PNG files into kegg.dir, but its native
# renderer writes the final image into getwd(). Shiny sessions share that
# working directory, so every invocation is moved to its own temporary
# directory and uses a deterministic output suffix.
plot_pathview <- function(
    ...,
    save_image = FALSE,
    .pathview_runner = pathview::pathview,
    .display_image = TRUE) {
  if (!"package:pathview" %in% search()) {
    suppressPackageStartupMessages(
      require("pathview", character.only = TRUE, quietly = TRUE)
    )
  }

  args <- list(...)
  pathway_id <- tolower(trimws(as.character(args$pathway.id)))
  species <- tolower(trimws(as.character(args$species)))

  if (length(pathway_id) != 1L || length(species) != 1L ||
      !species %in% c("hsa", "mmu")) {
    stop("A single valid species ('hsa' or 'mmu') is required.", call. = FALSE)
  }

  pathway_id <- sub(paste0("^", species), "", pathway_id)
  if (!grepl("^[0-9]{5}$", pathway_id)) {
    stop(
      "The KEGG pathway ID must contain exactly five digits, for example 04110.",
      call. = FALSE
    )
  }

  pathway_name <- paste0(species, pathway_id)
  output_suffix <- "scgenes"
  run_directory <- tempfile("scgenes-pathview-")
  if (!dir.create(run_directory, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create an isolated directory for the KEGG pathway.", call. = FALSE)
  }
  run_directory <- normalizePath(run_directory, winslash = "/", mustWork = TRUE)
  temporary_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  if (!startsWith(run_directory, paste0(temporary_root, "/"))) {
    stop("Refusing to use a KEGG working directory outside the R temporary directory.", call. = FALSE)
  }

  original_directory <- getwd()
  cleanup_directory <- !isTRUE(save_image)
  completed <- FALSE
  on.exit({
    try(setwd(original_directory), silent = TRUE)
    if ((!completed || cleanup_directory) && dir.exists(run_directory)) {
      unlink(run_directory, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  args$pathway.id <- pathway_id
  args$species <- species
  args$kegg.dir <- run_directory
  args$out.suffix <- output_suffix

  # pathview's native renderer writes relative to getwd(). This is synchronous
  # in Shiny's R process, and the original directory is restored by on.exit().
  setwd(run_directory)
  do.call(.pathview_runner, args)

  png_files <- list.files(
    run_directory,
    pattern = "[.]png$",
    full.names = TRUE
  )
  png_names <- basename(png_files)
  expected_prefix <- paste0(pathway_name, ".", output_suffix)
  output_files <- png_files[
    startsWith(png_names, expected_prefix) & endsWith(png_names, ".png")
  ]

  if (length(output_files) == 0L) {
    stop(
      paste0(
        "Pathview did not create the expected pathway image for ",
        pathway_name,
        ". Check the pathway ID, organism, and gene IDs."
      ),
      call. = FALSE
    )
  }

  if (length(output_files) > 1L) {
    output_files <- output_files[which.max(file.info(output_files)$mtime)]
  }

  image <- png::readPNG(output_files[[1]])
  if (isTRUE(.display_image)) {
    grid::grid.raster(image)
  }

  completed <- TRUE
  invisible(if (isTRUE(save_image)) output_files[[1]] else NULL)
}
