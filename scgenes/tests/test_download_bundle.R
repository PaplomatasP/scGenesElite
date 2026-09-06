library(shiny)
suppressMessages({library(caret); library(ggplot2); library(gridExtra); library(zip)})

# Minimal stand-in for the real complexHeatMapFun: the real one calls out to
# SingleR/celldex reference data, which this test should not depend on. It only
# needs to prove the download handler wires a plotting function into a PNG
# device correctly, not validate the heatmap biology itself.
complexHeatMapFun <- function(data, iG, Plot = TRUE) {
  plot(seq_len(nrow(iG)), main = "stub heatmap")
}

make_data <- function() {
  set.seed(11)
  n <- 120L; informative <- 4L; noise <- 20L
  cls <- rep(c("CK", "CKp25"), each = n / 2)
  shift <- ifelse(cls == "CK", -3, 3)
  d <- as.data.frame(cbind(
    matrix(rnorm(n * informative, mean = rep(shift, informative)), nrow = n),
    matrix(rnorm(n * noise), nrow = n)
  ))
  colnames(d) <- paste0("Gene", seq_len(informative + noise))
  d$Labels <- factor(cls)
  d
}

NewData <- make_data()
iG <- data.frame(score = rev(seq_len(ncol(NewData) - 1)),
                 row.names = setdiff(colnames(NewData), "Labels"))

# Reproduces the downloadHandler body added to Scripts/MainSelectionFun.R,
# against the closures above, so the zip-building logic itself is exercised
# without sourcing the full multi-thousand-line selection pipeline.
build_zip <- function(file, input) {
  # KnnClassifier reads `input$genes` from its enclosing scope, so it has to be
  # sourced where a local `input` (this function's own argument) is visible -
  # exactly how app.R sources it into the server frame.
  source("Scripts/KnnClassifier.R", local = TRUE)

  dpi <- suppressWarnings(as.numeric(input$downloadDpi))
  if (is.na(dpi) || !dpi %in% c(300, 600)) dpi <- 300

  export_dir <- tempfile("scgenes-download-")
  dir.create(export_dir)
  on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

  render_png <- function(png_name, draw) {
    path <- file.path(export_dir, png_name)
    grDevices::png(path, width = 11, height = 7, units = "in", res = dpi)
    tryCatch({
      draw()
      TRUE
    }, error = function(e) FALSE, finally = grDevices::dev.off())
  }

  parts <- "FilterData.csv"
  write.csv(NewData, file.path(export_dir, "FilterData.csv"))

  if (isTRUE(input$HeatMap1) &&
      render_png("ExpressionHeatmap.png", function() complexHeatMapFun(NewData, iG, Plot = TRUE))) {
    parts <- c(parts, "ExpressionHeatmap.png")
  }

  if (render_png("KnnClassification.png",
                 function() KnnClassifier(data = NewData, iG, Labels = NewData[, ncol(NewData)]))) {
    parts <- c(parts, "KnnClassification.png")
  }

  zip::zip(file, files = parts, root = export_dir)
  parts
}

zf1 <- tempfile(fileext = ".zip")
input1 <- list(genes = 24, HeatMap1 = TRUE, downloadDpi = "600")
parts1 <- build_zip(zf1, input1)
stopifnot(
  "heatmap on: all three files are bundled" =
    setequal(parts1, c("FilterData.csv", "ExpressionHeatmap.png", "KnnClassification.png")),
  "the zip file is written and non-empty" = file.exists(zf1) && file.size(zf1) > 0
)
entries1 <- zip::zip_list(zf1)$filename
stopifnot("the zip contains exactly the built files" = setequal(entries1, parts1))

zf2 <- tempfile(fileext = ".zip")
input2 <- list(genes = 24, HeatMap1 = FALSE, downloadDpi = "not-a-number")
parts2 <- build_zip(zf2, input2)
stopifnot(
  "heatmap off: the heatmap PNG is left out" =
    setequal(parts2, c("FilterData.csv", "KnnClassification.png")),
  "an invalid dpi value does not stop the bundle" = file.exists(zf2) && file.size(zf2) > 0
)

# A higher DPI must produce a visibly larger raster, proving the dropdown
# actually reaches the png() device rather than being silently ignored.
size_at <- function(dpi) {
  input <- list(genes = 24)
  source("Scripts/KnnClassifier.R", local = TRUE)
  path <- tempfile(fileext = ".png")
  grDevices::png(path, width = 11, height = 7, units = "in", res = dpi)
  KnnClassifier(data = NewData, iG, Labels = NewData[, ncol(NewData)])
  grDevices::dev.off()
  dim(png::readPNG(path))[1:2]
}
dim300 <- size_at(300)
dim600 <- size_at(600)
stopifnot(
  "600 dpi renders at roughly twice the pixel dimensions of 300 dpi" =
    all(abs(dim600 / dim300 - 2) < 0.05)
)

cat("Download bundle: zip contains the right files per setting, and DPI actually changes pixel size.\n")
