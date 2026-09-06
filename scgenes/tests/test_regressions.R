source("Scripts/InputValidation.R")

assert_true <- function(value, message) {
  if (!isTRUE(value)) {
    stop(message, call. = FALSE)
  }
}

# Real example files must parse to the same validated shape. The CSV contains
# a leading row-ID column, which should become row names rather than a gene.
rds_upload <- list(
  name = "ExampleData.rds",
  size = unname(file.info("data/ExampleData.rds")$size),
  datapath = "data/ExampleData.rds"
)
csv_upload <- list(
  name = "ExampleData.csv",
  size = unname(file.info("data/ExampleData.csv")$size),
  datapath = "data/ExampleData.csv"
)
rds_data <- read_uploaded_rds(rds_upload)
csv_data <- read_uploaded_csv(csv_upload)
assert_true(nrow(rds_data) == nrow(csv_data), "CSV/RDS row counts differ.")
assert_true(ncol(rds_data) == ncol(csv_data), "CSV row-ID normalization failed.")
assert_true(inspect_expression_dataset(rds_data)$valid, "Example RDS validation failed.")

small_data <- data.frame(
  GeneA = 1:4,
  GeneB = 5:8,
  Labels = c("control", "control", "case", "case"),
  check.names = FALSE
)
small_preview <- preview_expression_data(small_data)
assert_true(
  identical(dim(small_preview), c(4L, 3L)),
  "Preview must support datasets with fewer than 11 columns."
)

duplicate_genes <- small_data
colnames(duplicate_genes)[1:2] <- c("GeneA", "genea")
assert_true(
  !inspect_expression_dataset(duplicate_genes)$valid,
  "Case-insensitive duplicate gene IDs were accepted."
)

nonnumeric_gene <- small_data
nonnumeric_gene$GeneA <- c("x", "y", "z", "q")
assert_true(
  !inspect_expression_dataset(nonnumeric_gene)$valid,
  "A nonnumeric expression column was accepted."
)

empty_label <- small_data
empty_label$Labels[4] <- ""
assert_true(
  !inspect_expression_dataset(empty_label)$valid,
  "An empty label was accepted."
)

single_class <- small_data
single_class$Labels <- "control"
assert_true(
  !inspect_expression_dataset(single_class)$valid,
  "A single-class dataset was accepted."
)

assert_true(
  !inspect_expression_dataset(small_data, max_expression_values = 1)$valid,
  "The expression-value dimension limit was not enforced."
)

# Ensemble regression: a separate SHAP choice must execute and contribute to
# the vote alongside a wrapper method.
ensemble_env <- new.env(parent = globalenv())
ensemble_env$input <- list(
  ensembleVar = "NoMethod",
  ensemblePvalue = "NoMethod",
  ensembleWrapper = "rf",
  ensembleSHAP = "shap_rf",
  SHAP_importanceLimit = 0.01,
  importanceLimit = 0.01,
  VarFilter = "Unselect"
)
ensemble_env$shap_calls <- 0L
ensemble_env$SelectionFilter1 <- function(data, Labels, MLmethod) {
  list(newdata = data.frame(GeneA = data$GeneA, GeneB = data$GeneB, Labels = Labels))
}
ensemble_env$ShapValuesFilter <- function(data, Labels, MLmethod, importanceLimit) {
  ensemble_env$shap_calls <- ensemble_env$shap_calls + 1L
  list(newdata = data.frame(GeneB = data$GeneB, GeneC = data$GeneC, Labels = Labels))
}
ensemble_env$create_vote <- function(GenesList, xtype, candidate) GenesList
ensemble_env$borda_method <- function(vote, modified) {
  list(other_info = list(count_max = c(GeneA = 1, GeneB = 2, GeneC = 1)))
}
ensemble_env$shinyalert <- function(...) invisible(NULL)
sys.source("Scripts/EnsemleMethod.R", envir = ensemble_env)

ensemble_input <- data.frame(
  GeneA = 1:4,
  GeneB = 2:5,
  GeneC = 3:6
)
ensemble_result <- ensemble_env$EnsemleMethod(
  ensemble_input,
  c("control", "control", "case", "case")
)
assert_true(ensemble_env$shap_calls == 1L, "The ensemble SHAP method did not execute.")
assert_true(
  is.list(ensemble_result) && !is.null(ensemble_result$ig),
  "The ensemble SHAP result did not contribute to the vote."
)

# KEGG regression: a fake runner verifies cwd isolation and deterministic file
# selection without using the network.
source("Scripts/pathway.R")
original_directory <- normalizePath(getwd(), winslash = "/")
fake_pathview <- function(pathway.id, species, kegg.dir, out.suffix, ...) {
  png::writePNG(
    array(1, dim = c(10, 10, 4)),
    target = paste0(species, pathway.id, ".", out.suffix, ".png")
  )
  invisible(list())
}
isolated_image <- plot_pathview(
  gene.data = c("7157" = 1),
  pathway.id = "04110",
  species = "hsa",
  save_image = TRUE,
  .pathview_runner = fake_pathview,
  .display_image = FALSE
)
assert_true(file.exists(isolated_image), "The isolated KEGG image was not created.")
assert_true(
  normalizePath(getwd(), winslash = "/") == original_directory,
  "plot_pathview did not restore the working directory."
)
assert_true(
  normalizePath(dirname(isolated_image), winslash = "/") != original_directory,
  "The KEGG image was written into the shared application directory."
)
unlink(dirname(isolated_image), recursive = TRUE, force = TRUE)

app_source <- paste(readLines("app.R", warn = FALSE), collapse = "\n")
ensemble_source <- paste(readLines("Scripts/EnsemleMethod.R", warn = FALSE), collapse = "\n")
statistical_source <- paste(readLines("Scripts/StatisticalPvalueFilter.R", warn = FALSE), collapse = "\n")
assert_true(!grepl("Human_Phenoty pe_Ontology", app_source, fixed = TRUE), "Ontology typo remains.")
assert_true(!grepl("exists(\"input$", ensemble_source, fixed = TRUE), "Broken input exists() check remains.")
assert_true(!grepl("monocle_method", ensemble_source, fixed = TRUE), "Obsolete ensemble method mapping remains.")
assert_true(
  grepl('if (PvalueMethod == "DESeq2_method")', statistical_source, fixed = TRUE),
  "DESeq2 ensemble post-processing still checks the wrong input."
)

cat("All regression tests passed.\n")
