library(shiny)
source("Scripts/InputValidation.R")
source("Scripts/UploadPreview.R")

testServer(function(input, output, session) {
  register_upload_preview(input, output)
}, {
  stopifnot(identical(output$uploadState, "empty"))
  session$setInputs(
    file1 = data.frame(name = "ExampleData.csv", size = file.info("data/ExampleData.csv")$size,
                      datapath = normalizePath("data/ExampleData.csv")),
    header = TRUE, sep = ",", quote = '"', disp = "head"
  )
  stopifnot(identical(output$uploadState, "ready"))
  stopifnot(grepl("Labels", output$contents, fixed = TRUE))
  session$setInputs(file1 = NULL)
  stopifnot(identical(output$uploadState, "empty"))

  # Use a freshly serialized fixture to test the RDS path independently.
  fixture <- tempfile(fileext = ".rds")
  saveRDS(data.frame(GeneA = 1:4, GeneB = 5:8,
                     Labels = c("control", "control", "case", "case")), fixture)
  session$setInputs(rdsFile = data.frame(name = "fixture.rds", size = file.info(fixture)$size,
                                        datapath = fixture))
  stopifnot(identical(output$uploadState, "ready"))
  stopifnot(grepl("GeneA", output$Rvalue, fixed = TRUE))
  session$setInputs(rdsFile = NULL)
  stopifnot(identical(output$uploadState, "empty"))
  unlink(fixture)
})
cat("Upload preview checks passed for CSV, RDS, and cleared uploads.\n")
