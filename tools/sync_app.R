# Run from the repository root after editing the canonical application.
stopifnot(dir.exists("scgenes/Scripts"))
files <- c("app.R", "ui.R", unlist(lapply(c("Scripts", "data", "www"), function(d) {
  file.path(d, list.files(file.path("scgenes", d), recursive = TRUE))
})))
for (file in files) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.copy(file.path("scgenes", file), file, overwrite = TRUE))
}
message("Synchronized ", length(files), " application compatibility files.")
