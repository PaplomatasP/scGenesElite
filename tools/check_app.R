stopifnot(dir.exists("scgenes/Scripts"))
files <- c("app.R", "ui.R", unlist(lapply(c("Scripts", "data", "www"), function(d) {
  file.path(d, list.files(file.path("scgenes", d), recursive = TRUE))
})))
for (file in files) {
  stopifnot(file.exists(file))
  if (unname(tools::md5sum(file)) != unname(tools::md5sum(file.path("scgenes", file)))) {
    stop("Compatibility copy differs from scgenes/: ", file)
  }
}
rfiles <- list.files("scgenes", pattern = "\\.[Rr]$", recursive = TRUE, full.names = TRUE)
invisible(lapply(rfiles, parse))
message("Parsed ", length(rfiles), " R files; ", length(files), " compatibility files match.")
