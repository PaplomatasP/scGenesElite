# Each test runs in a fresh R process, from the canonical app directory.
stopifnot(dir.exists("scgenes/tests"))
setwd("scgenes")
rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
tests <- sort(list.files("tests", pattern = "^test_.*\\.R$", full.names = TRUE))
stopifnot(length(tests) > 0L)
for (test in tests) {
  message("Running ", test)
  status <- system2(rscript, c("--vanilla", shQuote(test)))
  if (status != 0L) stop(test, " failed with exit code ", status)
}
message("All ", length(tests), " regression scripts passed.")
