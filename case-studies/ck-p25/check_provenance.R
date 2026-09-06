# Run from the repository root. Optionally pass the public source file path.
library(data.table)
out <- "case-studies/ck-p25"
args <- commandArgs(trailingOnly = TRUE)
source_file <- if (length(args)) args[[1]] else file.path(out, "cache/GSE103334_FPKM_CKP25_TOPHAT.txt.gz")
if (!file.exists(source_file)) {
  if (length(args)) stop("Supplied source file does not exist: ", source_file)
  dir.create(dirname(source_file), recursive = TRUE, showWarnings = FALSE)
  options(timeout = max(600, getOption("timeout")))
  download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE103nnn/GSE103334/suppl/GSE103334_FPKM_CKP25_TOPHAT.txt.gz", source_file, mode = "wb")
}
src <- fread(source_file)
d <- read.csv(file.path(out, "input_ExampleData.csv"), row.names = 1, check.names = FALSE)
idx <- match(names(d)[-ncol(d)], make.names(src[[1]], unique = TRUE))
stopifnot(!anyNA(idx), all(rownames(d) %in% names(src)))
a <- as.matrix(d[, -ncol(d)])
b <- t(as.matrix(src[idx, rownames(d), with = FALSE]))
stopifnot(max(abs(a - b)) < 1e-12)
message("Verified ", length(a), " values against GSE103334; maximum difference = ", max(abs(a - b)))
