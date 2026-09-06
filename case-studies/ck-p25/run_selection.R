# Run the repository's SCMarker implementation on the user's supplied example.
# Run from the scgenes project root. No changes to application source are needed.
library(SCMarker)
out <- "case-studies/ck-p25"
source(file.path(out,"code_snapshot/InputValidation.R"))
source(file.path(out,"code_snapshot/Lexikon.R"))
source(file.path(out,"code_snapshot/VariableGenesMethods.R"))
set.seed(20260906)
d <- read_uploaded_csv(list(name = "ExampleData.csv",
  size = file.info(file.path(out, "input_ExampleData.csv"))$size,
  datapath = file.path(out, "input_ExampleData.csv")))
mapped <- LexikonFun(d, "Mouse", "SYMBOL")
stopifnot(all(vapply(names(mapped), function(g) identical(mapped[[g]], d[[g]]), logical(1))))
write.csv(data.frame(gene = setdiff(names(d), names(mapped))),
          file.path(out, "unmapped_input_genes.csv"), row.names = FALSE)
metadata <- data.frame(cell = rownames(mapped), group = mapped$Labels,
  week = sub("^.*_([0-9]+)w_.*$", "\\1", rownames(mapped)),
  sample = sub("_[A-H][0-9]+.*$", "", rownames(mapped)))
write.csv(metadata, file.path(out, "cell_metadata.csv"), row.names = FALSE)
# Values are already noninteger abundances. Do not re-normalize as UMI counts.
input <- list(VarFilter = "Unselect", Norm = "No_Normal")
for (scope in c("all_cells", "week2")) {
  subset <- if (scope == "week2") mapped[metadata$week == "2", , drop = FALSE] else mapped
  cat("Running", scope, "with", nrow(subset), "cells and", ncol(subset)-1, "genes\n")
  result <- SCMarkerfun(subset, GeneSK = 10, CellSK = 10, Labels = subset$Labels)
  stopifnot(nrow(result$ig) > 0, all(rownames(result$ig) %in% names(subset)))
  ranking <- data.frame(rank = seq_len(nrow(result$ig)), gene = rownames(result$ig),
                        score = result$ig[[1]])
  write.csv(ranking, file.path(out, paste0("scmarker_", scope, "_ranking.csv")), row.names = FALSE)
  saveRDS(list(input = subset, result = result, metadata = metadata[match(rownames(subset),metadata$cell),]),
          file.path(out, paste0("scmarker_", scope, ".rds")))
  print(head(ranking, 15))
}
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo.txt"))
