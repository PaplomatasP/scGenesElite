# Verify committed numerical outputs without modifying the case-study package.
out <- "case-studies/ck-p25"
d <- read.csv(file.path(out, "input_ExampleData.csv"), row.names = 1, check.names = FALSE)
stopifnot(identical(dim(d), c(384L, 2001L)))
metadata <- read.csv(file.path(out, "cell_metadata.csv"))
stopifnot(identical(metadata$cell, rownames(d)), all(table(metadata$sample) == 96L))
sets <- strsplit(readLines(file.path(out, "KEGG_2019_Mouse.gmt"), warn = FALSE), "\t", fixed = TRUE)
sets <- Filter(function(z) length(z) >= 3L && nzchar(trimws(z[[1]])), sets)
terms <- vapply(sets, `[[`, character(1), 1)
sets <- setNames(lapply(sets, function(z) unique(toupper(z[-c(1, 2)]))), terms)
stopifnot(length(sets) == 303L)
ranks <- list()
for (scope in c("all_cells", "week2")) {
  run <- readRDS(file.path(out, paste0("scmarker_", scope, ".rds")))
  rank <- read.csv(file.path(out, paste0("scmarker_", scope, "_ranking.csv")))
  expected <- if (scope == "all_cells") c(384L, 302L) else c(192L, 374L)
  stopifnot(nrow(run$input) == expected[1], nrow(rank) == expected[2], ncol(run$input) == 1861L)
  stopifnot(identical(rank$gene, rownames(run$result$ig)))
  stopifnot(isTRUE(all.equal(rank$score, run$result$ig[[1]])))
  genes <- names(run$input)[-ncol(run$input)]
  stopifnot(isTRUE(all.equal(as.matrix(run$input[, genes]), as.matrix(d[rownames(run$input), genes]), check.attributes = FALSE)))
  universe <- unique(toupper(genes))
  for (cutoff in c(20, 50, 100)) {
    selected <- toupper(head(rank$gene, cutoff))
    p <- vapply(sets, function(members) {
      members <- intersect(members, universe)
      phyper(length(intersect(selected, members)) - 1, length(members), length(universe) - length(members), length(selected), lower.tail = FALSE)
    }, numeric(1))
    fdr <- p.adjust(p, "BH")
    saved <- read.csv(file.path(out, paste0("kegg_", scope, "_top", cutoff, ".csv")))
    stopifnot(nrow(saved) == 303L, !any(saved$fdr < 0.05))
    stopifnot(isTRUE(all.equal(saved$p_value, unname(p[saved$term]), tolerance = 1e-12)))
    stopifnot(isTRUE(all.equal(saved$fdr, unname(fdr[saved$term]), tolerance = 1e-12)))
  }
  ranks[[scope]] <- rank
}
stopifnot(length(intersect(head(ranks$all_cells$gene, 50), head(ranks$week2$gene, 50))) == 7L)
run <- readRDS(file.path(out, "scmarker_all_cells.rds"))
genes <- head(ranks$all_cells$gene, 20)
scaled <- scale(log2(as.matrix(run$input[, genes]) + 1))
heat <- read.csv(file.path(out, "figure_expression_values.csv"))
computed <- mapply(function(g, s) mean(scaled[run$metadata$sample == s, g]), heat$gene, heat$sample)
stopifnot(isTRUE(all.equal(heat$value, unname(computed), tolerance = 1e-12)))
message("Verified both rankings, 768,000 input values, top-50 overlap, heatmap values and all six KEGG analyses.")
