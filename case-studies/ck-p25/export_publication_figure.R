# Regenerate the paper figure through the same helpers used by the live app.
library(shiny)
source('Scripts/KnnClassifier.R')
source('Scripts/PublicationFigures.R')
run <- readRDS('case-studies/ck-p25/scmarker_all_cells.rds')
result <- run$result
result$plot_context <- list(expression=run$input,input_genes=2000,
  mapped_genes=ncol(run$input)-1L,method='SCMarker')
metadata <- run$metadata
metadata$sample <- c('CK 0 weeks','CK-p25 0 weeks','CK 2 weeks','CK-p25 2 weeks')[
  match(metadata$sample,c('CK_0w_m1','CKp25_0w_m1','CK_2w_m2','CKp25_2w_m2'))]
metadata <- metadata[order(match(metadata$sample,c('CK 0 weeks','CK-p25 0 weeks','CK 2 weeks','CK-p25 2 weeks'))),]
write.csv(metadata[,c('cell','sample')],'case-studies/ck-p25/publication_metadata.csv',row.names=FALSE)
d <- publication_data(result,20,'log2',metadata,c('Il6ra','Sall1','Lcp1'))
fit <- KnnClassifier(result$newdata,result$ig,result$newdata[[ncol(result$newdata)]],
  genes_count=20,seed=20260908,fit_only=TRUE)
old <- attr(readRDS('case-studies/ck-p25/collage_knn.rds'),'confusionMatrix')
stopifnot(identical(as.vector(fit$confusionMatrix$table),as.vector(old$table)))
expected <- read.csv('case-studies/ck-p25/figure_expression_values.csv')
stopifnot(max(abs(d$heat$value-expected$value))<1e-12)
plots <- publication_plots(d,fit)
bundle <- list(data=d,plots=plots,classifier=fit,seed=20260908,
  filtered=result$newdata,ranking=data.frame(gene=rownames(result$ig),score=result$ig[[1]]))
export_publication_plot('case-studies/ck-p25/figures/scGenesFinder_results_collage.pdf',plots,n=20)
export_publication_plot('case-studies/ck-p25/figures/scGenesFinder_results_collage_600dpi.png',plots,format='png',n=20)
write_publication_bundle('case-studies/ck-p25/scGenesFinder_app_publication_bundle.zip',bundle,600)
writeLines(capture.output(sessionInfo()), 'case-studies/ck-p25/publication_sessionInfo.txt')
cat('App export exactly reproduces the saved heatmap values and 48/76 correct predictions.\n')
