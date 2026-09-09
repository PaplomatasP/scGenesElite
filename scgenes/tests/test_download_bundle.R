library(shiny)
source("Scripts/KnnClassifier.R")
source("Scripts/PublicationFigures.R")
d <- data.frame(A=1:40,B=seq(2,80,2),Labels=rep(c("a","b"),each=20))
r <- list(ig=data.frame(score=2:1,row.names=c("A","B")),newdata=d,
 plot_context=list(input_genes=4,mapped_genes=2,expression=d,method="Test"))
fit <- KnnClassifier(d,r$ig,d$Labels,genes_count=2,seed=9,fit_only=TRUE)
dat <- publication_data(r,2)
b <- list(data=dat,plots=publication_plots(dat,fit),classifier=fit,seed=9,
 filtered=d,ranking=data.frame(gene=c("A","B"),score=2:1))
archive <- tempfile(fileext=".zip")
write_publication_bundle(archive,b,300)
entries <- zip::zip_list(archive)$filename
expected <- c("composite.pdf","composite.png","selection.pdf","selection.png",
 "heatmap.pdf","heatmap.png","expression.pdf","expression.png","classification.pdf",
 "classification.png","heat.csv","boxes.csv","counts.csv","groups.csv","ranking.csv",
 "FilterData.csv","predictions.csv","confusion.csv","classification_metrics.txt","classification.rds","README.txt")
stopifnot(setequal(entries,expected))
directory <- tempfile();dir.create(directory);utils::unzip(archive,exdir=directory)
stopifnot(identical(readRDS(file.path(directory,"classification.rds"))$predictions,fit$predictions))
path <- tempfile(fileext=".png")
export_publication_plot(path,b$plots,"selection","png",600,2)
size300 <- dim(png::readPNG(file.path(directory,"selection.png")))[1:2]
size600 <- dim(png::readPNG(path))[1:2]
stopifnot(all(abs(size600/size300-2)<.01))
unlink(c(directory,archive,path),recursive=TRUE)
cat("Publication downloads: plots, cached predictions, tables, methods and 300/600 dpi passed.\n")
