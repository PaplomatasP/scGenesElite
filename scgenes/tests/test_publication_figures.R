library(shiny)
library(ggplot2)
source("Scripts/KnnClassifier.R")
source("Scripts/PublicationFigures.R")
set.seed(42)
expr <- data.frame(A=rep(c(0,2,8,14),each=10),B=rep(3,40),C=seq_len(40),
 Labels=rep(c("Control","Treated"),each=20),row.names=paste0("cell",1:40))
r <- list(ig=data.frame(score=3:1,row.names=c("A","B","C")),newdata=expr,
 plot_context=list(expression=expr,input_genes=5L,mapped_genes=3L,method="Test"))
meta <- data.frame(cell=rev(rownames(expr)),sample=rev(rep(c("s1","s2","s3","s4"),each=10)))
d <- publication_data(r,3,"log2",meta,c("C","A"))
stopifnot(d$counts$n==c(5,3,3,3),d$constant=="B",all(d$heat$value[d$heat$gene=="B"]==0),
 identical(as.character(d$boxes$cell[1:40]),rownames(expr)),
 identical(as.character(d$boxes$group[1:40]),rep(c("s1","s2","s3","s4"),each=10)))
expected <- mean(scale(log2(expr$A+1))[1:10])
stopifnot(abs(d$heat$value[d$heat$gene=="A" & d$heat$group=="s1"]-expected)<1e-12)
must_fail <- function(code) stopifnot(inherits(tryCatch(force(code),error=identity),"error"))
must_fail(publication_data(r,3,metadata=meta[-1,]))
must_fail(publication_data(r,3,metadata=rbind(meta,meta[1,])))
must_fail(publication_data(r,3,inspect_genes="missing"))
negative <- r;negative$plot_context$expression$A[1] <- -2
must_fail(publication_data(negative,3))
stopifnot(publication_data(negative,3,"none")$boxes$value[1]==-2)
table3 <- matrix(c(7,2,1,0,8,2,1,1,8),3,
 dimnames=list(Prediction=c("a","b","c"),Reference=c("a","b","c")))
cm <- caret::confusionMatrix(as.table(table3),mode="everything")
p <- ConfMatrixPlot(cm)
stopifnot(all(abs(tapply(p$data$fraction,p$data$Reference,sum)-1)<1e-12))
saved_seed <- .Random.seed
expr$B <- seq_len(nrow(expr))/10
r$newdata <- expr
fit <- KnnClassifier(expr,r$ig,expr$Labels,genes_count=3,seed=33,fit_only=TRUE)
stopifnot(identical(saved_seed,.Random.seed),length(fit$test_cells)==8L,
 !length(intersect(fit$train_cells,fit$test_cells)),fit$genes==c("A","B","C"))
fit2 <- KnnClassifier(expr,r$ig,expr$Labels,genes_count=3,seed=33,fit_only=TRUE)
stopifnot(identical(fit$predictions,fit2$predictions))
real_knn <- KnnClassifier
testServer(function(input,output,session) {
 counter <- new.env();counter$calls <- 0L
 KnnClassifier <- function(...) {counter$calls <- counter$calls+1L;real_knn(...)}
 source("Scripts/PublicationFigures.R",local=TRUE)
 state <- reactiveValues(ready=TRUE)
 rr <- reactiveVal(r)
 api <- register_publication_results(input,output,session,rr,state)
}, {
 session$setInputs(genes=3,publicationSeed=33,publicationTransform="log2",
  publicationGrouping="labels",publicationInspect=c("A","C"),downloadDpi="300")
 initial <- api$classifier();before <- counter$calls
 session$setInputs(publicationTransform="none",publicationInspect="C")
 api$bundle();api$bundle();stopifnot(counter$calls==before)
 session$setInputs(generatePublication=1);stopifnot(!is.null(api$snapshot()))
 session$setInputs(publicationTransform="log2");stopifnot(is.null(api$snapshot()),counter$calls==before)
 session$setInputs(publicationSeed=34);api$classifier();stopifnot(counter$calls==before+1L)
 state$ready <- FALSE;session$flushReact();stopifnot(is.null(api$snapshot()))
})
cat("Publication: transforms, metadata, multiclass percentages, cached fits and invalidation passed.\n")
