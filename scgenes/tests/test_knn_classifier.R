library(shiny)
suppressMessages({
  library(caret); library(ggplot2); library(gridExtra)
  library(foreach); library(doParallel)
})

# KnnClassifier reads `input$genes` from its enclosing scope, exactly as it does
# when app.R sources it into the server frame.
run_knn <- function(NewData, iG, genes) {
  input <- list(genes = genes)
  source("Scripts/KnnClassifier.R", local = TRUE)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  KnnClassifier(data = NewData, iG, Labels = NewData[, ncol(NewData)])
}

# Five informative genes, ranked top, buried in 200 noise genes.
make_data <- function(label_column = "Labels") {
  set.seed(7)
  n <- 200L; informative <- 5L; noise <- 200L
  cls <- rep(c("CK", "CKp25"), each = n / 2)
  shift <- ifelse(cls == "CK", -3, 3)
  d <- as.data.frame(cbind(
    matrix(rnorm(n * informative, mean = rep(shift, informative)), nrow = n),
    matrix(rnorm(n * noise), nrow = n)
  ))
  colnames(d) <- paste0("Gene", seq_len(informative + noise))
  d[[label_column]] <- factor(cls)
  d
}

ranked <- function(d, label_column = "Labels") {
  genes <- setdiff(colnames(d), label_column)
  data.frame(score = rev(seq_along(genes)), row.names = genes)
}

d  <- make_data()
iG <- ranked(d)

res <- run_knn(d, iG, genes = 5)
stopifnot(
  "only the selected genes reach the model" =
    identical(attr(res, "genes"), paste0("Gene", 1:5)),
  "the label column is never used as a predictor" =
    !("Labels" %in% attr(res, "genes"))
)

cm <- attr(res, "confusionMatrix")
stopifnot(
  "the five informative genes must separate the classes" =
    unname(cm$overall[["Accuracy"]]) > 0.9
)

# The ensemble path names its label column after the expression that built it.
# That column must not be handed to the model as a feature.
d2 <- make_data(label_column = "as.factor.EnseLabels.")
res2 <- run_knn(d2, ranked(d2, "as.factor.EnseLabels."), genes = 5)
stopifnot(
  "a differently named label column stays out of the predictors" =
    identical(attr(res2, "genes"), paste0("Gene", 1:5))
)

# Regression for "task 1 failed": an earlier step can leave foreach registered
# against a cluster it already stopped. k-NN must not dispatch to it.
cl <- parallel::makeCluster(2)
doParallel::registerDoParallel(cl)
parallel::stopCluster(cl)
stopifnot("the test needs a stale registration" = foreach::getDoParRegistered())
res3 <- tryCatch(run_knn(d, iG, genes = 5),
                 error = function(e) stop("k-NN broke on a stale parallel backend: ",
                                          conditionMessage(e), call. = FALSE))
stopifnot(
  "the run still produces a confusion matrix" =
    !is.null(attr(res3, "confusionMatrix"))
)
foreach::registerDoSEQ()

cat("k-NN: trains on the selected genes only, keeps labels out of the features,",
    "and survives a stale parallel backend.\n")
