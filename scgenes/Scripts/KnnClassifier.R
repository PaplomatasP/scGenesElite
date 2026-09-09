
#Create a Knn model , and plot the confusion Matrix
KnnClassifier = function(data,iG,
                            Labels, genes_count = input$genes, seed = NULL,
                            fit_only = FALSE
                            ) {

  Labels <- droplevels(as.factor(Labels))
  genes <- rownames(head(iG, genes_count))
  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
            else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
              rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }

  # `rownames(by_iG) %in% colnames(data)` is as long as the gene list, not as
  # long as the frame, so using it to index columns recycled it: when every
  # selected gene was present it was all TRUE and kept *every* column, and the
  # model was trained on the full matrix instead of the selected genes.
  selected <- intersect(genes, colnames(data))

  validate(
    need(length(selected) >= 2,
         "k-NN needs at least two of the selected genes to be present in the data."),
    need(nlevels(Labels) >= 2,
         "k-NN needs at least two classes in the label column."),
    need(min(table(Labels)) >= 2,
         "k-NN needs at least two cells per class to build a train/test split.")
  )

  # Rebuild the frame from the selected genes only. A label column travelling in
  # `data` under another name (the ensemble path names it after the expression
  # that built it) would otherwise be handed to the model as a predictor.
  data <- data[, selected, drop = FALSE]
  data$Labels <- Labels

partitionData <-
  caret::createDataPartition(data$Labels, p = 0.8, list = FALSE)
trainData <- data[partitionData, , drop = FALSE]
testData  <- data[-partitionData, , drop = FALSE]
shiny::validate(shiny::need(nrow(testData) > 0L &&
  all(table(factor(testData$Labels, levels=levels(Labels))) > 0L),
  "Each class needs enough cells for a non-empty 20% test split."))

# Never ask for more folds than the smallest class can fill.
folds <- max(2L, min(5L, min(table(trainData$Labels))))

# allowParallel = FALSE: resampling a k-NN on a handful of genes is cheap, and
# it keeps this panel from dispatching to whatever foreach backend an earlier
# step happened to leave registered.
trainControl <-
  caret::trainControl(method = "cv",
                      number = folds,
                      p = 0.8,
                      allowParallel = FALSE)
preProcMethod <- c("center", "scale")
model <-
  caret::train(
    Labels ~ .,
    data = trainData,
    method = "knn",
    metric = "Accuracy",
    preProc = preProcMethod,
    trControl = trainControl,
    na.action = na.omit
  )


pred = predict(model, newdata=testData)

ConfMatrix<-caret::confusionMatrix(data=pred, testData$Labels, mode = "everything")

if (fit_only) return(list(confusionMatrix = ConfMatrix, genes = selected,
  seed = seed, k = model$bestTune$k, folds = folds,
  train_cells = rownames(trainData), test_cells = rownames(testData),
  predictions = data.frame(cell=rownames(testData), observed=testData$Labels,
                           predicted=pred)))

Confusion_Matrix=ConfMatrixPlot(ConfMatrix)
# Carried along so the run can be inspected without re-reading the plot.
attr(Confusion_Matrix, "confusionMatrix") <- ConfMatrix
attr(Confusion_Matrix, "genes") <- selected
return(Confusion_Matrix)

}

ConfMatrixPlot <- function(ConfMatrix, title = "D  Within-dataset classification", subtitle = NULL) {
  d <- as.data.frame(ConfMatrix$table)
  totals <- tapply(d$Freq, d$Reference, sum)
  d$fraction <- d$Freq / totals[as.character(d$Reference)]
  d$label <- sprintf("%d cells\n%.1f%%", d$Freq, 100*d$fraction)
  if (is.null(subtitle)) subtitle <- sprintf("k-NN | accuracy %.1f%%", 100*ConfMatrix$overall['Accuracy'])
  ggplot2::ggplot(d, ggplot2::aes(Prediction, Reference, fill=fraction)) +
    ggplot2::geom_tile(colour="white", linewidth=1) +
    ggplot2::geom_text(ggplot2::aes(label=label), size=3, colour="#18243B") +
    ggplot2::scale_fill_gradient(low="#F2F5FA", high="#8DA6D4", limits=c(0,1), guide="none") +
    ggplot2::scale_x_discrete(expand=c(0,0)) +
    ggplot2::scale_y_discrete(limits=rev(levels(d$Reference)), expand=c(0,0)) +
    ggplot2::labs(title=title, subtitle=subtitle, x="Predicted label", y="Observed label",
      caption=paste0(sum(d$Freq), " held-out cells; percentages within observed class\n",
                     "Cell-level split; independent samples were not held out")) +
    ggplot2::theme_minimal(base_size=8) +
    ggplot2::theme(text=ggplot2::element_text(colour="#18243B"),
      panel.grid=ggplot2::element_blank(),plot.title=ggplot2::element_text(size=10,face="bold"),
      plot.subtitle=ggplot2::element_text(size=7,colour="#526178",margin=ggplot2::margin(b=8)),
      plot.caption=ggplot2::element_text(size=7,hjust=0),plot.margin=ggplot2::margin(8,12,8,8))
}
