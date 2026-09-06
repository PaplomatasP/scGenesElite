
#Create a Knn model , and plot the confusion Matrix
KnnClassifier = function(data,iG,
                            Labels
                            ) {

  Labels <- droplevels(as.factor(Labels))
  genes <- rownames(head(iG, input$genes))

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

Confusion_Matrix=ConfMatrixPlot(ConfMatrix)
# Carried along so the run can be inspected without re-reading the plot.
attr(Confusion_Matrix, "confusionMatrix") <- ConfMatrix
attr(Confusion_Matrix, "genes") <- selected
return(Confusion_Matrix)

}


ConfMatrixPlot=function(ConfMatrix){
  
  cm_d <- as.data.frame(ConfMatrix$table)

  byClass <- ConfMatrix[["byClass"]]
  if (is.matrix(byClass)) {
    # Three or more classes: caret returns one row per class, so keep the
    # classes as columns instead of indexing the matrix as a flat vector.
    metrics <- as.data.frame(round(t(byClass[, -c(1, 2), drop = FALSE]), 2))
  } else {
    metrics <- data.frame(round(byClass[-c(1, 2)], 2))
    colnames(metrics) <- "k-NN Metrics"
  }

  cm_p <- as.data.frame(prop.table(ConfMatrix$table))
  cm_d$Perc <- round(cm_p$Freq*100,2)
  cm_st_p <-  gridExtra::tableGrob(metrics)

  # Colour by whether the cell sits on the diagonal. The old gradient ran on the
  # count alone, which painted a large pile of misclassifications green and a
  # small number of correct calls dark red.
  cm_d$Outcome <- factor(
    ifelse(cm_d$Prediction == cm_d$Reference, "Correct", "Misclassified"),
    levels = c("Correct", "Misclassified")
  )

  cm_d_p <- ggplot2::ggplot(data = cm_d, aes(x = Prediction  , y =  Reference))+
    ggplot2::geom_tile(aes(fill = Outcome, alpha = Perc), colour = "white") +
    ggplot2::scale_fill_manual(values = c(Correct = "#0fbe0e",
                                          Misclassified = "#900700")) +
    ggplot2::scale_alpha_continuous(range = c(0.15, 1), limits = c(0, 100)) +
    # Put the diagonal top-left, the way the printed caret table reads.
    ggplot2::scale_y_discrete(limits = rev(levels(cm_d$Reference))) +
    ggplot2::geom_text(aes(label = paste(Freq, "\n", Perc,"%")), color = 'black', size = 7)+  # Adjusted size and added newline
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = element_text(vjust = 12)) +
    ggplot2::guides(fill="none", alpha="none")
  gridExtra::grid.arrange(cm_d_p, cm_st_p,nrow = 1, ncol = 2, 
                          top= grid::textGrob("Confusion Matrix",x = 0.5, y = 0.6, just = "center", gp=grid::gpar(fontsize=21,font=1)))
  
}
