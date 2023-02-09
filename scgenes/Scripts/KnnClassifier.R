
#Create a Knn model , and plot the confusion Matrix
KnnClassifier = function(data,
                            Labels
                            ) {
data$Labels = as.factor(Labels)

partitionData <-
  caret::createDataPartition(data$Labels, p = 0.8, list = FALSE)
trainData <- data[partitionData, ]
testData  <- data[-partitionData, ]


trainControl <-
  caret::trainControl(method = "cv",
                      number = 5,
                      p = 0.8)
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

ConfMatrix<-confusionMatrix(data=pred, testData$Labels, mode = "everything")
# 
# Confusion_Matrix<-as.data.frame(round(ConfMatrix[["byClass"]],2) )
# Confusion_Matrix$Metrics=rownames(Confusion_Matrix)
# colnames(Confusion_Matrix)[1]=""
# colnames(Confusion_Matrix)[2]=" k-NN Classifier"
# Confusion_Matrix=Confusion_Matrix[-c(1,2),]
Confusion_Matrix=ConfMatrixPlot(ConfMatrix)
return(Confusion_Matrix)

}


ConfMatrixPlot=function(ConfMatrix){
  
  cm_d <- as.data.frame(ConfMatrix$table)
  metrics <-as.data.frame(round(ConfMatrix[["byClass"]][-c(1,2)],2) )

  colnames(metrics)="ConfMatrix.metrics"
  metrics$confMatrix.metrics <- round(metrics[,1], digits = 2)
  metrics<-metrics[-1]
  colnames(metrics)="k-NN Metrics"
 

  cm_p <- as.data.frame(prop.table(ConfMatrix$table))
  cm_d$Perc <- round(cm_p$Freq*100,2)
  cm_st_p <-  gridExtra::tableGrob(metrics)
  
  
  cm_d_p <- ggplot2::ggplot(data = cm_d, aes(x = Prediction  , y =  Reference, fill = Freq))+
    ggplot2::geom_tile(plot.title = element_text(vjust = 12)) + 
   
    ggplot2::scale_fill_gradient(low = "#900700",
                                 high = "#0fbe0e",
                                 space = "Lab",
                                 na.value = "grey50",
                                 guide = "colourbar",
                                 aesthetics = "fill")+
    ggplot2::geom_text(aes(label = paste("",Freq,",",Perc,"%")), color = 'black', size = 7)+ 
    ggplot2::theme_bw() +
    ggplot2::guides(fill=FALSE) 
  gridExtra::grid.arrange(cm_d_p, cm_st_p,nrow = 1, ncol = 2, 
                          top= grid::textGrob("Confusion Matrix",x = 0.5, y = 0.6, just = "center", gp=grid::gpar(fontsize=19,font=1)))
  
}

