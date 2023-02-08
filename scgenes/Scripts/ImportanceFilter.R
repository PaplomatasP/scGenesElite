

#calculate the requested model and find the most importance genes according the running model .
#Return a expression matrix with the isolated genes .
ImportanceFilter = function(data,
                            Labels,
                            MLmethod, importanceLimit) {
  set.seed(1016)
  # Make Valid Column Names 
  colnames(data) <- make.names(colnames(data))
  
  MLlist = c(
    "rf",
    "knn",
    "svmRadial",
    "rpart",
    "C5.0",
    "glmnet",
    "lda",
    "bartMachine",
    "adaboost",
    "treebag",
    "xgbTree"
  )
  # Set ML parameters, based on the ML method chosen
  print(MLmethod)
  if (MLmethod %in% MLlist) {
    preProcMethod <- c("zv")
    
  } else
  {
    message("ML method input is not valid, exiting...")
    return(-1)
  }
  
 
  data$Labels = as.factor(Labels)
  
  partitionData <-
    caret::createDataPartition(data$Labels, p = 0.8, list = FALSE)
  trainData <- data[partitionData,]
  testData  <- data[-partitionData,]
  
  trainControl <-
    caret::trainControl(method = "repeatedcv",
                        number = 10,
                        repeats = 1,
                        p = 0.8,
                        preProcOptions=list(thresh=0.95,na.remove=TRUE,verbose=TRUE)
                        )
  

  
  
  model <-
    caret::train(
      Labels ~ .,
      data = trainData,
      method = MLmethod,  #MLmethod
      metric = "Accuracy",
      preProc = preProcMethod,
      trControl = trainControl,
      na.action = na.omit
    )
  
  df_imps = varImp(model)
  
  
  df_imps1 <- df_imps[["importance"]][1]
  df_imps1 <- subset(df_imps1, df_imps1[, 1] > input$importanceLimit)
  df_imps1 = data.frame(df_imps1[order(df_imps1, decreasing = TRUE), drop = FALSE,])
  newdata <- na.omit(df_imps1)
  
  
  
  iG <<- newdata
  
  newdata = data[, colnames(data) %in% rownames(newdata)]
  newdata <- cbind(as.data.frame(newdata), as.data.frame(Labels))
  newdata$Labels = as.factor(newdata$Labels)
  print("Finish")
  
  return(newdata)
}