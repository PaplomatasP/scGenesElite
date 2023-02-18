#calculate the requested model and find the most importance genes according the running model .
#Return a expression matrix with the isolated genes .
ImportanceFilter = function(data,
                            Labels,
                            MLmethod, importanceLimit) {
  
  
  
  MLlist = c(
    "rf",
    "knn",
    "svmRadial",
    "rpart",
    "C5.0",
    "glmnet",
    "lda",
    "bartMachine",
    "treebag",
    "xgbTree"
  )
  # Set ML parameters, based on the ML method chosen
  print(MLmethod)
  if (MLmethod %in% MLlist) {
    preProcMethod <- c("knnImpute")
    
  } else
  {
    message("ML method input is not valid, exiting...")
    return(-1)
  }
  
  #| MLmethod=="rpart"
  if (MLmethod=="rpart" | MLmethod=="treebag" ) {  
    
    # Make Valid Column Names 
    colnames(data) <- make.names(colnames(data))
    
  }
  
  
  data$Labels = as.factor(Labels)
  
  partitionData <-
    caret::createDataPartition(data$Labels, p = 0.8, list = FALSE)
  trainData <- data[partitionData,]
  testData  <- data[-partitionData,]
  
  trainControl <-
    caret::trainControl(method = "repeatedcv",
                        number = 3,
                        repeats = 1,
                        p = 0.75
    )
  
  
  
  
  model<-
    caret::train(
      Labels ~ .,
      data = trainData,
      method = MLmethod,  #MLmethod
      metric = "Accuracy",
      preProc = preProcMethod,
      trControl = trainControl,
      na.action = na.omit
    )

  df_imps <- varImp(model)
  

   
  df_imps1 <- df_imps[["importance"]][1]
  df_imps1 <- subset(df_imps1, df_imps1[, 1] > input$importanceLimit)
  print("here2")
  
  
  df_imps1_mat <- as.matrix(df_imps1)
  df_imps1_mat <- df_imps1_mat[order(df_imps1_mat[,1], decreasing = TRUE),, drop = FALSE]
  df_imps1 <- as.data.frame(df_imps1_mat)
  
  newdata <- na.omit(df_imps1)

  
  
  
  iG <<- newdata
  
  newdata = data[, colnames(data) %in% rownames(newdata)]
  newdata <- cbind(as.data.frame(newdata), as.data.frame(Labels))
  newdata$Labels = as.factor(newdata$Labels)
  print("Finish")
  
  return(newdata)
}