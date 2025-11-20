#calculate the requested model and find the most importance genes according the running model .
#Return a expression matrix with the isolated genes .
ImportanceFilter = function(data,
                            Labels,
                            MLmethod, importanceLimit) {
  
  
  MLlist = c(
    "rf",
    "rpart",
    "C5.0",
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
  
  if (input$GENEid=="EnsemblGenes"){
    colnames(data) <- make.names(colnames(data))
  }
  
  
  
  data$Labels = as.factor(Labels)
  
  
  trainControl <-
    caret::trainControl(method = "repeatedcv",
                        number = 3,
                        repeats = 1
    )
  
  print("before train")
  
  if (MLmethod == "xgbTree") {
    # Only include hyperparameters that are being varied
    tune_grid <- expand.grid(
      nrounds = c(50),               
      max_depth = c(2),
      eta = c(0.1),                   # Constant value
      gamma = c(0),                   # Constant value
      colsample_bytree = c(1),        # Constant value
      min_child_weight = c(1),        # Constant value
      subsample = c(1)                # Constant value
    )
    print("iam xgbTree")
    
    model <- caret::train(
      Labels ~ .,
      data = data,
      method = "xgbTree",
      tuneGrid = tune_grid,
      metric = "Accuracy",
      trControl = trainControl(verboseIter = FALSE) # Control verbosity
    )
    
  } else if (MLmethod == "rf") {
    
    # Define the tuning grid
    tune_grid <- expand.grid(mtry = c(3, 6))
    print("iam rf")
    # Train the model 
    model <- caret::train(
      Labels ~ .,
      data = data,
      method = "rf",
      tuneGrid = tune_grid,
      ntree = 100,
      nodesize = 10,  # Minimum size of terminal nodes
      trace = FALSE   # Suppress messages
    )
    
  }else if (MLmethod == "C5.0") {
    print("iam C5.0")
    # Set up the parameter grid
    tune_grid <- expand.grid(model = "tree", # tree or rules
                             trials = seq(1, 20, 2), # number of boosting iterations
                             winnow = FALSE ) # whether to use winnowing
    
    model <- caret::train(
      Labels ~ .,
      data = data,
      method = "C5.0",
      tuneGrid = tune_grid
    )
  } else {
    print("iam without")
    model <- caret::train(
      Labels ~ .,
      data = data,
      method = MLmethod
    )
  }
  
  print("after train")
  
  df_imps <- varImp(model)
  
  print("after after")
  
   
  df_imps1 <- df_imps[["importance"]][1]
  df_imps1 <- subset(df_imps1, df_imps1[, 1] > input$importanceLimit)
  print("here2")
  
  
  df_imps1_mat <- as.matrix(df_imps1)
  df_imps1_mat <- df_imps1_mat[order(df_imps1_mat[,1], decreasing = TRUE),, drop = FALSE]
  df_imps1 <- as.data.frame(df_imps1_mat)
  
  newdata <- na.omit(df_imps1)

  
  
  
  iG <- newdata
  
  newdata = data[, colnames(data) %in% rownames(newdata)]
  newdata <- cbind(as.data.frame(newdata), as.data.frame(Labels))
  newdata$Labels = as.factor(newdata$Labels)
  print("Finish")
  
  return(list(ig=iG, newdata=newdata))
}