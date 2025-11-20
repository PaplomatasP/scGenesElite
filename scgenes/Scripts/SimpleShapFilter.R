# Simple SHAP-like Feature Selection (Fallback)
# A simplified version that uses Random Forest importance when SHAP fails

SimpleShapFilter = function(data, Labels, MLmethod, importanceLimit) {
  
  print(paste("Simple SHAP-like calculation for method:", MLmethod))
  
  # Make Valid Column Names
  colnames(data) <- make.names(colnames(data))
  
  # Prepare data
  data$Labels = as.factor(Labels)
  
  # Check data size
  if (nrow(data) < 5) {
    message("Insufficient data for analysis.")
    return(NULL)
  }
  
  tryCatch({
    if (MLmethod == "rf") {
      library(randomForest)
      
      # Train Random Forest model
      model <- randomForest(
        Labels ~ .,
        data = data,
        ntree = 50,  # Reduced number of trees
        importance = TRUE
      )
      
      # Get variable importance
      rf_importance <- importance(model)
      mean_importance <- rf_importance[, "MeanDecreaseAccuracy"]
      feature_names <- names(mean_importance)
      
    } else if (MLmethod == "xgbTree") {
      library(xgboost)
      
      # Prepare data for XGBoost
      train_matrix <- xgb.DMatrix(
        data = as.matrix(data[, -ncol(data)]),
        label = as.numeric(data$Labels) - 1
      )
      
      # Train XGBoost model with reduced complexity
      params <- list(
        objective = "multi:softprob",
        num_class = length(unique(data$Labels)),
        max_depth = 3,
        eta = 0.3,
        nthread = 1
      )
      
      model <- xgb.train(
        params = params,
        data = train_matrix,
        nrounds = 20,  # Reduced rounds
        verbose = 0
      )
      
      # Get feature importance
      importance_matrix <- xgb.importance(model = model)
      mean_importance <- rep(0, ncol(data) - 1)
      names(mean_importance) <- colnames(data[, -ncol(data)])
      mean_importance[importance_matrix$Feature] <- importance_matrix$Gain
      feature_names <- names(mean_importance)
    }
    
    # Create importance dataframe
    importance_df <- data.frame(
      Feature = feature_names,
      Importance = mean_importance,
      stringsAsFactors = FALSE
    )
    
    # Sort by importance
    importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
    
    # Filter by importance threshold
    filtered_importance <- importance_df[importance_df$Importance > importanceLimit, ]
    
    if (nrow(filtered_importance) == 0) {
      message("No features meet the importance threshold. Using top 10 features.")
      filtered_importance <- importance_df[1:min(10, nrow(importance_df)), ]
    }
    
    # Create gene importance matrix
    iG <- data.frame(
      Importance = filtered_importance$Importance,
      row.names = filtered_importance$Feature
    )
    
    # Select genes from original data
    selected_genes <- rownames(iG)
    available_genes <- colnames(data)[colnames(data) %in% selected_genes]
    
    if (length(available_genes) == 0) {
      message("No genes found in the dataset matching the selected features.")
      return(NULL)
    }
    
    # Filter data to include only selected genes
    newdata <- data[, available_genes, drop = FALSE]
    newdata$Labels <- as.factor(Labels)
    
    # Update iG to only include available genes
    iG <- iG[available_genes, , drop = FALSE]
    
    print(paste("Simple SHAP-like analysis completed. Selected", nrow(iG), "genes."))
    
    return(list(ig = iG, newdata = newdata))
    
  }, error = function(e) {
    message("Error in Simple SHAP calculation: ", e$message)
    return(NULL)
  })
}
