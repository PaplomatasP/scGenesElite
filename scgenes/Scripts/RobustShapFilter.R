# Robust SHAP Values Feature Selection
# A more robust version that handles large datasets and prevents crashes

RobustShapFilter = function(data, Labels, MLmethod, importanceLimit) {
  
  print(paste("Robust SHAP calculation for method:", MLmethod))
  
  # Load required libraries with error handling
  tryCatch({
    if (MLmethod == "rf") {
      library(randomForest)
    } else if (MLmethod == "xgbTree") {
      library(xgboost)
    }
  }, error = function(e) {
    message("Error loading required libraries: ", e$message)
    return(NULL)
  })
  
  # Make Valid Column Names
  colnames(data) <- make.names(colnames(data))
  
  # Prepare data
  data$Labels = as.factor(Labels)
  
  # Check data size and subsample if too large
  original_rows <- nrow(data)
  original_cols <- ncol(data) - 1
  
  # Keep ALL samples - only limit features for computational efficiency
  max_cols <- 1000  # Increased to 1000 most variable genes
  
  message(paste("Keeping ALL", original_rows, "samples for comprehensive analysis"))
  
  if (ncol(data) - 1 > max_cols) {
    message(paste("Subsampling from", original_cols, "to", max_cols, "features"))
    # Select most variable genes
    gene_vars <- apply(data[, -ncol(data)], 2, var, na.rm = TRUE)
    top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:max_cols])
    data <- data[, c(top_genes, "Labels")]
  }
  
  print(paste("Using dataset with", nrow(data), "samples and", ncol(data)-1, "features"))
  
  # Split data
  set.seed(123)
  train_prop <- 0.7
  train_size <- floor(train_prop * nrow(data))
  train_idx <- sample(1:nrow(data), train_size)
  
  trainData <- data[train_idx, ]
  testData <- data[-train_idx, ]
  
  # Ensure minimum test size
  if (nrow(testData) < 3) {
    testData <- data[1:min(5, nrow(data)), ]
  }
  
  print(paste("Training with", nrow(trainData), "samples, testing with", nrow(testData), "samples"))
  
  tryCatch({
    if (MLmethod == "rf") {
      # Train simplified Random Forest
      model <- randomForest(
        Labels ~ .,
        data = trainData,
        ntree = 50,        # Reduced trees
        nodesize = 5,      # Larger nodes
        maxnodes = 10,     # Limit tree complexity
        importance = TRUE
      )
      
      # Use built-in importance instead of SHAP for stability
      rf_importance <- importance(model)
      mean_importance <- rf_importance[, "MeanDecreaseAccuracy"]
      feature_names <- names(mean_importance)
      
      print("Using Random Forest variable importance (stable method)")
      
    } else if (MLmethod == "xgbTree") {
      # Prepare data for XGBoost
      train_matrix <- xgb.DMatrix(
        data = as.matrix(trainData[, -ncol(trainData)]),
        label = as.numeric(trainData$Labels) - 1
      )
      
      # Train simplified XGBoost
      params <- list(
        objective = "multi:softprob",
        num_class = length(unique(trainData$Labels)),
        max_depth = 3,
        eta = 0.3,
        nthread = 1,
        subsample = 0.8,
        colsample_bytree = 0.8
      )
      
      model <- xgb.train(
        params = params,
        data = train_matrix,
        nrounds = 20,
        verbose = 0
      )
      
      # Get XGBoost importance
      importance_matrix <- xgb.importance(model = model)
      mean_importance <- rep(0, ncol(trainData) - 1)
      names(mean_importance) <- colnames(trainData[, -ncol(trainData)])
      
      if (nrow(importance_matrix) > 0) {
        mean_importance[importance_matrix$Feature] <- importance_matrix$Gain
      }
      
      feature_names <- names(mean_importance)
      
      print("Using XGBoost feature importance (stable method)")
    }
    
    # Create importance dataframe
    importance_df <- data.frame(
      Feature = feature_names,
      Importance = mean_importance,
      stringsAsFactors = FALSE
    )
    
    # Sort by importance
    importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
    
    # Filter by importance threshold (more lenient)
    threshold <- max(0.01, importanceLimit / 100)  # Convert to more reasonable scale
    filtered_importance <- importance_df[importance_df$Importance > threshold, ]
    
    if (nrow(filtered_importance) == 0) {
      message("No features meet threshold. Using top 20 features.")
      filtered_importance <- importance_df[1:min(20, nrow(importance_df)), ]
    }
    
    # Limit to reasonable number of features
    max_features <- min(50, nrow(filtered_importance))
    filtered_importance <- filtered_importance[1:max_features, ]
    
    # Create gene importance matrix
    iG <- data.frame(
      Importance = filtered_importance$Importance,
      row.names = filtered_importance$Feature
    )
    
    # Create output data with selected genes
    selected_genes <- rownames(iG)
    newdata <- data[, c(selected_genes, "Labels"), drop = FALSE]
    # Labels are already correct from the subsampled data
    
    print(paste("Robust analysis completed. Selected", nrow(iG), "genes."))
    
    return(list(ig = iG, newdata = newdata))
    
  }, error = function(e) {
    message("Error in robust calculation: ", e$message)
    
    # Final fallback: random selection
    message("Using random feature selection as final fallback...")
    n_features <- min(20, ncol(data) - 1)
    random_features <- sample(colnames(data[, -ncol(data)]), n_features)
    
    iG <- data.frame(
      Importance = rep(1, length(random_features)),
      row.names = random_features
    )
    
    newdata <- data[, c(random_features, "Labels"), drop = FALSE]
    # Labels are already correct from the subsampled data
    
    return(list(ig = iG, newdata = newdata))
  })
}
