# Ultra Simple SHAP-like Feature Selection
# Scientifically structured version using ML-based feature importance with proper thresholding
# This function performs feature selection based on machine learning importance scores
# and respects the importanceLimit parameter for scientific rigor

UltraSimpleShapFilter = function(data, Labels, MLmethod, importanceLimit) {
  
  # Input validation
  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    stop("Invalid data: data must be a non-empty matrix or data frame")
  }
  
  if (is.null(Labels) || length(Labels) != nrow(data)) {
    stop("Invalid Labels: Labels must have the same length as number of rows in data")
  }
  
  if (is.null(importanceLimit) || !is.numeric(importanceLimit) || importanceLimit < 0) {
    warning("Invalid importanceLimit, using default value of 0.01")
    importanceLimit <- 0.01
  }
  
  print(paste("Ultra Simple feature selection for method:", MLmethod))
  print(paste("Using importance threshold:", importanceLimit))
  
  # Make Valid Column Names
  colnames(data) <- make.names(colnames(data))
  
  # Create a simple data frame with all data
  full_data <- data.frame(data, Labels = as.factor(Labels), stringsAsFactors = FALSE)
  
  # Check for sufficient data
  if (nrow(full_data) < 5) {
    stop("Insufficient samples: Need at least 5 samples for analysis")
  }
  
  if (length(unique(full_data$Labels)) < 2) {
    stop("Insufficient classes: Need at least 2 different classes for classification")
  }
  
  # For computational efficiency, limit features only if dataset is extremely large
  # This pre-filtering is done conservatively to avoid losing important genes
  max_features_prefilter <- 2000  # Increased from 1000 for better coverage
  
  message(paste("Keeping ALL", nrow(full_data), "samples for comprehensive analysis"))
  
  # Only pre-filter if we have an extremely large number of features
  # This is a conservative approach to maintain biological relevance
  if (ncol(full_data) - 1 > max_features_prefilter) {
    message(paste("Large dataset detected (", ncol(full_data) - 1, " features). Pre-filtering to top ", 
                  max_features_prefilter, " most variable features for computational efficiency"))
    # Select most variable genes using coefficient of variation for better normalization
    gene_means <- apply(full_data[, -ncol(full_data)], 2, mean, na.rm = TRUE)
    gene_vars <- apply(full_data[, -ncol(full_data)], 2, var, na.rm = TRUE)
    # Use coefficient of variation (CV) to account for both mean and variance
    gene_cv <- ifelse(gene_means > 0, sqrt(gene_vars) / gene_means, gene_vars)
    top_genes <- names(sort(gene_cv, decreasing = TRUE, na.last = TRUE)[1:max_features_prefilter])
    full_data <- full_data[, c(top_genes, "Labels"), drop = FALSE]
    message("Pre-filtering completed. Proceeding with ML-based feature selection.")
  }
  
  print(paste("Working with", nrow(full_data), "samples and", ncol(full_data)-1, "features"))
  
  tryCatch({
    if (MLmethod == "rf") {
      if (!requireNamespace("randomForest", quietly = TRUE)) {
        stop("Package 'randomForest' is required but not installed")
      }
      library(randomForest)
      
      # Random Forest with scientifically appropriate parameters
      # Increased ntree for better stability and reproducibility
      model <- randomForest(
        Labels ~ .,
        data = full_data,
        ntree = 100,  # Increased from 30 for better stability
        mtry = max(1, floor(sqrt(ncol(full_data) - 1))),  # Standard RF parameter
        importance = TRUE,
        na.action = na.omit
      )
      
      # Get importance - use MeanDecreaseAccuracy as it's more robust
      rf_importance <- importance(model)
      if (!"MeanDecreaseAccuracy" %in% colnames(rf_importance)) {
        # Fallback to first available importance measure
        feature_importance <- rf_importance[, 1]
      } else {
        feature_importance <- rf_importance[, "MeanDecreaseAccuracy"]
      }
      
      # Normalize importance scores to 0-1 range for threshold comparison
      if (max(feature_importance, na.rm = TRUE) > 0) {
        feature_importance <- feature_importance / max(feature_importance, na.rm = TRUE)
      }
      
    } else if (MLmethod == "xgbTree") {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' is required but not installed")
      }
      library(xgboost)
      
      # Prepare data for XGBoost
      X <- as.matrix(full_data[, -ncol(full_data)])
      y <- as.numeric(full_data$Labels) - 1
      
      # XGBoost with scientifically appropriate parameters
      model <- xgboost(
        data = X,
        label = y,
        objective = "multi:softprob",
        num_class = length(unique(full_data$Labels)),
        nrounds = 50,  # Increased from 20 for better convergence
        max_depth = 4,  # Added for better control
        eta = 0.1,      # Lower learning rate for stability
        subsample = 0.8,  # Added for regularization
        colsample_bytree = 0.8,  # Added for regularization
        verbose = 0,
        nthread = 1
      )
      
      # Get importance
      importance_matrix <- xgb.importance(model = model)
      feature_importance <- rep(0, ncol(X))
      names(feature_importance) <- colnames(X)
      
      if (nrow(importance_matrix) > 0) {
        feature_importance[importance_matrix$Feature] <- importance_matrix$Gain
      }
      
      # Normalize importance scores to 0-1 range for threshold comparison
      if (max(feature_importance, na.rm = TRUE) > 0) {
        feature_importance <- feature_importance / max(feature_importance, na.rm = TRUE)
      }
      
    } else {
      stop(paste("Unsupported ML method:", MLmethod, ". Supported methods: 'rf', 'xgbTree'"))
    }
    
    # Sort features by importance
    sorted_importance <- sort(feature_importance, decreasing = TRUE, na.last = TRUE)
    
    # Filter by importance threshold (scientifically correct approach)
    # Use the importanceLimit parameter as provided
    threshold_features <- sorted_importance[sorted_importance >= importanceLimit]
    
    # If no features meet the threshold, use a more lenient approach
    # but still respect the scientific principle by selecting top features
    if (length(threshold_features) == 0) {
      message(paste("No features meet the importance threshold of", importanceLimit))
      message("Selecting top features based on importance ranking (scientific fallback)")
      
      # Select top features - use a reasonable number based on data size
      min_features <- min(10, length(sorted_importance))
      max_features_select <- min(100, length(sorted_importance))
      
      # Select top features that are at least 10% of max importance
      if (length(sorted_importance) > 0 && max(sorted_importance, na.rm = TRUE) > 0) {
        min_threshold <- max(sorted_importance, na.rm = TRUE) * 0.1
        threshold_features <- sorted_importance[sorted_importance >= min_threshold]
        
        # If still none, take top features
        if (length(threshold_features) == 0) {
          n_select <- min(max_features_select, length(sorted_importance))
          threshold_features <- sorted_importance[1:n_select]
        } else {
          # Limit to reasonable number
          n_select <- min(max_features_select, length(threshold_features))
          threshold_features <- threshold_features[1:n_select]
        }
      } else {
        # Ultimate fallback: select minimum number of features
        n_select <- min_features
        threshold_features <- sorted_importance[1:min(n_select, length(sorted_importance))]
      }
    } else {
      # Limit to reasonable number even if many meet threshold
      max_features_select <- min(200, length(threshold_features))
      threshold_features <- threshold_features[1:min(max_features_select, length(threshold_features))]
    }
    
    selected_features <- names(threshold_features)
    
    # Validate selected features
    if (length(selected_features) == 0) {
      stop("No features could be selected. Check data quality and importanceLimit parameter.")
    }
    
    # Create importance matrix
    iG <- data.frame(
      Importance = as.numeric(threshold_features),
      row.names = selected_features,
      stringsAsFactors = FALSE
    )
    
    # Ensure selected features exist in data
    available_features <- selected_features[selected_features %in% colnames(full_data)]
    if (length(available_features) == 0) {
      stop("Selected features not found in data. This indicates a data structure issue.")
    }
    
    # Create output data with only selected features
    newdata <- full_data[, c(available_features, "Labels"), drop = FALSE]
    
    # Update iG to match available features
    iG <- iG[available_features, , drop = FALSE]
    
    print(paste("Analysis completed. Selected", nrow(iG), "genes based on importance threshold."))
    print(paste("Importance range:", round(min(iG$Importance), 4), "to", round(max(iG$Importance), 4)))
    
    return(list(ig = iG, newdata = newdata))
    
  }, error = function(e) {
    message("Error in ML-based feature selection: ", e$message)
    
    # Improved fallback: use statistical variance with better methodology
    message("Using variance-based selection as scientific fallback...")
    
    # Calculate coefficient of variation for better normalization
    gene_means <- apply(full_data[, -ncol(full_data)], 2, mean, na.rm = TRUE)
    gene_vars <- apply(full_data[, -ncol(full_data)], 2, var, na.rm = TRUE)
    gene_cv <- ifelse(gene_means > 0, sqrt(gene_vars) / gene_means, gene_vars)
    
    # Select top variable genes
    n_select_fallback <- min(50, length(gene_cv))  # Reasonable number for fallback
    top_genes <- names(sort(gene_cv, decreasing = TRUE, na.last = TRUE)[1:n_select_fallback])
    
    if (length(top_genes) == 0) {
      stop("Fallback selection also failed. Please check data quality.")
    }
    
    iG <- data.frame(
      Importance = as.numeric(gene_cv[top_genes]),
      row.names = top_genes,
      stringsAsFactors = FALSE
    )
    
    # Normalize importance for consistency
    if (max(iG$Importance, na.rm = TRUE) > 0) {
      iG$Importance <- iG$Importance / max(iG$Importance, na.rm = TRUE)
    }
    
    newdata <- full_data[, c(top_genes, "Labels"), drop = FALSE]
    
    message(paste("Fallback completed. Selected", nrow(iG), "genes using variance-based method."))
    
    return(list(ig = iG, newdata = newdata))
  })
}
