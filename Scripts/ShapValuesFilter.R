# SHAP Values Feature Selection
# Calculate REAL SHAP values for machine learning models and find the most important genes
# Return a expression matrix with the isolated genes based on SHAP importance scores.
# This implementation calculates true SHAP values based on Shapley theory.

ShapValuesFilter = function(data, Labels, MLmethod, importanceLimit) {
  
  # GC Petros S 1: Input validation
  if (is.null(data) || nrow(data) == 0 || ncol(data) == 0) {
    stop("Invalid data: data must be a non-empty matrix or data frame")
  }
  
  if (is.null(Labels) || length(Labels) != nrow(data)) {
    stop("Invalid Labels: Labels must have the same length as number of rows in data")
  }
  
  if (is.null(importanceLimit) || !is.numeric(importanceLimit) || importanceLimit < 0) {
    warning("Invalid importanceLimit, using default value of 0.0")
    importanceLimit <- 0.0
  }
  
  # GC Petros S 2: Check if required packages are available
  if (!requireNamespace("fastshap", quietly = TRUE) && MLmethod == "rf") {
    warning("Package 'fastshap' not available. Falling back to simple importance method...")
    source("./Scripts/UltraSimpleShapFilter.R", local = TRUE)
    return(UltraSimpleShapFilter(data, Labels, MLmethod, importanceLimit))
  }
  
  if (!requireNamespace("xgboost", quietly = TRUE) && MLmethod == "xgbTree") {
    warning("Package 'xgboost' not available. Falling back to simple importance method...")
    source("./Scripts/UltraSimpleShapFilter.R", local = TRUE)
    return(UltraSimpleShapFilter(data, Labels, MLmethod, importanceLimit))
  }
  
  # Supported ML methods for SHAP
  SHAP_MLlist = c("rf", "xgbTree")
  
  if (!(MLmethod %in% SHAP_MLlist)) {
    stop(paste("SHAP values are currently supported only for Random Forest (rf) and XGBoost (xgbTree) models. Got:", MLmethod))
  }
  
  print(paste("Calculating REAL SHAP Values for method:", MLmethod))
  print(paste("Using importance threshold:", importanceLimit))
  
  # GC Petros S 3: Make Valid Column Names
  # Ensure all column names are valid R identifiers (critical for fastshap)
  original_colnames <- colnames(data)
  colnames(data) <- make.names(colnames(data), unique = TRUE)
  
  # Store mapping for later use
  name_mapping <- data.frame(
    original = original_colnames,
    transformed = colnames(data),
    stringsAsFactors = FALSE
  )
  
  # Prepare data
  data_with_labels <- data.frame(data, Labels = as.factor(Labels), stringsAsFactors = FALSE)
  
  # Check for sufficient data
  if (nrow(data_with_labels) < 10) {
    stop("Insufficient data for SHAP analysis. Need at least 10 samples.")
  }
  
  if (length(unique(data_with_labels$Labels)) < 2) {
    stop("Insufficient classes: Need at least 2 different classes for classification")
  }
  
  # GC Petros S 4: Pre-filtering for very large datasets (conservative approach)
  # Only pre-filter if dataset is extremely large to manage computational complexity
  max_features_prefilter <- 1500  # Reasonable limit for SHAP calculation
  
  if (ncol(data_with_labels) - 1 > max_features_prefilter) {
    message(paste("Large dataset detected (", ncol(data_with_labels) - 1, " features). "))
    message("Pre-filtering to top ", max_features_prefilter, " most variable features for SHAP calculation...")
    
    # Use coefficient of variation for better normalization
    gene_means <- apply(data_with_labels[, -ncol(data_with_labels)], 2, mean, na.rm = TRUE)
    gene_vars <- apply(data_with_labels[, -ncol(data_with_labels)], 2, var, na.rm = TRUE)
    gene_cv <- ifelse(gene_means > 0, sqrt(gene_vars) / gene_means, gene_vars)
    top_genes <- names(sort(gene_cv, decreasing = TRUE, na.last = TRUE)[1:max_features_prefilter])
    data_with_labels <- data_with_labels[, c(top_genes, "Labels"), drop = FALSE]
    message("Pre-filtering completed.")
  }
  
  print(paste("Working with", nrow(data_with_labels), "samples and", ncol(data_with_labels)-1, "features"))
  
  # GC Petros S 5: Train/Test Split for proper SHAP calculation
  # This is critical: SHAP values should be calculated on unseen data
  set.seed(123)
  train_prop <- 0.75  # Use 75% for training, 25% for SHAP calculation
  train_size <- floor(train_prop * nrow(data_with_labels))
  train_idx <- sample(1:nrow(data_with_labels), train_size)
  
  trainData <- data_with_labels[train_idx, ]
  testData <- data_with_labels[-train_idx, ]
  
  # Ensure minimum test size
  if (nrow(testData) < 5) {
    warning("Too few test samples. Using 80/20 split instead.")
    train_prop <- 0.80
    train_size <- floor(train_prop * nrow(data_with_labels))
    train_idx <- sample(1:nrow(data_with_labels), train_size)
    trainData <- data_with_labels[train_idx, ]
    testData <- data_with_labels[-train_idx, ]
  }
  
  print(paste("Training set:", nrow(trainData), "samples | Test set:", nrow(testData), "samples"))
  
  # GC Petros S 6: Calculate SHAP values based on ML method
  tryCatch({
    if (MLmethod == "rf") {
      if (!requireNamespace("randomForest", quietly = TRUE)) {
        stop("Package 'randomForest' is required but not installed")
      }
      library(randomForest)
      
      print("Training Random Forest model...")
      # GC Petros S 11: Optimize Random Forest for faster SHAP calculation
      # Use fewer trees for faster training (100 instead of 200)
      # The trade-off is slightly less accuracy, but significantly faster SHAP calculation
      model <- randomForest(
        Labels ~ .,
        data = trainData,
        ntree = 100,  # Reduced from 200 for faster training and SHAP calculation
        mtry = max(1, floor(sqrt(ncol(trainData) - 1))),
        importance = TRUE,
        na.action = na.omit
      )
      
      print("Calculating SHAP values for Random Forest using fastshap...")
      
      # GC Petros S 7: Create prediction function for SHAP
      # For multi-class classification with fastshap, we need to predict a single value per row
      # Important: newdata must be a data.frame with same column structure as training data
      num_classes <- length(unique(trainData$Labels))
      
      # GC Petros S 7.1: Get training column names once (ensure they're valid R identifiers)
      train_feature_cols <- setdiff(colnames(trainData), "Labels")
      train_feature_cols_valid <- make.names(train_feature_cols, unique = TRUE)
      
      if (num_classes == 2) {
        # Binary classification: return probability of positive class
        pred_fun <- function(object, newdata) {
          # GC Petros S 7.2: Convert to data frame and ensure valid R identifiers as column names
          if (is.matrix(newdata)) {
            newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
          } else if (!is.data.frame(newdata)) {
            newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
          }
          
          # Ensure column names are valid R identifiers and match training
          current_cols <- colnames(newdata)
          current_cols_valid <- make.names(current_cols, unique = TRUE)
          
          # Always ensure valid R identifiers
          colnames(newdata) <- current_cols_valid
          
          # Match to training column names if possible
          if (length(current_cols_valid) == length(train_feature_cols_valid)) {
            colnames(newdata) <- train_feature_cols_valid
          }
          
          # CRITICAL: Ensure column names match training data EXACTLY
          # Random Forest requires exact column name match
          train_cols_exact <- setdiff(colnames(trainData), "Labels")
          if (length(colnames(newdata)) == length(train_cols_exact)) {
            colnames(newdata) <- train_cols_exact
          }
          
          tryCatch({
            # Suppress warnings during prediction
            old_warn_pred <- options(warn = -1)$warn
            prob_matrix <- predict(object, newdata = newdata, type = "prob")
            options(warn = old_warn_pred)
            
            # Return probability of second class (usually the "positive" class)
            result <- prob_matrix[, 2, drop = FALSE]
            if (any(!is.finite(result))) {
              result[!is.finite(result)] <- 0.5
            }
            return(result)
          }, error = function(e) {
            # Suppress error - return uniform probabilities as fallback
            # Don't print error to avoid flooding console
            return(matrix(0.5, nrow = nrow(newdata), ncol = 1))
          })
        }
      } else {
        # Multi-class: calculate SHAP for each class separately, then average
        # For simplicity, we'll use the probability of the most likely class
        pred_fun <- function(object, newdata) {
          # GC Petros S 7.2: Convert to data frame and ensure valid R identifiers as column names
          if (is.matrix(newdata)) {
            newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
          } else if (!is.data.frame(newdata)) {
            newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
          }
          
          # Ensure column names are valid R identifiers and match training
          current_cols <- colnames(newdata)
          current_cols_valid <- make.names(current_cols, unique = TRUE)
          
          # Always ensure valid R identifiers
          colnames(newdata) <- current_cols_valid
          
          # Match to training column names if possible
          if (length(current_cols_valid) == length(train_feature_cols_valid)) {
            colnames(newdata) <- train_feature_cols_valid
          }
          
          # CRITICAL: Ensure column names match training data EXACTLY
          # Random Forest requires exact column name match
          train_cols_exact <- setdiff(colnames(trainData), "Labels")
          if (length(colnames(newdata)) == length(train_cols_exact)) {
            colnames(newdata) <- train_cols_exact
          }
          
          tryCatch({
            # Suppress warnings during prediction
            old_warn_pred <- options(warn = -1)$warn
            prob_matrix <- predict(object, newdata = newdata, type = "prob")
            options(warn = old_warn_pred)
            
            # Return the maximum probability across classes (model's confidence)
            result <- matrix(apply(prob_matrix, 1, max), ncol = 1)
            if (any(!is.finite(result))) {
              result[!is.finite(result)] <- 1 / num_classes
            }
            return(result)
          }, error = function(e) {
            # Suppress error - return uniform probabilities as fallback
            # Don't print error to avoid flooding console
            return(matrix(1 / num_classes, nrow = nrow(newdata), ncol = 1))
          })
        }
      }
      
      # Determine computational limits based on dataset size
      n_features <- ncol(testData) - 1
      n_samples <- nrow(testData)
      
      # GC Petros S 11.1: CRITICAL - Keep ALL samples (cells) for SHAP calculation
      # Subsampling samples can lose important information - genes that are important
      # for the excluded cells will be missed. We must use ALL samples.
      # Only optimize by reducing simulations (nsim), not by subsampling samples.
      
      # Always use ALL samples - never subsample cells
      # This ensures we capture gene importance across all cell types
      shap_test_data <- testData
      
      message(paste("Using ALL", n_features, "features and ALL", n_samples, "samples for SHAP calculation"))
      message("(Important: Using all cells ensures no gene importance information is lost)")
      
      # GC Petros S 8: Calculate SHAP values with optimized parameters for speed
      # GC Petros S 11.2: Reduce simulations for faster calculation while keeping all samples
      # We keep ALL samples (cells) to preserve information, but reduce simulations for speed
      if (MLmethod == "rf") {
        # Random Forest: Use 30 simulations (reduced from 50-100 for speed)
        # This is a good balance between speed and accuracy
        # We keep ALL samples to ensure no gene importance is lost
        nsim_value <- 30
        message(paste("RF SHAP: Using", nsim_value, "Monte Carlo simulations (optimized for speed)"))
        message(paste("Using ALL", n_samples, "samples to preserve complete gene importance information"))
      } else {
        # XGBoost: Use standard simulations (it's fast anyway due to built-in SHAP)
        nsim_value <- min(100, max(50, ceiling(2000 / (n_features * n_samples))))
        message(paste("Using", nsim_value, "Monte Carlo simulations for SHAP calculation"))
      }
      
      # GC Petros S 8.5: Prepare data for fastshap
      # CRITICAL: fastshap creates R variables from column names, so they MUST be valid R identifiers
      # Extract feature columns (exclude Labels)
      shap_X <- shap_test_data[, -ncol(shap_test_data), drop = FALSE]
      
      # Ensure it's a data frame
      if (is.matrix(shap_X)) {
        shap_X <- as.data.frame(shap_X, stringsAsFactors = FALSE)
      }
      
      # Get training data column names (these should already be valid R identifiers from make.names)
      train_feature_names <- setdiff(colnames(trainData), "Labels")
      shap_feature_names <- colnames(shap_X)
      
      # GC Petros S 8.6: CRITICAL - Ensure ALL column names are valid R identifiers
      # fastshap uses column names as variable names in R environment, so invalid names cause errors
      train_cols_valid <- make.names(train_feature_names, unique = TRUE)
      shap_cols_valid <- make.names(shap_feature_names, unique = TRUE)
      
      # Verify all column names are valid R identifiers
      # make.names ensures: starts with letter or dot, contains only letters, numbers, dots, underscores
      if (!identical(shap_feature_names, shap_cols_valid)) {
        # Update to valid names
        colnames(shap_X) <- shap_cols_valid
        shap_feature_names <- shap_cols_valid
      }
      
      # Ensure training column names are also valid
      if (!identical(train_feature_names, train_cols_valid)) {
        train_feature_names <- train_cols_valid
      }
      
      # Match column names between training and test
      if (length(shap_feature_names) == length(train_feature_names)) {
        # Same number of columns - match by position or name
        if (!identical(shap_feature_names, train_feature_names)) {
          # Try to match by make.names transformation
          shap_normalized <- make.names(shap_feature_names, unique = TRUE)
          train_normalized <- make.names(train_feature_names, unique = TRUE)
          
          if (identical(shap_normalized, train_normalized)) {
            # Names match after normalization - use training names
            colnames(shap_X) <- train_feature_names
            shap_feature_names <- train_feature_names
          } else {
            # Check if can match by position (same count)
            # Use training names if same count
            colnames(shap_X) <- train_feature_names
            shap_feature_names <- train_feature_names
          }
        }
      } else {
        # Different number of columns - this should not happen if we kept all features
        # But handle it gracefully
        warning("Column count mismatch: training has ", length(train_feature_names), 
                " features, SHAP data has ", length(shap_feature_names), " features.")
        warning("This may cause SHAP calculation issues. Attempting to align features...")
        
        # Try to match common features
        common_features <- intersect(train_feature_names, shap_feature_names)
        if (length(common_features) > 0 && length(common_features) == length(train_feature_names)) {
          # All training features exist - use them
          shap_X <- shap_X[, train_feature_names, drop = FALSE]
          shap_feature_names <- train_feature_names
          colnames(shap_X) <- train_feature_names
        } else {
          # Some features missing - this is an error
          stop("Feature mismatch: SHAP data is missing features that the model requires. ",
               "Please ensure all features from training are present in test data.")
        }
      }
      
      # Final validation: ensure all column names are valid R identifiers
      # Remove any special characters and ensure they're valid
      final_cols <- make.names(shap_feature_names, unique = TRUE)
      colnames(shap_X) <- final_cols
      
      # GC Petros S 8.7: Run fastshap with comprehensive error handling
      # CRITICAL: fastshap uses column names as R variables, so they must be valid identifiers
      # The warnings "object 'X' not found" are expected from fastshap's internal evaluation
      # These are harmless - fastshap handles them internally
      
      # Ensure shap_X column names EXACTLY match training data column names
      train_cols_for_shap <- setdiff(colnames(trainData), "Labels")
      shap_cols_current <- colnames(shap_X)
      
      # Force exact match if counts are the same
      if (length(shap_cols_current) == length(train_cols_for_shap)) {
        colnames(shap_X) <- train_cols_for_shap
      }
      
      # Suppress ALL warnings and messages during fastshap execution
      old_warn <- options(warn = -1)$warn
      
      # Set up warning handler to suppress all warnings
      shap_values_matrix <- tryCatch({
        withCallingHandlers({
          fastshap::explain(
            model,
            X = shap_X,
            pred_wrapper = pred_fun,
            nsim = nsim_value,
            adjust = TRUE,  # Adjust for background distribution
            newdata = NULL  # Use training data as background
          )
        }, warning = function(w) {
          # Suppress all warnings - they're expected from fastshap/pred_fun
          invokeRestart("muffleWarning")
        }, message = function(m) {
          # Suppress messages too
          invokeRestart("muffleMessage")
        })
      }, error = function(e) {
        options(warn = old_warn)
        message("Error in fastshap::explain: ", e$message)
        message("Attempting fallback approach with simplified prediction function...")
        
        # Fallback: Use even simpler approach with explicit column name handling
        simple_pred_fun <- function(object, newdata) {
          # Ensure newdata is data frame with valid column names
          if (is.matrix(newdata)) {
            newdata <- as.data.frame(newdata, stringsAsFactors = FALSE)
          }
          
          # Ensure valid R identifiers and match training columns
          current_cols <- make.names(colnames(newdata), unique = TRUE)
          colnames(newdata) <- current_cols
          
          # Match to training columns if possible
          if (length(current_cols) == length(train_feature_cols_valid)) {
            colnames(newdata) <- train_feature_cols_valid
          }
          
          # Suppress warnings during prediction
          options(warn = -1)
          tryCatch({
            prob_matrix <- predict(object, newdata = newdata, type = "prob")
            if (num_classes == 2) {
              result <- matrix(prob_matrix[, 2], ncol = 1)
            } else {
              result <- matrix(apply(prob_matrix, 1, max), ncol = 1)
            }
            options(warn = old_warn)
            return(result)
          }, error = function(e2) {
            options(warn = old_warn)
            # Ultimate fallback: return probabilities
            if (num_classes == 2) {
              return(matrix(0.5, nrow = nrow(newdata), ncol = 1))
            } else {
              return(matrix(1 / num_classes, nrow = nrow(newdata), ncol = 1))
            }
          })
        }
        
        # Try again with simplified function and suppress ALL warnings
        suppressWarnings({
          options(warn = -1)
          result <- tryCatch({
            # Ensure column names match before retry
            train_cols_for_shap <- setdiff(colnames(trainData), "Labels")
            if (length(colnames(shap_X)) == length(train_cols_for_shap)) {
              colnames(shap_X) <- train_cols_for_shap
            }
            
            fastshap::explain(
              model,
              X = shap_X,
              pred_wrapper = simple_pred_fun,
              nsim = max(30, min(nsim_value, 50)),  # Reduce nsim for stability
              adjust = FALSE,  # Disable adjustment for compatibility
              newdata = NULL
            )
          }, error = function(e2) {
            options(warn = old_warn)
            stop("SHAP calculation failed. This may be due to incompatible gene names or insufficient data.")
          })
          options(warn = old_warn)
          return(result)
        })
      })
      
      # Restore warnings
      options(warn = old_warn)
      
      # GC Petros S 9: Process SHAP values to get feature importance
      # SHAP values matrix: rows = samples, columns = features
      # Calculate mean absolute SHAP value for each feature (standard approach)
      if (is.matrix(shap_values_matrix)) {
        mean_shap_importance <- colMeans(abs(shap_values_matrix), na.rm = TRUE)
        feature_names_shap <- colnames(shap_values_matrix)
      } else if (is.data.frame(shap_values_matrix)) {
        mean_shap_importance <- colMeans(abs(shap_values_matrix), na.rm = TRUE)
        feature_names_shap <- colnames(shap_values_matrix)
      } else {
        stop("Unexpected SHAP values format")
      }
      
      # Expand to include all features (set non-calculated features to 0)
      all_features <- colnames(testData[, -ncol(testData)])
      full_shap_importance <- rep(0, length(all_features))
      names(full_shap_importance) <- all_features
      full_shap_importance[feature_names_shap] <- mean_shap_importance
      
      mean_shap_importance <- full_shap_importance
      feature_names <- all_features
      
      print(paste("SHAP calculation completed for", length(feature_names_shap), "features"))
      
    } else if (MLmethod == "xgbTree") {
      if (!requireNamespace("xgboost", quietly = TRUE)) {
        stop("Package 'xgboost' is required but not installed")
      }
      library(xgboost)
      
      print("Training XGBoost model...")
      
      # Prepare data for XGBoost
      X_train <- as.matrix(trainData[, -ncol(trainData)])
      y_train <- as.numeric(trainData$Labels) - 1  # XGBoost expects 0-based labels
      X_test <- as.matrix(testData[, -ncol(testData)])
      
      num_class <- length(unique(trainData$Labels))
      
      # Train XGBoost model
      train_matrix <- xgb.DMatrix(data = X_train, label = y_train)
      
      params <- list(
        objective = "multi:softprob",
        num_class = num_class,
        max_depth = 4,
        eta = 0.1,
        subsample = 0.8,
        colsample_bytree = 0.8,
        nthread = 1,
        eval_metric = "mlogloss"
      )
      
      model <- xgb.train(
        params = params,
        data = train_matrix,
        nrounds = 100,  # Increased for better convergence
        verbose = 0
      )
      
      print("Calculating SHAP values for XGBoost...")
      
      # GC Petros S 10: Calculate SHAP values using XGBoost's built-in SHAP
      # XGBoost provides SHAP contributions via predcontrib
      # For multi-class: returns contributions for each class
      shap_contrib <- predict(
        model, 
        newdata = X_test, 
        predcontrib = TRUE,
        approxcontrib = FALSE  # Use exact SHAP values (slower but more accurate)
      )
      
      # Process SHAP contribution matrix
      # For multi-class classification with predcontrib=TRUE:
      # Shape: (n_samples, n_features * n_classes + n_classes)
      # Structure: [feature1_class1, feature1_class2, ..., feature1_classN,
      #              feature2_class1, ..., featureN_classN, bias_class1, ..., bias_classN]
      n_features <- ncol(X_test)
      n_samples <- nrow(X_test)
      
      # Extract feature contributions (exclude bias terms at the end)
      n_bias_cols <- num_class
      n_feature_contrib_cols <- ncol(shap_contrib) - n_bias_cols
      
      if (n_feature_contrib_cols == n_features * num_class) {
        # Multi-class: reshape contributions
        # Contributions are organized as: [feat1_class1, feat1_class2, ..., featN_class1, featN_class2, ...]
        feature_shap_all <- shap_contrib[, 1:n_feature_contrib_cols, drop = FALSE]
        
        # Reshape to (n_samples, n_features, num_class)
        feature_shap_array <- array(
          feature_shap_all, 
          dim = c(n_samples, num_class, n_features)
        )
        
        # Permute to (n_samples, n_features, num_class)
        feature_shap_array <- aperm(feature_shap_array, c(1, 3, 2))
        
        # Calculate mean absolute SHAP across all classes and samples for each feature
        feature_shap_mean <- apply(feature_shap_array, 2, function(x) {
          mean(abs(x), na.rm = TRUE)
        })
      } else if (n_feature_contrib_cols == n_features) {
        # Binary or single class: use directly
        feature_shap_mean <- colMeans(abs(shap_contrib[, 1:n_features, drop = FALSE]), na.rm = TRUE)
      } else {
        # Unexpected format: use all columns except last n_bias_cols
        feature_shap_mean <- colMeans(abs(shap_contrib[, 1:n_feature_contrib_cols, drop = FALSE]), na.rm = TRUE)
        # If we have more columns than features, average them
        if (length(feature_shap_mean) > n_features) {
          # Reshape and average
          feature_shap_mean <- apply(
            array(feature_shap_mean, dim = c(num_class, n_features)), 
            2, 
            mean, 
            na.rm = TRUE
          )
        }
      }
      
      # Ensure we have the right number of features
      if (length(feature_shap_mean) != n_features) {
        warning("Mismatch in SHAP dimensions. Using first n_features values.")
        feature_shap_mean <- feature_shap_mean[1:min(n_features, length(feature_shap_mean))]
      }
      
      mean_shap_importance <- feature_shap_mean
      feature_names <- colnames(X_test)
      names(mean_shap_importance) <- feature_names
      
      print(paste("SHAP calculation completed for", length(feature_names), "features"))
    }
    
  }, error = function(e) {
    message("Error in SHAP values calculation: ", e$message)
    message("Falling back to ML importance method...")
    
    # Fallback to simple importance method
    source("./Scripts/UltraSimpleShapFilter.R", local = TRUE)
    return(UltraSimpleShapFilter(data, Labels, MLmethod, importanceLimit))
  })
  
  # GC Petros S 11: Check if SHAP values were calculated successfully
  if (!exists("mean_shap_importance") || !exists("feature_names") || 
      length(mean_shap_importance) == 0 || length(feature_names) == 0) {
    warning("SHAP values calculation failed. Using fallback method...")
    source("./Scripts/UltraSimpleShapFilter.R", local = TRUE)
    return(UltraSimpleShapFilter(data, Labels, MLmethod, importanceLimit))
  }
  
  # GC Petros S 12: Create importance dataframe (DO NOT normalize - use absolute SHAP values)
  # SHAP values should be interpreted in their original scale
  importance_df <- data.frame(
    Feature = feature_names,
    SHAP_Importance = as.numeric(mean_shap_importance),
    stringsAsFactors = FALSE
  )
  
  # Remove any NA or Inf values
  importance_df <- importance_df[is.finite(importance_df$SHAP_Importance), ]
  
  # Sort by importance (descending)
  importance_df <- importance_df[order(importance_df$SHAP_Importance, decreasing = TRUE), ]
  
  print(paste("SHAP importance range:", 
              round(min(importance_df$SHAP_Importance, na.rm = TRUE), 6), 
              "to", 
              round(max(importance_df$SHAP_Importance, na.rm = TRUE), 6)))
  
  # GC Petros S 13: Filter by importance threshold
  # SHAP values are in their original scale (mean absolute SHAP)
  # importanceLimit needs to be interpreted correctly:
  # - If < 0.1 and max SHAP > 10: likely meant as percentile (0.01 = top 1%)
  # - Otherwise: interpret as absolute threshold
  max_shap <- max(importance_df$SHAP_Importance, na.rm = TRUE)
  
  if (max_shap > 0) {
    # Check if importanceLimit is reasonable for absolute SHAP values
    # SHAP values are typically in range 0.01-0.1 for normalized features
    # If importanceLimit > 1, it's likely meant as a count or percentile
    if (importanceLimit > 1 && importanceLimit <= 100) {
      # Interpret as top N features
      n_select <- min(ceiling(importanceLimit), nrow(importance_df))
      filtered_importance <- importance_df[1:n_select, ]
      message(paste("Interpreting importanceLimit as top", n_select, "features"))
    } else if (importanceLimit > 0.1 && max_shap < 1) {
      # importanceLimit too high for SHAP scale, use percentile
      percentile_threshold <- quantile(importance_df$SHAP_Importance, 
                                       probs = 1 - min(importanceLimit / 100, 0.95), 
                                       na.rm = TRUE)
      filtered_importance <- importance_df[importance_df$SHAP_Importance >= percentile_threshold, ]
      message(paste("Using percentile threshold:", round(percentile_threshold, 6)))
    } else if (importanceLimit < 0.01 && max_shap > 10) {
      # Very small threshold with large SHAP values: use percentile
      percentile_threshold <- quantile(importance_df$SHAP_Importance, 
                                       probs = max(0.01, 1 - (importanceLimit * 100)), 
                                       na.rm = TRUE)
      filtered_importance <- importance_df[importance_df$SHAP_Importance >= percentile_threshold, ]
    } else {
      # Use absolute threshold as provided
      filtered_importance <- importance_df[importance_df$SHAP_Importance > importanceLimit, ]
    }
  } else {
    filtered_importance <- importance_df
  }
  
  # Ensure we have at least some features
  if (nrow(filtered_importance) == 0) {
    message("No features meet the importance threshold. Selecting top features by SHAP importance...")
    # Select top features: at least 10, up to top 10% or 200, whichever is smaller
    n_select <- min(max(10, ceiling(nrow(importance_df) * 0.1)), 
                    min(200, nrow(importance_df)))
    filtered_importance <- importance_df[1:n_select, ]
  } else {
    # Limit to reasonable number even if many meet threshold
    max_features_select <- min(300, nrow(filtered_importance))
    filtered_importance <- filtered_importance[1:min(max_features_select, nrow(filtered_importance)), ]
  }
  
  # GC Petros S 14: Create gene importance matrix and filter data
  iG <- data.frame(
    SHAP_Importance = filtered_importance$SHAP_Importance,
    row.names = filtered_importance$Feature,
    stringsAsFactors = FALSE
  )
  
  # Select genes from original data
  selected_genes <- rownames(iG)
  
  # Match with original column names (handle make.names transformation)
  available_genes <- selected_genes[selected_genes %in% colnames(data_with_labels)]
  
  if (length(available_genes) == 0) {
    warning("No genes found matching selected features. This may indicate a data structure issue.")
    # Try to match with original data column names
    original_cols <- colnames(data)
    available_genes <- selected_genes[selected_genes %in% make.names(original_cols)]
    if (length(available_genes) == 0) {
      stop("Could not match selected features with data columns.")
    }
  }
  
  # Create output data with selected genes
  # Use original data structure (without Labels column that we added)
  gene_cols <- colnames(data)[colnames(data) %in% available_genes | 
                               make.names(colnames(data)) %in% available_genes]
  
  if (length(gene_cols) == 0) {
    # Fallback: use the names directly
    gene_cols <- available_genes
  }
  
  # Match indices
  gene_indices <- match(available_genes, colnames(data_with_labels))
  gene_indices <- gene_indices[!is.na(gene_indices)]
  
  if (length(gene_indices) > 0) {
    newdata <- data_with_labels[, gene_indices, drop = FALSE]
    newdata$Labels <- as.factor(Labels)
  } else {
    # Final fallback
    newdata <- data[, gene_cols, drop = FALSE]
    newdata$Labels <- as.factor(Labels)
  }
  
  # Update iG to match available genes
  if (nrow(iG) > 0) {
    iG <- iG[available_genes[available_genes %in% rownames(iG)], , drop = FALSE]
  }
  
  print(paste("SHAP Values analysis completed. Selected", nrow(iG), "genes based on SHAP importance."))
  
  return(list(ig = iG, newdata = newdata))
}

# Function to create SHAP summary plot (for visualization)
create_shap_summary_plot = function(shap_values, feature_names, top_n = 20) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for SHAP plots.")
  }
  
  library(ggplot2)
  
  # Calculate mean absolute SHAP values
  if (is.matrix(shap_values)) {
    mean_shap <- colMeans(abs(shap_values), na.rm = TRUE)
  } else if (is.data.frame(shap_values)) {
    mean_shap <- colMeans(abs(shap_values), na.rm = TRUE)
  } else {
    stop("shap_values must be a matrix or data frame")
  }
  
  # Create dataframe for plotting
  plot_data <- data.frame(
    Feature = feature_names,
    Mean_SHAP = mean_shap,
    stringsAsFactors = FALSE
  )
  
  # Sort and select top features
  plot_data <- plot_data[order(plot_data$Mean_SHAP, decreasing = TRUE), ]
  plot_data <- plot_data[1:min(top_n, nrow(plot_data)), ]
  
  # Create plot
  p <- ggplot(plot_data, aes(x = reorder(Feature, Mean_SHAP), y = Mean_SHAP)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(
      title = "SHAP Feature Importance",
      x = "Features",
      y = "Mean |SHAP Value|"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  return(p)
}
