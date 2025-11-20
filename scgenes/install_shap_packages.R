# Installation script for SHAP Values dependencies
# Run this script to install the required packages for SHAP functionality

# Install required packages for SHAP values
if (!requireNamespace("fastshap", quietly = TRUE)) {
  install.packages("fastshap")
}

if (!requireNamespace("xgboost", quietly = TRUE)) {
  install.packages("xgboost")
}

if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}

# Load libraries to verify installation
library(fastshap)
library(xgboost)
library(randomForest)

cat("SHAP Values packages installed successfully!\n")
cat("Required packages:\n")
cat("- fastshap: For SHAP values calculation\n")
cat("- xgboost: For XGBoost model training\n")
cat("- randomForest: For Random Forest model training\n")
