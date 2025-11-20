# Install required packages for SHAP Values functionality
# Run this script to install the necessary packages

cat("Installing required packages for SHAP Values...\n")

# Install CRAN packages
required_packages <- c("fastshap", "xgboost", "randomForest")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, dependencies = TRUE)
  } else {
    cat(pkg, "is already installed.\n")
  }
}

# Load libraries to verify installation
cat("\nLoading libraries to verify installation...\n")
library(fastshap)
library(xgboost)
library(randomForest)

cat("\nSHAP Values packages installed successfully!\n")
cat("Required packages:\n")
cat("- fastshap: For SHAP values calculation\n")
cat("- xgboost: For XGBoost model training\n")
cat("- randomForest: For Random Forest model training\n")

cat("\nYou can now use SHAP Values in scGenesElite!\n")
