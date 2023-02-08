InstallCran=function(){
  
  CranPackages=c("shinythemes","caret","enrichR",
                 "gridExtra","shinyWidgets","shinydashboard","Seurat",
                 "visNetwork","grid","png","DT","shinyjs", "votesys",
                 "shinycustomloader","igraph", "BiocManager", "shinyalert", "randomForest", "xgboost", "C50", "Matrix", "glmnet")

install.packages(CranPackages)

}

InstallCran()
