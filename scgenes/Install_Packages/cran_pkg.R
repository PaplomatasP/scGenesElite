InstallCran=function(){
  
  CranPackages=c("shinythemes","caret","enrichR","fastshap","SeuratObject",
                 "gridExtra","shinyWidgets","shinydashboard","Seurat",
                 "visNetwork","grid","png","DT","shinyjs", "votesys",
                 "shinycustomloader","igraph", "BiocManager", "shinyalert", "randomForest", "xgboost", "C50", "Matrix", "glmnet",
                 "doParallel", "foreach")

install.packages(CranPackages)

}

InstallCran()
