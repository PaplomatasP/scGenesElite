library(shiny)
library(ggiraph)
source('Scripts/InputValidation.R')
# Exercise both upload branches without rerunning a biological selector.
fixture <- data.frame(A=1:40,B=2:41,C=3:42,Labels=rep(c('x','y'),each=20))
testServer(function(input,output,session) {
  read_uploaded_csv <- function(...) fixture
  read_uploaded_rds <- function(...) fixture
  LexikonFun <- function(...) fixture[,c('A','B','Labels')]
  SCMarkerfun <- function(data,...) list(ig=data.frame(score=c(5,4),row.names=c('A','B')),newdata=data)
  source('Scripts/MainSelectionFun.R',local=TRUE)
}, {
  session$setInputs(click=1,file1=list(name='test.csv'),GENEid='SYMBOL',organismus='Mouse',genes=20,
    VariableM='SCMarker',P_method='Empty',ML_Method='Empty',SHAP_Method='Empty',
    ensembleVar='NoMethod',ensemblePvalue='NoMethod',ensembleWrapper='NoMethod')
  result <- MethodData()
  stopifnot(result$plot_context$input_genes==3L,result$plot_context$mapped_genes==2L,
    result$plot_context$method=='SCMarker',nrow(result$plot_context$expression)==40L)
  session$setInputs(file1=NULL,rdsFile=list(name='test.rds'))
  result <- MethodData()
  stopifnot(result$plot_context$input_genes==3L,result$plot_context$mapped_genes==2L)
})
cat('Publication input context: CSV and RDS paths preserve input/mapping counts.\n')
