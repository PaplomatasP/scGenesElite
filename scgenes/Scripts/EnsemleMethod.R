
# Thus functin give the opportunity to user, to create an ensemble method choosing of all available methods.
# return a expression matrix of the isolated genes.
#using the borda voting method for ensemble purpose.
EnsemleMethod = function(obj, EnseLabels) {
  GenesList = list()
  ####Variable METHODS
  
  if (input$ensembleVar == "SCMarker") {
    DATA <<- SCMarkerfun(
      obj,
      Labels = EnseLabels,
      GeneSK = input$geneK,
      CellSK = input$cellK
    )
    GenesList = append(GenesList, list(DATA = colnames(DATA)))
  }
  if (input$ensembleVar == "SelfE") {
    DATA1 <- SelfEGenes(obj,
                        Labels =EnseLabels,
                        input$n)
    GenesList = append(GenesList, list(DATA1 = colnames(DATA1)))
  }
  if (input$ensembleVar == "DUBStepR") {
    DATA2 <<- DUBStepRfun(obj,Labels = EnseLabels)
    
    GenesList = append(GenesList, list(DATA2 = colnames(DATA2)))
  }
  if (input$ensembleVar == "ScPNMF") {
    DATA3 <<- scPNMFfun(obj,
                       Labels = EnseLabels,
                       DM <- input$distMethod)
    GenesList = append(GenesList, list(DATA3 = colnames(DATA3)))
  }
  if (input$ensembleVar == "M3Dropfun") {
    DATA8 <<- M3Dropfun(obj,
                       Labels = EnseLabels)
    GenesList = append(GenesList, list(DATA8 = colnames(DATA8)))
  }
  
  ####PVALUE METHODS
  
  if (input$ensemblePvalue == "Seurat_method") {
    DATA4 <<- SelectionFilter(
      obj,
      Labels = EnseLabels,
      PvalueNum = input$PvalueNum
    )
    GenesList = append(GenesList, list(DATA4 = colnames(DATA4)[-ncol(DATA4)]))
  }
  if (input$ensemblePvalue == "BPSC_metchod") {
    DATA5 <<- SelectionFilter(
      obj,
      Labels = EnseLabels,
      PvalueNum = input$PvalueNum
    )
    GenesList = append(GenesList, list(DATA5 = colnames(DATA5)[-ncol(DATA5)]))
  }
  if (input$ensemblePvalue == "MAST_method") {
    DATA6 <<- SelectionFilter(
       obj,
      Labels = EnseLabels,
      PvalueNum = input$PvalueNum
    )
    GenesList = append(GenesList, list(DATA6 = colnames(DATA6)[-ncol(DATA6)]))
  }
  if (input$ensemblePvalue == "monocle_method") {
    DATA7 <<- SelectionFilter(
       obj,
      Labels = EnseLabels,
      PvalueNum = input$PvalueNum
    )
    GenesList = append(GenesList, list(DATA7 = colnames(DATA7)[-ncol(DATA7)]))
  }
  
  ####ML METHODS
  MLlist = c(
    "rf",
    "knn",
    "svmRadial",
    "rpart",
    "C5.0",
    "glmnet",
    "lda",
    "bartMachine",
    "adaboost",
    "treebag",
    "xgbTree"
  )
  if ((input$ensembleWrapper %in% MLlist)) {
    #MLmethodiS<-ifelse(input$Wrapper_ML_Method=="Empty", input$ML_Method,input$Wrapper_ML_Method)
    
    DATA8 <<- SelectionFilter1(
      data = obj,
      Labels = EnseLabels,
      MLmethod = input$ensembleWrapper
    )
    GenesList = append(GenesList, list(DATA8 = colnames(DATA8)))
  }
  
  if (input$ensembleMLBased %in% MLlist) {
    DATA9 <<- SelectionFilter1(
      data = obj,
      Labels = EnseLabels,
      MLmethod =  input$ensembleMLBased
    )
    GenesList = append(GenesList, list(DATA9 = colnames(DATA9)))
  }
  
  
  GenesList <<- GenesList
  

    if (length(GenesList)>=2){
      
      
      
      #
      vote <-
        create_vote(GenesList, xtype = 3, candidate = c(unlist(GenesList)))
      y <<- borda_method(vote, modified = TRUE)
      # 
      # RankingResults=as.data.frame(unique(cbind(y[["other_info"]][["count_max"]])) )
      # rownames(RankingResults)=unique(y[["candidate"]])
      
      RankingResults=as.data.frame(cbind(sort(y[["other_info"]][["count_max"]],decreasing = TRUE ) )  )
      # RankingResultsDF=data.frame(RankingResults[1:200,])
      # rownames(RankingResultsDF)=rownames(RankingResults)[1:200]
      
     # RankingResults <<- unique(data.frame( y[["other_info"]][["count_max"]]))
      #genesnames=unique(names(y[["other_info"]][["count_max"]]))
      #rownames(RankingResults) = make.names(names(y[["other_info"]][["count_max"]]), unique = TRUE)
      #rownames(RankingResults) =unique(genesnames)
      
      # RankingResultsDF <- RankingResults %>%
      #   tibble::rownames_to_column() %>%
      #   arrange(desc(RankingResults[, 1])) %>%
      #   mutate(Ranking = 1:nrow(RankingResults))
     
      
   print("here")
      newdata <- as.data.frame(RankingResults[1:input$genes,]) 
      row.names(newdata) = rownames(RankingResults)[1:input$genes]

      rownames(newdata)=ifelse(substr(rownames(newdata), 1, 1) == "X",substr(rownames(newdata), 2, nchar(rownames(newdata))), rownames(newdata))
      iG <- newdata
      genesNames=colnames(obj)[which(colnames(obj) %in% rownames(iG) ==TRUE)]
      Genes=rownames(iG)[which(rownames(iG) %in% genesNames) ]
      
      newdata <- obj[,  colnames(obj) %in% rownames(iG)]
    
      iG=data.frame(iG[which(rownames(iG) %in% genesNames),] )
      rownames(iG)=Genes
      newdata <- cbind(as.data.frame(newdata), as.data.frame(as.factor(EnseLabels)) )
      print("here!!!!")
      
      #newdata$Labels = as.factor(newdata$EnseLabels)
      print("Ensemble Method is Executed!!!")
      
      
      FilterData<<-newdata
      
      return(FilterData)
      
    }else {
      
      shinyalert(title = "Message",type = "error",
                 text = "No method has been selected! Multiple methods must be selected in order to use the Ensemble method", closeOnClickOutside = TRUE) 
   
      
    }
  EnsembledData=EnsembledData
  return(EnsembledData)
 
}
