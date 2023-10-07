
# Thus functin give the opportunity to user, to create an ensemble method choosing of all available methods.
# return a expression matrix of the isolated genes.
#using the borda voting method for ensemble purpose.
EnsemleMethod = function(obj, EnseLabels) {
  GenesList = list()
  ####Variable METHODS
 
  if (input$ensembleVar == "SCMarker") {
    DATA <- SCMarkerfun(
      obj,
      Labels = EnseLabels,
      GeneSK = input$geneK,
      CellSK = input$cellK
    )
    DATA <- DATA$newdata
    GenesList = append(GenesList, list(DATA = colnames(DATA)[-ncol(DATA)]))
  }
  if (input$ensembleVar == "SelfE") {
    DATA1 <- SelfEGenes(obj,
                        Labels =EnseLabels,
                        input$n)
    DATA1 <-DATA1$newdata
    GenesList = append(GenesList, list(DATA1 = colnames(DATA1)[-ncol(DATA1)]))
  }
  if (input$ensembleVar == "DUBStepR") {
    DATA2 <- DUBStepRfun(obj,Labels = EnseLabels)
    DATA2 <-DATA2$newdata
    GenesList = append(GenesList, list(DATA2 = colnames(DATA2)[-ncol(DATA2)]))
  }
  if (input$ensembleVar == "ScPNMF") {
    DATA3 <- scPNMFfun(obj,
                       Labels = EnseLabels,
                       DM <- input$distMethod)
    DATA3 <-DATA3$newdata
    GenesList = append(GenesList, list(DATA3 = colnames(DATA3)[-ncol(DATA3)]))
  } 
  if (input$ensembleVar == "M3Drop") {
    print( "i am here")
    DATA8 <- M3Dropfun(obj,
                       Labels = EnseLabels)
    DATA8 <-DATA8$newdata
    GenesList = append(GenesList, list(DATA8 = colnames(DATA8)[-ncol(DATA8)]))
  }
  
  ####PVALUE METHODS
  
  if (input$ensemblePvalue == "Seurat_method") {
    DATA4 <- SelectionFilter(
      obj,
      Labels = EnseLabels,
      PvalueNum = input$PvalueNum,logfc = input$logfc
    )
    DATA4 <- DATA4$newdata
    GenesList = append(GenesList, list(DATA4 = colnames(DATA4)[-ncol(DATA4)]))
  }
  if (input$ensemblePvalue == "BPSC_metchod") {
    DATA5 <- SelectionFilter(
      obj,
      Labels = EnseLabels,
      PvalueNum = input$PvalueNum,logfc = input$logfc
    )
    DATA5 <- DATA5$newdata
    GenesList = append(GenesList, list(DATA5 = colnames(DATA5)[-ncol(DATA5)]))
  }
  if (input$ensemblePvalue == "MAST_method") {
    DATA6 <- SelectionFilter(
       obj,
      Labels = EnseLabels,
      PvalueNum = input$PvalueNum,logfc = input$logfc
    )
    DATA6 <- DATA6$newdata
    GenesList = append(GenesList, list(DATA6 = colnames(DATA6)[-ncol(DATA6)]))
  }
  if (input$ensemblePvalue == "monocle_method") {
    DATA7 <- SelectionFilter(
       obj,
      Labels = EnseLabels,
      PvalueNum = input$PvalueNum
    )
    DATA7 <- DATA7$newdata
    GenesList = append(GenesList, list(DATA7 = colnames(DATA7)[-ncol(DATA7)]))
  }
  
  ####ML METHODS
  MLlist = c(
    "rf",
    "rpart",
    "C5.0",
    "glmnet",
    "bartMachine",
    "adaboost",
    "treebag",
    "xgbTree"
  )
  if ((input$ensembleWrapper %in% MLlist)) {
    #MLmethodiS<-ifelse(input$Wrapper_ML_Method=="Empty", input$ML_Method,input$Wrapper_ML_Method)
    
    DATA9 <- SelectionFilter1(
      data = obj,
      Labels = EnseLabels,
      MLmethod = input$ensembleWrapper
    )
    DATA9 <- DATA9$newdata
    GenesList = append(GenesList, list(DATA9 = colnames(DATA9)[-ncol(DATA9)]))
  }
  
  # if (input$ensembleMLBased %in% MLlist) {
  #   DATA9 <<- SelectionFilter1(
  #     data = obj,
  #     Labels = EnseLabels,
  #     MLmethod =  input$ensembleMLBased
  #   )
  #   GenesList = append(GenesList, list(DATA9 = colnames(DATA9)))
  # }
  
  
  GenesList <- GenesList
  # Combine all lists into one vector of unique gene names
  all_genes <- unique(unlist(GenesList))
  
    if (length(GenesList)>=2){
      
      
      # Create a vote object
      vote <- create_vote(GenesList, xtype = 3, candidate = all_genes)
      
      # Conduct Borda count
      y <- borda_method(vote, modified = TRUE)
       
      
      RankingResults=as.data.frame(cbind(sort(y[["other_info"]][["count_max"]],decreasing = TRUE ) )  )
 
      newdata <- as.data.frame(RankingResults) 
      row.names(newdata) = rownames(RankingResults)

      rownames(newdata)=ifelse(substr(rownames(newdata), 1, 1) == "X",substr(rownames(newdata), 2, nchar(rownames(newdata))), rownames(newdata))
      iG <- newdata
      genesNames=colnames(obj)[which(colnames(obj) %in% rownames(iG) ==TRUE)]
      Genes=rownames(iG)[which(rownames(iG) %in% genesNames) ]
      
      newdata <- obj[,  colnames(obj) %in% rownames(iG)]
    
      iG=data.frame(iG[which(rownames(iG) %in% genesNames),] )
      rownames(iG)=Genes
      iG <- iG
      newdata <- cbind(as.data.frame(newdata), as.data.frame(as.factor(EnseLabels)) )
      print("here!!!!")
      
      print("Ensemble Method is Executed!!!")
      
      
      FilterData <- newdata
      
      return(list(ig=iG, newdata=newdata))
      
    }else {
      
      shinyalert(title = "Message",type = "error",
                 text = "No method has been selected! Multiple methods must be selected in order to use the Ensemble method", closeOnClickOutside = TRUE) 
   
      
    }
  EnsembledData=EnsembledData
  return(EnsembledData)
}