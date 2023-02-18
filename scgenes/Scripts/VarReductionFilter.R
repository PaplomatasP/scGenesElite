#Function to find the variance each variable and keep only them according the threshold.

VarReductionFilter = function(data, Labels) {
  # FIlter according the FindVariableFeatures function of seraut . Keep only the higher variable genes according the threshold.
  
  
  
  if (input$VarFilter == "Lightweight_Filter") {
    
    colnames(data) <- make.names(colnames(data))

    Ndata = as.data.frame(t(data[, -ncol(data)]))
    dim(Ndata)
    Seraut_obj <- CreateSeuratObject(Ndata)
    DataSeurat <-
      FindVariableFeatures(Seraut_obj,
                           selection.method = "vst",
                           nfeatures = input$nfeatures)
    Genes <- head(VariableFeatures(DataSeurat), input$nfeatures)
    neudata <- as.data.frame(data[, colnames(data) %in% Genes])
    neudata$Labels <- as.factor(Labels)
  }
  # Filter according the zero variance filter using a cutoff threshold.
  else {
    print("Running Strict_Filter")
    
    
    MSET_VarFilter = data[, -caret::nearZeroVar(data, uniqueCut = input$uniqueCut, )[-ncol(data)]]
    
    neudata = as.data.frame(MSET_VarFilter)
    
    neudata$Labels = as.factor(Labels)
    
  }
  return(neudata)
}
