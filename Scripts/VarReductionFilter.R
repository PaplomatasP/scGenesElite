#Function to find the variance each variable and keep only them according the threshold.

VarReductionFilter = function(data, Labels) {
  # FIlter according the FindVariableFeatures function of seraut . Keep only the higher variable genes according the threshold.

    print("Running Strict_Filter")
    
    
    MSET_VarFilter = data[, -caret::nearZeroVar(data, uniqueCut = input$uniqueCut, )[-ncol(data)]]
    
    neudata = as.data.frame(MSET_VarFilter)
    
    neudata$Labels = as.factor(Labels)
    
  
  return(neudata)
}
