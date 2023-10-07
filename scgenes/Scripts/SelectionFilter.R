
#Function to run the ImportanceFilter function and return  a dataframe the expression matrix of the isolated genes.

SelectionFilter = function(data,
                           Labels,
                           PvalueNum ,logfc
                           ) {
  if (input$VarFilter!="Unselect" ){print("Variance Filter Activated")
    MSET_VarFilter <- VarReductionFilter(data, Labels)
    #MSET_VarFilter = MSET_VarFilter[,-(ncol(MSET_VarFilter))]
    
    MSET_PvalueFilter <- StatisticalPvalueFilter(MSET_VarFilter,Labels, threshold = PvalueNum,logfc = logfc)
    #MSET_PvalueFilter = MSET_PvalueFilter[,-(ncol(MSET_PvalueFilter))]
   
  }
  else{
    MSET_PvalueFilter <- StatisticalPvalueFilter(data,Labels, threshold = PvalueNum,logfc= logfc)
    
  }
  
  return(MSET_PvalueFilter)
}

SelectionFilter1 = function(data,
                           Labels,
                           MLmethod) {
  
  if (input$VarFilter!="Unselect") {print("Variance Filter Activated")
   
  MSET_VarFilter <- VarReductionFilter(data, Labels)
  
  df_impFilter = ImportanceFilter(MSET_VarFilter,
                                  Labels,
                                  MLmethod)
  
  }
else{
  
  df_impFilter = ImportanceFilter(data,
                                  Labels,
                                  MLmethod)
  
}
return(df_impFilter)
}

