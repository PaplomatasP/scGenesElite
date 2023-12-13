

#Function to create SingleCellExperiment  object
data_fitting = function(data, group, normal.Method) {
  if (missing(normal.Method)) {
    normcounts <- data
    
  }
  else{
    normcounts <- normaliz_method(data, normal.Method)
  }
  SingleCell <-
    SingleCellExperiment::SingleCellExperiment(assays = list(counts = data, normcounts = normcounts))
  gene_df <- data.frame(Gene = rownames(SingleCell))
  cell_df <- data.frame(label = group, cell = colnames(SingleCell))
  rownames(gene_df) <- gene_df$Gene
  rownames(cell_df) <- cell_df$cell
  SingleCell <-
    SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = data, normcounts = normcounts),
      colData = cell_df,
      rowData = gene_df
    )
}



#Calculates the P-value according to the requested method

PvalueCalc = function(data, Pvaluemethod) {
  if (Pvaluemethod == "Seurat_method") {
    object_Seurat <-
      as.matrix(SummarizedExperiment::assay(data, "counts"))
    
    tLabel <- data$label
    
    names(tLabel) <- colnames(data)
    MetaData <- data.frame(groups = tLabel)
    
    Input <-
      Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = MetaData)
    Input <-NormalizeData(Input)
   
    # Determine the number of levels in Labels
    num_levels <- length(levels(as.factor(tLabel)))
    
    # Create a list of ident values based on the number of levels
    ident_values <- levels(as.factor(tLabel))[1:num_levels]
    
    # Create a list of ident arguments for Seurat::FindMarkers
    ident_args <- paste0("ident.", 1:num_levels, " = ident_values[", 1:num_levels, "]")
    
    # Combine the ident arguments into a single string
    ident_string <- paste(ident_args, collapse = ", ")
    
    # Use the generated ident_string in Seurat::FindMarkers
    result <- eval(parse(text = paste("Seurat::FindMarkers(object = Input, 
                                    group.by = 'groups', logfc.threshold = -Inf, 
                                    test.use = 'wilcox', only.pos = FALSE, verbose = FALSE,
                                    min.cells.group = 1, ", ident_string, ")")))
    
    
    results_Seurat <- list(
      gene_names = row.names(result),
      pvalue = result$p_val,
      FDR = result$p_val_adj,
      logFC = result$avg_log2FC
    )
    
    return(results_Seurat)
  }
  
  
  
  
  
  else if (Pvaluemethod == "BPSC_metchod") {
    object_BPSC <-
      as.matrix(SummarizedExperiment::assay(data, "normcounts"))
    
    controllds <- which(data$label == levels(factor(data$label)[1]))
    LABELS <- data$label
    unique_conditions <- unique(data$label)
    norm_counts <- as.data.frame(t(object_BPSC))
    norm_counts$labels <- as.factor(LABELS)
    DiffMEan <- data.frame(t(aggregate(. ~ labels, data = norm_counts, FUN = mean)))
    colnames(DiffMEan) <- DiffMEan[1,]
    DiffMEan <- DiffMEan[-1,]
    DiffMEan <- as.data.frame(lapply(DiffMEan, as.numeric) )
    logFC_values <- log2(DiffMEan[, 2] / DiffMEan[, 1])
    
    design <- model.matrix( ~ data$label)
    resbp <-
      BPSC::BPglm(
        data = object_BPSC,
        controlIds = controllds,
        design = design,
        coef = 2,
        estIntPar = TRUE,
        useParallel = TRUE
      )
    FDR <- p.adjust(resbp$PVAL, method = "BH")
    result_BPSC <-
      list(
        gene_names = names(resbp$PVAL),
        pvalue = resbp$PVAL,
        FDR = FDR,
        logFC = logFC_values
      )
    return(result_BPSC)
  }
  
  
  
  else if (Pvaluemethod == "MAST_method") {
    
    options(warn = -1)
    object_MAST <-
      as.matrix(SummarizedExperiment::assay(data, "normcounts"))
    
    grp <- data$label
    names(grp) <- colnames(object_MAST)
    
    sca <- MAST::FromMatrix(
      exprsArray = log2(object_MAST + 1 + 1e-6),  # add a small pseudocount
      cData = data.frame(wellKey = names(grp),
                         grp = grp)
    )
    
    zlmdata <-
      MAST::zlm( ~ grp, sca, method = "bayesglm", parallel = TRUE)
    
    # Fetch coefficients
    coef_df <- summary(zlmdata)$datatable
    coef_df <- subset(coef_df, contrast != "(Intercept)")
    mast <- MAST::lrTest(zlmdata, "grp")
    mast_order <- match(rownames(mast), coef_df$primerid)
    coef_df <- coef_df[mast_order, ]
    # Get log fold changes
    logFC_values <- coef_df$coef
    
    
    FDR  <- p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH")
    
    result_MAST <-
      list(gene_names = names(mast[, "hurdle", "Pr(>Chisq)"]),
           pvalue = mast[, "hurdle", "Pr(>Chisq)"],
           FDR = FDR,
           logFC = logFC_values)
    return(result_MAST)
  }
  
  
  else if (Pvaluemethod == "DESeq2_method") {
    # perform DESeq2
    
    options(warn = -1)
    object_DESeq2 <- as.matrix(SummarizedExperiment::assay(data, "counts"))
    genes=rownames(object_DESeq2)
    object_DESeq2 <- data.frame(apply(object_DESeq2, 2, function(x) as.numeric(as.character(x))))
    rownames(object_DESeq2)=genes
    object_DESeq2=na.omit(object_DESeq2)
  
    object_DESeq2 <- DESeq2::DESeqDataSetFromMatrix(countData = round(object_DESeq2 + 1),
                                                    colData = data.frame(condition = factor(data$label)),
                                                    design = ~ condition)
    
    
    object_DESeq2 <- DESeq2::DESeq(object_DESeq2, test = "LRT", parallel = TRUE, betaPrior = FALSE, fitType = "parametric", reduced = ~ 1)
    
    
    res_cpm <- DESeq2::results(object_DESeq2)
    result_DESeq2 <- list(gene_names = rownames(res_cpm),
                          pvalue = res_cpm$pvalue,
                          FDR = res_cpm$padj,
                          logFC = res_cpm$log2FoldChange)
  
    return(result_DESeq2)
  }
  
  
  
  
}

#Calculates the P-value according to the requested method and returns a data frame with the expression table of the isolated genes according to the threshold.

StatisticalPvalueFilter = function(data, Labels, threshold,logfc) {
  options(scipen = 999)
  obj=data
  #Labels=data[,ncol(data)]
  data = as.data.frame(t(data[, -ncol(data)]))
  data1 = data_fitting(data, Labels)
  
  
  PvalueMethod <-
    ifelse(input$P_method == "Empty",
           input$ensemblePvalue,
           input$P_method)
  print(PvalueMethod)
  PvalueData <- PvalueCalc(data1, Pvaluemethod = PvalueMethod) #PvalueMethod
  
  PvalueData  <- as.data.frame(PvalueData)
  if (PvalueMethod == "MAST_method" |
      PvalueMethod == "BPSC_metchod" | PvalueMethod == "Seurat_method") {
    print("MAST or BPSC or Seurat")
    rownames(PvalueData) = PvalueData[, 1]
    count = 0
    PvalueTreshold = list()
    for (i in seq(1, nrow(PvalueData))) {
      if (!is.na(PvalueData[, 3][i]) && PvalueData[, 3][i] <= threshold && abs(PvalueData[, 4][i]) >= logfc )  {
        count = count + 1
        
        PvalueTreshold[[count]] = rownames(PvalueData)[i]
      }
    }
    PvalueTreshold = as.data.frame(cbind(PvalueTreshold))
    if (length(PvalueTreshold) != 0){
    row.names(PvalueTreshold) = PvalueTreshold[, 1]
    
    PvalueData1 <-
      data.frame(PvalueData[order(PvalueData$pvalue, decreasing = FALSE), drop = FALSE,])
    iG1 = subset(PvalueData1 ,
                 rownames(PvalueData) %in% rownames(PvalueTreshold))
    iG = cbind((iG1[, -c(1, 2,4)]))
    rownames(iG) = rownames(iG1)
    
    iG <- as.data.frame(round(iG, digits = 4))
    data = as.data.frame(t(data))
    neudata = data[, colnames(data) %in% row.names(PvalueTreshold)]
    
    neudata = as.data.frame(neudata)
    
    neudata$Labels = as.factor(Labels)
    return(list(ig=iG, newdata=neudata))
    }else{ showModal(modalDialog(
      title = "Message",
      "Please consider increasing the P-value threshold!!!",
      easyClose = TRUE,footer = modalButton("OK")
    ))
    }
  }
  if (input$P_method == "DESeq2_method") {
    print("DESeq2_method")
    count = 0
    PvalueTreshold = list()
    print(nrow(PvalueData))
   
    for (i in 1:nrow(PvalueData)) {
     
     
      if (!is.na(PvalueData[, 3][i]) && PvalueData[, 3][i] <= threshold  && abs(PvalueData[, 4][i]) >= logfc  ) {  #threshold logfc
        print(i)
        count = count + 1

        PvalueTreshold[[count]] = rownames(PvalueData)[i]
      }
    }
    PvalueTreshold <- as.data.frame(cbind(PvalueTreshold))
    if (length(PvalueTreshold) != 0){
     
    PvalueData <-as.data.frame(PvalueData)
    
    
    PvalueData1 <- data.frame(PvalueData[order(PvalueData$pval, decreasing = FALSE), drop = FALSE,])
    
    iG1 = subset(PvalueData1 , rownames(PvalueData1) %in% PvalueTreshold[, 1])
    iG = cbind((iG1[, -c(1, 2, 4)]))
    rownames(iG) = rownames(iG1)
    
    iG <- as.data.frame(round(iG, digits = 4))

    neudata <- obj[, colnames(obj) %in% row.names(iG)]
    neudata<-as.data.frame(neudata)

    
    neudata$Labels = as.factor(Labels)
   
    return(list(ig=iG, newdata=neudata))
    }else{ showModal(modalDialog(
      title = "Message",
      "Please consider increasing the P-value threshold!!!",
      easyClose = TRUE,footer = modalButton("OK")
    ))
    }
    
  }

}
