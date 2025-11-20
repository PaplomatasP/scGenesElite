

###SCMarker####
# Use the SCMarker feature selection method  and return a dataframe the expression matrix of the isolated genes.
SCMarkerfun = function(data, GeneSK, CellSK, Labels) {
  
  obj = as.matrix(data[, -ncol(data)])
  if (input$VarFilter != "Unselect") {
    print("Variance Filter Activated")
    MSET_VarFilter <- VarReductionFilter(obj, Labels)
    obj = MSET_VarFilter[, -ncol(MSET_VarFilter)]
    if (input$Norm == "Normal") {
      scale.factor <- mean(colSums(obj))
      
      obj <- Seurat::LogNormalize(obj,
                                  scale.factor = scale.factor)
    }
    
  }
  else{
    if (input$Norm == "Normal") {
      scale.factor <- mean(colSums(obj))
      
      obj <- Seurat::LogNormalize(obj,
                                  scale.factor = scale.factor)
    }
    
  }
  
  
  obj <- as.data.frame(t(obj))
 # colnames(obj) = Labels

  res <-
    ModalFilter(
      data = obj,
      geneK = 10,
      cellK = 10,
      cutoff = 2,
      width = 1
    )# default width = 1 for UMI data, width =2 for TPM data.
  

  
  res = SCMarker::GeneFilter(obj = res)
  res= SCMarker::getMarker(obj=res,k=150,n=20)
  genes = res[["marker"]]
  # Score = res[["geneSumm"]][["count"]]
  Score=cbind(res[["geneSumm"]][["gene"]],res[["geneSumm"]][["count"]])
  Score=Score[which(Score[,1] %in% genes),]
  
  iG <- cbind(as.numeric(Score[,2] ) )
  rownames(iG) = Score[,1]
  iG = iG[order(iG[, 1], decreasing = TRUE), ]
  
  
  iG <- as.data.frame(iG)
  
  
  newdata <- data[, colnames(data) %in%  genes]
  newdata$Labels <- as.factor(Labels)
  
  return(list(ig=iG, newdata=newdata))
}


#### VST ##### # not the selfE

SelfEGenes = function(obj, n, Labels) {
  obj = obj[, -ncol(obj)]
  if (input$VarFilter != "Unselect") {
    print("Variance Filter Activated")
    MSET_VarFilter <- VarReductionFilter(obj, Labels)
    obj = MSET_VarFilter[, -ncol(MSET_VarFilter)]
    print("var")
    if (input$Norm == "Normal") {
      scale.factor <- median(colSums(obj))
      
      obj <- Seurat::LogNormalize(obj,
                                  scale.factor = scale.factor)
      print("Norm")
    }
    
  }
  else {
    if (input$Norm == "Normal") {
      scale.factor <- mean(colSums(obj))
      
      obj <- Seurat::LogNormalize(obj,
                                  scale.factor = scale.factor)
      print("no Var")
      
    }
  }

  obj <- as.data.frame(t(obj) )
  #rownames(obj) <- make.names(rownames(obj))

  Seraut_obj <- CreateSeuratObject(obj)

  DataSeurat <-
    Seurat::FindVariableFeatures(Seraut_obj,
                         selection.method = "vst")

  variable_feature_data <- DataSeurat@assays[["RNA"]]@meta.data
  variable_feature_data <- variable_feature_data[order(variable_feature_data$vf_vst_counts_variance.standardized,decreasing = TRUE),]
  variable_feature_data <- head(variable_feature_data, n)
  Genes <- variable_feature_data$var.features
  
  iG <- data.frame(cbind(Genes,round(variable_feature_data$vf_vst_counts_variance.standardized,3) ) )
  # Set row names to be the gene names and remove the gene name column
  rownames(iG) <- iG$Genes
  iG$Genes <- NULL
  # Convert the remaining column to numeric
  iG[[1]] <- as.numeric(iG[[1]])
  # Rename the column
  names(iG) <- "vst.variance.standardized"
  iG <- iG
  obj <- as.data.frame(t(obj))
  neudata <- as.data.frame(obj[, colnames(obj) %in% Genes])
  neudata$Labels <- as.factor(Labels)
  
 
  return(list(ig=iG, newdata=neudata))
}

# 
# # Use the DUBStepR feature selection method and return a dataframe the expression matrix of the isolated genes.
# # DUBStepR replace  from scran , let just name same!
# DUBStepRfun = function(data, Labels) {
#   obj = as.matrix(data[, -ncol(data)])
#   a <-obj
#   obj<<-a
#   Labels <- Labels
#   l <<-Labels
#   if (input$VarFilter != "Unselect") {
#     print("Variance Filter Activated")
#     
#     
#     MSET_VarFilter <- VarReductionFilter(obj, Labels)
#     obj <- MSET_VarFilter[, -ncol(MSET_VarFilter)]
#     print("var")
#     if (input$Norm == "Normal") {
#       data1 = data_fitting(as.data.frame(t(obj)), Labels)
#       object_Seurat <-
#         as.matrix(SummarizedExperiment::assay(data1, "counts"))
#       
#       
#       MetaData <- data.frame(groups = Labels)
#       
#       Input <-
#         Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = MetaData)
#       seuratObj <-
#         Seurat::NormalizeData(object  = Input, normalization.method = "LogNormalize")
#       dubstepR.out <-
#         DUBStepR::DUBStepR(
#           input.data = seuratObj@assays$RNA@data,
#           min.cells = 0.05 * ncol(Input),
#           optimise.features = T,
#           k = input$k,
#           num.pcs = input$np,
#           error = 0
#         )
#       seuratObj@assays$RNA@var.features <-
#         dubstepR.out$optimal.feature.genes
#       genes = seuratObj@assays$RNA@var.features
#       print("Norm")
#     }
#     else{
#       print("her is the problme")
#       data1 = data_fitting(as.data.frame(t(obj)), Labels)
#       object_Seurat <-
#         as.matrix(SummarizedExperiment::assay(data1, "counts"))
#       
#       MetaData <- data.frame(groups = Labels)
#       
#       Input <-
#         Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = MetaData)
#       
#      # Input@assays$RNA$data <- ScaleData(Input@assays$RNA$counts )
#   
#       # Input <- Seurat::NormalizeData(object = Input, normalization.method = "LogNormalize", scale.factor = 10000)
#       # # Normalize the data
#       # seuratObj <- NormalizeData(seuratObj)
#       # 
#       # # Then scale the data
#       # seuratObj <- ScaleData(seuratObj)
#       Input <-ScaleData(Input)
#       
#       dubstepR.out <-
#         DUBStepR::DUBStepR(
#           input.data = Input@assays$RNA$scale.data,
#           min.cells = 0.05 * ncol(Input),
#           optimise.features = T,
#           k = input$k,
#           num.pcs = input$np,
#           error = 0
#         )
#       
#       Input@assays$RNA@var.features <-
#         dubstepR.out$optimal.feature.genes
#       genes = Input@assays$RNA@var.features
#       print("no normal")
#     }
#   } else if (input$Norm == "Normal") {
#     print("no Var")
#     
#     data1 = data_fitting(as.data.frame(t(obj)), Labels)
#     object_Seurat <-
#       as.matrix(SummarizedExperiment::assay(data1, "counts"))
#     
#     
#     MetaData <- data.frame(groups = Labels)
#     
#     Input <-
#       Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = MetaData)
#     seuratObj <-
#       Seurat::NormalizeData(object  = Input, normalization.method = "LogNormalize")
#     dubstepR.out <-
#       DUBStepR::DUBStepR(
#         input.data = seuratObj@assays$RNA@data,
#         min.cells = 0.05 * ncol(Input),
#         optimise.features = T,
#         k = input$k,
#         num.pcs = input$np,
#         error = 0
#       )
#     seuratObj@assays$RNA@var.features <-
#       dubstepR.out$optimal.feature.genes
#     genes = seuratObj@assays$RNA@var.features
#     print("finish")
#   } else{
#     print("no var no normal")
#     
#     data1 = data_fitting(as.data.frame(t(obj)), Labels)
#     object_Seurat <-
#       as.matrix(SummarizedExperiment::assay(data1, "counts"))
#     
#     MetaData <- data.frame(groups = Labels)
#     Input <-
#       Seurat::CreateSeuratObject(counts = object_Seurat, meta.data = MetaData)
#     dubstepR.out <-
#       DUBStepR::DUBStepR(
#         input.data = Input@assays$RNA@data,
#         min.cells = 0.05 * ncol(Input),
#         optimise.features = T,
#         k = input$k,
#         num.pcs = input$np,
#         error = 0
#       )
#     Input@assays$RNA@var.features <-
#       dubstepR.out$optimal.feature.genes
#     genes = Input@assays$RNA@var.features
#     
#   }
#   
#   print("maybe here_____??")
#   
#  
#   
#   iG1 <- dubstepR.out[["corr.info"]]
#   iG = as.data.frame(iG1[, -1])
#   iG = as.data.frame(iG[order(iG[, 1], decreasing = TRUE), ])
#   
#   
#   iG <- data.frame(iG[which(iG1[, 1] %in% genes), ])
#   
#   rownames(iG) = genes
#   
#   iG <- iG
#   
#   newdata = data[, colnames(data) %in% genes]
#   newdata$Labels = as.factor(Labels)
#   return(list(ig=iG, newdata=newdata))
#   
# }



#Use the scPNMF feature selection method and return a dataframe the expression matrix of the isolated genes.

scPNMFfun = function(data, Labels, DM) {
  if (input$VarFilter != "Unselect") {
    print("Variance Filter Activated")
    MSET_VarFilter <- VarReductionFilter(data, Labels)
    obj <- MSET_VarFilter[, -ncol(MSET_VarFilter)]
    if (input$Norm == "Normal") {
      scale.factor <- mean(colSums(obj))

      obj <- Seurat::LogNormalize(obj ,
                                  scale.factor = scale.factor)

      obj=as.matrix(obj)
      obj=t(obj)
      print("i am here")

      res_pnmf <- scPNMF::PNMFfun(
        X = obj,
        K = 15,
        method = DM,
        tol = 1e-3,
        maxIter = 500,
        verboseN = TRUE
      )
      print("Norm and Var")
    }
    else{
      print("var no normal")

      obj = t(obj)

      res_pnmf <- scPNMF::PNMFfun(
        X = obj,
        K = 15,
        method = DM,
        tol = 1e-3,
        maxIter = 500,
        verboseN = TRUE
      )

    }

  }
  else{
    if (input$Norm == "Normal") {
      scale.factor <- mean(colSums(data[, -ncol(data)]))

      obj <- Seurat::LogNormalize(data[, -ncol(data)] ,
                                  scale.factor = scale.factor)
      obj = as.data.frame(obj)
      obj = t(obj)

      res_pnmf <- scPNMF::PNMFfun(
        X = obj,
        K = 15,
        method = DM,
        tol = 1e-3,
        maxIter = 500,
        verboseN = TRUE
      )
      print("Norm")
    }
    else{
      print(" no Var no normal")
      obj <- as.data.frame(data[, -ncol(data)])
      obj = t(obj)

      res_pnmf <- scPNMF::PNMFfun(
        X =obj,
        K = 15,
        method = DM,
        tol = 1e-3,
        maxIter = 500,
        verboseN = TRUE
      )

    }

  }

  W <- res_pnmf$Weight
  S <- res_pnmf$Score

  W_select <- scPNMF::basisSelect(
    W = W,
    S = S,
    X = obj,
    toTest = TRUE,
    toAnnotate = FALSE,
    mc.cores = 1
  )

  Score = as.data.frame(apply(W, 1, sum))
  W_select=as.data.frame(W_select)
  ig <-
    scPNMF::getInfoGene(
      W_select,
      M =  input$gM,
      by_basis = FALSE,
      return_trunW = TRUE,
      dim_use = NULL
    )
    iG = as.data.frame(Score[rownames(Score) %in%  ig[["InfoGene"]]  ,])
  rownames(iG) = ig[["InfoGene"]]
  iG <- iG[order(iG[, 1], decreasing = TRUE), , drop = FALSE]

  newdata = data[, colnames(data) %in%  ig[["InfoGene"]]]
  newdata$Labels = as.factor(Labels)
  return(list(ig=iG, newdata=newdata))
}


DUBStepRfun = function(data, Labels) {

  obj = as.matrix(data[, -ncol(data)])
  
  if (input$VarFilter != "Unselect") {
    print("Variance Filter Activated")
    
    
    MSET_VarFilter <- VarReductionFilter(obj, Labels)
    obj <- MSET_VarFilter[, -ncol(MSET_VarFilter)]
    print("var")
    if (input$Norm == "Normal") {
      data1 = data_fitting(as.data.frame(t(obj)), Labels)
      sce <- computeSumFactors(data1)
      sce <- logNormCounts(sce)
      dec <- modelGeneVar(sce)
      topHVGs <- getTopHVGs(dec, n = input$np)
      genes = topHVGs
      print("Norm")
    }
    else{
      print("her is the problme")
      
      data1 = data_fitting(as.data.frame(t(obj)), Labels)
      sce <- logNormCounts(data1)
      dec <- modelGeneVar(sce)
      topHVGs <- getTopHVGs(dec, n = input$np)
      genes = topHVGs
      print("no normal")
    }
  } else if (input$Norm == "Normal") {
    print("no Var")
    data1 = data_fitting(as.data.frame(t(obj)), Labels)
    sce <- computeSumFactors(data1)
    sce <- logNormCounts(sce)
    dec <- modelGeneVar(sce)
    topHVGs <- getTopHVGs(dec, n = input$np)
    genes = topHVGs
    print("Norm")
    print("finish")
  } else{
    print("no var no normal")
    
    data1 = data_fitting(as.data.frame(t(obj)), Labels)
    sce <- logNormCounts(data1)
    dec <- modelGeneVar(sce)
    topHVGs <- getTopHVGs(dec, n = input$np)
    genes = topHVGs
    print("no normal")
    
  }
  
  print("maybe here_____??")
  
  # Creating a dataframe with gene names and their scores
  Score <- data.frame(dec$bio)
  # df and Ordering the genes from highest to lowest score
  iG <- data.frame(row.names = rownames(dec), Score = Score[order(-dec$bio),] )
  
  # Filtering iG to include only rows with gene names that are in the 'genes' vector, keeping it as a data frame
  iG <- iG[which(rownames(iG) %in% genes), , drop = FALSE]
  
  newdata = as.data.frame(data[, colnames(data) %in% genes] )
  newdata$Labels = as.factor(Labels)
  return(list(ig=iG, newdata=newdata))
  
}
# Use the M3Drop feature selection method and return a dataframe the expression matrix of the isolated genes.

M3Dropfun = function(data, Labels) {
  if (input$VarFilter != "Unselect") {
    print("Variance Filter Activated")
    MSET_VarFilter <- VarReductionFilter(data, Labels)
    obj <- MSET_VarFilter[, -ncol(MSET_VarFilter)]
    if (input$Norm == "Normal") {
      print("eimai edw")
      
      scale.factor <- mean(colSums(obj))
      
      obj <- Seurat::LogNormalize(obj ,
                                  scale.factor = scale.factor)
      obj = as.data.frame(obj)
      obj = t(obj)
      #colnames(obj) = Labels
      DropsGenes <-
        M3DropFeatureSelection(
          obj,
          mt_method = input$M3Method,
          mt_threshold = input$M3dropthreshold,
          suppress.plot = TRUE
        )
      print("Norm and Var")
    } else{
      print("var no normal")
      obj = as.data.frame(t(obj))
      # obj=t(obj)
      
     # colnames(obj) = Labels
      DropsGenes <-
        M3DropFeatureSelection(
          obj,
          mt_method = input$M3Method,
          mt_threshold = input$M3dropthreshold,
          suppress.plot = TRUE
        )
      
    }
    
  }
  else{
    if (input$Norm == "Normal") {
      scale.factor <- mean(colSums(data[, -ncol(data)]))
      
      obj <- Seurat::LogNormalize(data[, -ncol(data)] ,
                                  scale.factor = scale.factor)
      obj = as.data.frame(obj)
      obj = t(obj)
     # colnames(obj) = Labels
      DropsGenes <-
        M3DropFeatureSelection(
          obj,
          mt_method = input$M3Method,
          mt_threshold = input$M3dropthreshold,
          suppress.plot = TRUE
        )
      
      print("Norm")
    } else{
      print(" no Var no normal")
      
      obj <- as.data.frame(data[, -ncol(data)])
      obj = t(obj)
     # colnames(obj) = Labels
      DropsGenes <-
        M3DropFeatureSelection(
          obj,
          mt_method = input$M3Method,
          mt_threshold = input$M3dropthreshold,
          suppress.plot = TRUE
        )
      
      
    }
    
  }
  
  
  iG <-
    DropsGenes[order(DropsGenes[, 4], decreasing = FALSE), 4 , drop = FALSE]
  
  newdata = data[, colnames(data) %in%  rownames(iG)]
  newdata$Labels = as.factor(Labels)
  newdata <- newdata
  return(list(ig=iG, newdata=newdata))
}
