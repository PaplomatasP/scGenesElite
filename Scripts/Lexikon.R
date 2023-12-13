#Translate the genes Id to genes Symbol

LexikonFun = function(data, organismus, id) {
  # 1. Convert from ensembl.gene to gene.symbol
  if (organismus == "Human") {
    if (id == "EnsemblGenes") {
      Genes <- c(colnames(data[, -ncol(data)]))
      geneIDS <-
        ensembldb::select(
          EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79,
          keys = Genes,
          keytype = "GENEID",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
      num = which(colnames(data) %in% geneIDS[, 2])
      newdata = data[, num]
      colnames(newdata) = unique(geneIDS[, 1])
      CleanCN = as.data.frame(colnames(newdata))
      CleanCN = na.omit(CleanCN)
      newdata = newdata[, colnames(newdata) %in% CleanCN[, 1]]
      newdata$Labels = data[, ncol(data)]
    }
    # 2. Convert from gene.symbol to ensembl.gene
    if (id == "ENTREZID") {
      Genes <- c(colnames(data[, -ncol(data)]))
      geneIDS <-
        ensembldb::select(
          EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79,
          keys = Genes,
          keytype = "ENTREZID",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
      num = which(colnames(data) %in% geneIDS[, 3])
      newdata = data[, num]
      colnames(newdata) = unique(geneIDS[, 1])
      CleanCN = as.data.frame(colnames(newdata))
      CleanCN = na.omit(CleanCN)
      newdata = newdata[, colnames(newdata) %in% CleanCN[, 1]]
      newdata$Labels = data[, ncol(data)]
      
    }
    if (id == "SYMBOL") {
      Genes <- c(colnames(data[, -ncol(data)]))
      geneIDS <-
        ensembldb::select(
          EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79,
          keys = Genes,
          keytype = "SYMBOL",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
      num = which(colnames(data) %in% geneIDS[, 1])
      newdata = data[, num]
      colnames(newdata) = unique(geneIDS[, 1])
      CleanCN = as.data.frame(colnames(newdata))
      CleanCN = na.omit(CleanCN)
      newdata = newdata[, colnames(newdata) %in% CleanCN[, 1]]
      
      newdata$Labels = data[, ncol(data)]
    }
    
    
  }
  
  if (organismus == "Mouse") {
    if (id == "EnsemblGenes") {
      Genes <- c(colnames(data[, -ncol(data)]))
      geneIDS <-
        ensembldb::select(
          EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
          keys = Genes,
          keytype = "GENEID",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
      num = which(colnames(data) %in% geneIDS[, 2])
      newdata = data[, num]
      colnames(newdata) = unique(geneIDS[, 1])
      CleanCN = as.data.frame(colnames(newdata))
      CleanCN = na.omit(CleanCN)
      newdata = newdata[, colnames(newdata) %in% CleanCN[, 1]]
      
      newdata$Labels = data[, ncol(data)]
    }
    # 2. Convert from gene.symbol to ensembl.gene
    if (id == "ENTREZID") {
      Genes <- c(colnames(data[, -ncol(data)]))
      geneIDS <-
        ensembldb::select(
          EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
          keys = Genes,
          keytype = "ENTREZID",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
      num = which(colnames(data) %in% geneIDS[, 3])
      newdata = data[, num]
      colnames(newdata) = unique(geneIDS[, 1])
      CleanCN = as.data.frame(colnames(newdata))
      CleanCN = na.omit(CleanCN)
      newdata = newdata[, colnames(newdata) %in% CleanCN[, 1]]
      
      newdata$Labels = data[, ncol(data)]
      
    }
    if (id == "SYMBOL") {
      Genes <- c(colnames(data[, -ncol(data)]))
      geneIDS <-
        ensembldb::select(
          EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
          keys = Genes,
          keytype = "SYMBOL",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
      num = which(colnames(data) %in% geneIDS[, 1])
      newdata = data[, num]
      colnames(newdata) = unique(geneIDS[, 1])
      CleanCN = as.data.frame(colnames(newdata))
      CleanCN = na.omit(CleanCN)
      newdata = newdata[, colnames(newdata) %in% CleanCN[, 1]]
      newdata$Labels = data[, ncol(data)]
    }
    
    
  }
  return(newdata)
}

# iGlexikon = function(iG, id,organismus) {
#   if (Species == "hsa"){
#     Genes <- c(colnames(data[, -ncol(data)]))
#     Genes
#     
#   }
#   Genes <- c(colnames(as.data.frame(t(iG))))
#   geneIDS <-
#     ensembldb::select(
#       EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
#       keys = Genes,
#       keytype = "SYMBOL",
#       columns = c("GENENAME", "GENEID", "ENTREZID")
#     )
#   if (nrow(geneIDS[-which(is.na(geneIDS[, 4][1:nrow(geneIDS)])), ]) >
#       0) {
#     geneIDS = geneIDS[-which(is.na(geneIDS[, 4][1:nrow(geneIDS)])), ]
#   }
#   geneIDS = geneIDS[!duplicated(geneIDS[, 1]),]
#   num = which(colnames(as.data.frame(t(iG))) %in% geneIDS[, 1])
#   
#   newdata = as.data.frame(t(iG)[, num])
#   newdata = as.data.frame(t(newdata))
#   colnames(newdata) = geneIDS[, 4]
#   newdata = as.data.frame(t(newdata))
#   Scaler <-
#     caret::preProcess(newdata , rangeBounds = c(-1, 1), method = "range")
#   ScaledLogiG <- predict(Scaler, newdata)
#   
#   #}
#   return(ScaledLogiG)
# }


iGlexikon = function(iG, id,organismus) {
  if (organismus == "Human") {
    if (id == "EnsemblGenes") {
      Genes <-rownames(iG)
      geneIDS <-
        ensembldb::select(
          EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79,
          keys = Genes,
          keytype = "GENEID",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
    }
    # 2. Convert from gene.symbol to ensembl.gene
    if (id == "ENTREZID") {
      Genes <- rownames(iG)
      geneIDS <-
        ensembldb::select(
          EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79,
          keys = Genes,
          keytype = "ENTREZID",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
    }
    if (id == "SYMBOL") {
      Genes <- rownames(iG)
      geneIDS <-
        ensembldb::select(
          EnsDb.Hsapiens.v79::EnsDb.Hsapiens.v79,
          keys = Genes,
          keytype = "SYMBOL",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )  
    }
  }
  if (organismus == "Mouse") {
    if (id == "EnsemblGenes") {
      Genes <- rownames(iG)
      geneIDS <-
        ensembldb::select(
          EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
          keys = Genes,
          keytype = "GENEID",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
    }
    # 2. Convert from gene.symbol to ensembl.gene
    if (id == "ENTREZID") {
      Genes <- rownames(iG)
      geneIDS <-
        ensembldb::select(
          EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
          keys = Genes,
          keytype = "ENTREZID",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
  
    }
    if (id == "SYMBOL") {
      Genes <- rownames(iG)
      geneIDS <-
        ensembldb::select(
          EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
          keys = Genes,
          keytype = "SYMBOL",
          columns = c("GENENAME", "GENEID", "ENTREZID")
        )
    }
  }
  if (nrow(geneIDS[-which(is.na(geneIDS[, 4][1:nrow(geneIDS)])), ]) >
      0) {
    geneIDS = geneIDS[-which(is.na(geneIDS[, 4][1:nrow(geneIDS)])), ]
  }
  geneIDS = geneIDS[!duplicated(geneIDS[, 1]),]
  num = which(colnames(as.data.frame(t(iG))) %in% geneIDS[, 1])
  
  newdata = as.data.frame(t(iG)[, num])
  newdata = as.data.frame(t(newdata))
  colnames(newdata) = geneIDS[, 4]
  newdata = as.data.frame(t(newdata))
  Scaler <-
    caret::preProcess(newdata , rangeBounds = c(-1, 1), method = "range")
  ScaledLogiG <- predict(Scaler, newdata)
  
  #}
  return(ScaledLogiG)
}