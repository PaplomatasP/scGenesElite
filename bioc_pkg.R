InstallBioc=function(){

  install.packages("BiocManager")

  BiocManagerPackages=c("M3Drop","ComplexHeatmap","pathview","ensembldb","celldex","org.Mm.eg.db",
                        "AnnotationFilter","AnnotationDbi","twoddpcr","EnsDb.Mmusculus.v79","EnsDb.Hsapiens.v79",
                        "Biobase","BiocFileCache","BiocGenerics","BiocParallel","BiocStyle",
                        "SingleR","STRINGdb","MAST" ,  "DESeq2")


  for ( Biopackage in BiocManagerPackages){
    BiocManager::install(Biopackage)
  }
}

InstallBioc()


