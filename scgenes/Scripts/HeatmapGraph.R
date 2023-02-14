

#Create a complexHeatmap of the isolated genes and plot it.
complexHeatMapFun = function(data,Plot=TRUE) {
  withProgress(message = 'Please wait........', value = 0, {
    {
      incProgress(2 / 10)
      Sys.sleep(0.10)

      Labels = as.factor(data[, ncol(data)])
      data <- data[, -ncol(data)]
      print(dim(data))
      #Cells=rownames(data)
    }

    {
      incProgress(6 / 10)
      Sys.sleep(0.10)
      if (input$organismus == "Human") {
        print("Human")
        ref <- celldex::HumanPrimaryCellAtlasData()
        
        data = as.data.frame(t(data))
        results <-
          SingleR(test = data,
                  ref = ref,
                  labels = ref$label.main)
        data = as.data.frame(t(data))
        data$singlr_labels <- results$labels
        
        
        LabelsCellType = as.factor(data$singlr_labels)
        data = aggregate(data[, 1:length(data[, -ncol(data)])], list(data$singlr_labels), FUN =
                           mean)
        
        data = as.matrix(t(data))
        colnames(data) = data[1, ]
        data = data[-1, ]
        data <-
          data.frame(apply(data, 2, function(x)
            as.numeric(as.character(x))))
        rownames(data) = rownames(iG)
        
      } else{
        ref <- celldex::MouseRNAseqData()
        print("Mouse")
        data = as.data.frame(t(data))
        results <-
          SingleR(test = data,
                  ref = ref,
                  labels = ref$label.main)
        data = as.data.frame(t(data))
        
        if (input$Split == "CellType") {
          data$singlr_labels <- results$labels
          
          
          HeatmapLabel = as.factor(data$singlr_labels)
          data = aggregate(data[, 1:length(data[, -ncol(data)])], list(data$singlr_labels), FUN =
                             mean)
          
        } else{
          HeatmapLabel = as.factor(Labels)
          data$singlr_labels <- HeatmapLabel
          data = aggregate(data[, 1:length(data[, -ncol(data)])], list(data$singlr_labels), FUN =
                             mean)
        }
        
        
        data <- as.matrix(t(data))
        colnames(data) = data[1, ]
        data = data[-1, ]
        genes=rownames(data)
        data <-
          data.frame(apply(data, 2, function(x)
            as.numeric(as.character(x))))
        
        rownames(data) = genes  #rownames(iG)
        data=data[order(rowMeans(data) ,decreasing = TRUE),]
        
      }
      
      
      
      
      my_palette <-
        colorRampPalette(c(
          "darkred",
                   "red",
                   "#222222",
                   "blue4",
                   "orange",
                   "purple",
                   "#296e28"
        ))(n = length(levels(HeatmapLabel)))
      
      column_ha <- HeatmapAnnotation(cluster = anno_block(
        #
        gp = gpar(fill = my_palette[1:length(levels(HeatmapLabel))]),
        # <- here controls the filled color
        labels = levels(HeatmapLabel),
        labels_gp = gpar(col = "white", fontsize = 10)
      ))
      
      if (input$genes < nrow(iG)) {
        GenesN = input$genes
        
      } else{
        GenesN = nrow(iG)
      }
     
      
      if (input$clustering == "cluster_between_groups") {
        dend1 = cluster_between_groups(data[1:input$genes, ], levels(HeatmapLabel))
        
        print("cluster_between_groups")
        
        
        
        
        
        HeatMAP = Heatmap(
          as.matrix(data[1:GenesN, ]),
          cluster_columns = dend1,
          cluster_rows =FALSE,
          column_split = length(levels(HeatmapLabel)),
          row_labels = rownames(data)[1:GenesN],
          top_annotation = column_ha,
          show_column_dend = FALSE,
          row_title = "",
          row_names_gp = gpar(fontsize = 5),
          show_row_dend = FALSE,
          column_names_gp = gpar(fontsize = 2.5),
          show_row_names = TRUE,
          show_column_names = FALSE,
          use_raster = FALSE,
          raster_by_magick = FALSE,
          column_title = NULL,
          heatmap_legend_param = list(title = "Expression Matrix")
        )
        
        
      } else{
        print("None Clusterin")
       # library(scales)
        
        # data$CK <- rescale(data$CK, to = c(1, 2))
        # data$CKp25<- rescale(data$CKp25, to = c(1, 2))
        
        HeatMAP = Heatmap(
          as.matrix(data[1:GenesN, ]),
          top_annotation = column_ha,
          cluster_columns = FALSE,
          cluster_rows =FALSE,
          row_labels = rownames(data)[1:GenesN],
          show_column_dend = FALSE,
          column_split = levels(HeatmapLabel) ,
          row_title = "",
          row_names_gp = gpar(fontsize = 5),
          show_row_dend = FALSE,
          column_names_gp = gpar(fontsize = 2.5),
          show_row_names = TRUE,
          show_column_names = FALSE,
          use_raster = FALSE,
          raster_by_magick = FALSE,
          column_title = NULL,
          heatmap_legend_param = list(title = "Expression Matrix")
        )
        
        
      }
      if(Plot==TRUE){ 
      draw(HeatMAP, heatmap_legend_side = "left")
      }
    }
  })
  s1 = lapply(GenesN,
              function(i)
                data[1:i,])
  s1 = as.data.frame(s1)
  S = as.data.frame(row.names(s1))
  colnames(S) = paste("Heatmap Genes:")
  
  return(S)
}
#Create a Similarity graph and plot it.

GraphsFun = function(data) {
  withProgress(message = 'Please wait........', value = 0, {
    {
      incProgress(2 / 10)
      Sys.sleep(0.10)
      
      Labels = data[, ncol(data)]
      data = data[, -ncol(data)]
      if (input$Genes < nrow(iG)) {
        GenesN = input$Genes
        
      } else{
        GenesN = nrow(iG)
      }
      data=data[colnames(data) %in% rownames(iG)[1:GenesN]]
      newdata = as.matrix(data)
      
      
      # Create a graph adjacency based on correlation distances between genes in  pairwise fashion.
      g <- graph.adjacency(
        as.matrix(as.dist(
          cor(newdata, use = "pairwise.complete.obs", method = "pearson")
        )),
        mode = "plus",
        weighted = TRUE,
        diag = FALSE
      )
    }
    
    {
      incProgress(6 / 10)
      Sys.sleep(0.10)
      # Simplfy the adjacency object
      g <- simplify(g)
      E(g)$weight
      # E(g)[is.na(E(g)$weight)]<-0.0001
      # subset( E(g), (E(g)$weight=="NA"))
      # E(g)$weight=na.omit(E(g)$weight)
      #
      # E(g)[is.na(E(g)$weight)] <- 1
      # Colour negative correlation edges as blue
      E(g)[which(E(g)$weight < 0)]$color <- "darkblue"
        
      # Colour positive correlation edges as red
      E(g)[which(E(g)$weight > 0)]$color <- "darkred"
        
      # Convert edge weights to absolute values
      E(g)$weight <- abs(E(g)$weight)
      
      # Change arrow size
      # For directed graphs only
      #E(g)$arrow.size <- 1.0
      
      # Remove edges below absolute Pearson correlation 0.8
      g <-
        delete_edges(g, E(g)[which(E(g)$weight < input$Pearson_correlation)])
      
      
      
      # Assign names to the graph vertices (optional)
      V(g)$name <- V(g)$name
      
      # Change shape of graph vertices
      V(g)$shape <- "sphere"
      
      # Change colour of graph vertices
      V(g)$color <- "skyblue"
        
      # Change colour of vertex frames
      V(g)$vertex.frame.color <- "white"
        
      # Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
      # Multiply scaled vales by a factor of 10
      scale01 <- function(x) {
        (x - min(x)) / (max(x) - min(x))
      }
      vSizes <- (scale01(apply(t(newdata), 1, mean)) + 1.0) * 10
      
      # Amplify or decrease the width of the edges
      edgeweights <- E(g)$weight * 2.0
      
      # Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
      mst <- mst(g, algorithm = "prim")
      
      E(g)$weight
      
      mst.communities <-
        edge.betweenness.community(mst, weights = NULL, directed = FALSE)
      mst.clustering <-
        make_clusters(mst, membership = mst.communities$membership)
      V(mst)$color <- mst.communities$membership + 1
      
      
      visIgraph(mst,
                idToLabel = TRUE,
                layout = "layout_nicely",
                type = "full") %>%
        visIgraphLayout() %>%
        visNodes(
          shape = "dot",
          color = list(
            background = "#0085AF",
            border = "#013848",
            highlight = "#FF8000"
          ),
          shadow = list(
            enabled = TRUE,
            width = 12,
            size = 12
          )
        ) %>%
        visEdges(
          width = 12,
          shadow = TRUE,
          color = list(color = "#d12d71", highlight = "red")
        ) %>%
        visOptions(
          width = NULL,
          height = NULL,
          highlightNearest = TRUE,
          nodesIdSelection = TRUE,
          selectedBy = NULL,
          collapse = FALSE,
          autoResize = TRUE,
          clickToUse = TRUE,
          manipulation = NULL
        )  %>%
        visLayout(randomSeed = 11)
    }
  })
  
}

#function to create a PPI network analysis to the isolated genes.
PPInetwork = function(Genes) {
  withProgress(message = 'Please wait........', value = 0, {
    {
      incProgress(2 / 10)
      Sys.sleep(0.10)
      
      
      if (input$organismus == "Human") {
        PPI_species = 9606
      } else{
        PPI_species = 10090
      }
      
      string_db <- STRINGdb$new(
        version = "11.5",
        species = PPI_species,
        score_threshold = input$Score_Threshold_PPI,
        network_type = "full",
        input_directory = ""
      )
    }
    {
      incProgress(8 / 10)
      Sys.sleep(0.10)
      if (input$Genes < nrow(iG)) {
        GenesN = input$Genes
        
      } else{
        GenesN = nrow(iG)
      }
      
      dfGenes<- data.frame(cbind(Genes[,1][1:GenesN],row.names(Genes)[1:GenesN]) )
   
      colnames(dfGenes)[1] = "Score"
      colnames(dfGenes)[2] = "Genes"
      
      Mapped <-
        string_db$map(dfGenes, "Genes", removeUnmappedRows = TRUE)
      hits <- Mapped$STRING_id
      string_db$plot_network(hits)
      
    }
  })
}

