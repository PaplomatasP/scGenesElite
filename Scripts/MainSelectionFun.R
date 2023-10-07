                                #Here is the main function that performs the feature selection method and visualizes the isolated genes.

MethodData = function() {
  #check if the input data are rds. and do the all process ........
  if (length(input$file1) == 0) {

    
    
    if (input$GENEid == "EnsemblGenes") {
    
      RDS_file <- readRDS(input$rdsFile$datapath)
      # Convert to title case
    #  colnames(RDS_file) <- tools::toTitleCase(tolower(colnames(RDS_file)))
      RDS_file1 <- LexikonFun(RDS_file, input$organismus, input$GENEid)
  
       
      validate(need(
        ncol(RDS_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
      
    }
    if (input$GENEid == "ENTREZID") {
    
      RDS_file <- readRDS(input$rdsFile$datapath)
     # colnames(RDS_file) <- tools::toTitleCase(tolower(colnames(RDS_file)))
      
      RDS_file1 <- LexikonFun(RDS_file, input$organismus, input$GENEid)
      validate(need(
        ncol(RDS_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
      
    }
    if (input$GENEid == "SYMBOL") {
      
      RDS_file <- readRDS(input$rdsFile$datapath)
    #  colnames(RDS_file) <- tools::toTitleCase(tolower(colnames(RDS_file)))
      
      RDS_file1 <- LexikonFun(RDS_file, input$organismus, input$GENEid)
      validate(need(
        ncol(RDS_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
    }
    

    withProgress(message = 'Please wait........', value = 0
                 , {
                   {
                     incProgress(3 / 10)
                     Sys.sleep(0.10)
                   }
                   if (input$P_method == "Empty" &
                       input$VariableM == "NoMethod"  &
                       input$ML_Method == "Empty") {
                     print("EnsemleMethod")
                     FilterData <- EnsemleMethod(obj = RDS_file1[,-ncol(RDS_file1)],
                                                  EnseLabels = RDS_file1[, ncol(RDS_file1)])
                   }
                   
                   
                   else if (input$P_method == "Empty" &
                          input$VariableM == "NoMethod"   &
                          input$ensembleVar=="NoMethod" & 
                           input$ensemblePvalue=="NoMethod"
                         ) {
                     print("ML_Method")
                     FilterData <- SelectionFilter1(
                       data = RDS_file1[,-ncol(RDS_file1)],
                       Labels = RDS_file1[, ncol(RDS_file1)],
                       MLmethod = input$ML_Method
                       
                     )
                   }
                   
                   else if (input$VariableM == "SCMarker" &  input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
                            & input$ensembleWrapper=="NoMethod" & 
                            input$ML_Method=="Empty") {
                     colnames(RDS_file1) <- tools::toTitleCase(tolower(colnames(RDS_file1)))
                     FilterData <- SCMarkerfun(
                       RDS_file1,
                       Labels = RDS_file1[, ncol(RDS_file1)],
                       GeneSK = input$geneK,
                       CellSK = input$cellK
                       
                     )
                    
                   }
                   
                   else if (input$VariableM == "SelfE" &  input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
                            & input$ensembleWrapper=="NoMethod" & 
                             input$ML_Method=="Empty") {
                     FilterData <- SelfEGenes(RDS_file1,
                                              Labels = RDS_file1[, ncol(RDS_file1)],
                                              n = input$n)
                     
                     
                   }
                   else if (input$VariableM == "DUBStepR" &  input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
                            & input$ensembleWrapper=="NoMethod" & 
                          input$ML_Method=="Empty") {
                     
                     FilterData <- DUBStepRfun(RDS_file1,
                                               Labels = RDS_file1[, ncol(RDS_file1)])
                     
                     
                   }
                   
                   else if (input$VariableM == "ScPNMF" &  input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
                            & input$ensembleWrapper=="NoMethod" & 
                            input$ML_Method=="Empty") {
                     FilterData <- scPNMFfun(RDS_file1,
                                             Labels = RDS_file1[, ncol(RDS_file1)],
                                             DM <- input$distMethod)
                     
                     
                   }
                   else if (input$VariableM == "M3Drop"  &  input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
                            & input$ensembleWrapper=="NoMethod" & 
                             input$ML_Method=="Empty") {
                     
                     FilterData <- M3Dropfun(RDS_file1,
                                             Labels = RDS_file1[, ncol(RDS_file1)])
                   }
                   else if  (input$VariableM == "NoMethod"  & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
                        & input$ensembleWrapper=="NoMethod" & input$P_method!="Empty"& 
                         input$ML_Method=="Empty"){
                     FilterData <- SelectionFilter( 
                       data = RDS_file1[,-ncol(RDS_file1)],
                       Labels = RDS_file1[, ncol(RDS_file1)],
                       PvalueNum = input$PvalueNum , logfc= input$logfc
                     )
                   }else {
                    shinyalert(title = "Message",type = "error",
                               text = "Please choose the appropriate options. It seems that either data was not uploaded or the uploaded data is incorrect, or the necessary conditions were not met.", closeOnClickOutside = TRUE) 
                   }
                   { 
                     incProgress(7 / 10)
                     Sys.sleep(0.20)
                   }
                   
                   tryCatch({
                    
                     iG <- FilterData$ig
                    if (exists("iG")) {
                   output$TheBarPlot <- renderPlot(execOnResize = FALSE,{
                     
                         
                    
                     dfbar = as.data.frame(head(iG, input$genes) ) #input$genes
                    
                     ColorFun <-
                       colorRampPalette(c("#CCCCCC" , "#104E8B"))
                     ColorPaleta <- ColorFun(n = nrow(x = dfbar))
                     
                     dfbar$Color <-
                       as.character(x = cut(
                         x = rank(x = dfbar[, 1])  # used to assign order in the event of ties
                         ,
                         breaks = nrow(x = dfbar)  # same as the 'n' supplied in ColorFun
                         ,
                         labels = ColorPaleta  # label the groups with the color in ColorPaleta
                       ))
                   
                     
                     par(mar = c(7, 4.2, 4.1, 3))
                     barplot(
                       height = dfbar[, 1],
                       names.arg = rownames(dfbar),
                       las = 2,
                       col = dfbar$Color,
                       border = NA,
                       main = "Most Important Genes which Act as Potential Biomarkers for the given case-study",
                       cex.main = 1.2,
                       #xlab = "Genes ID",
                       cex.names = 0.5
                     )
                     mtext(
                       "Genes ID",
                       side = 1,
                       line = 3,
                       cex = 1.2,
                       font = 2,
                       col = "black",
                       family = "Calibri Light",
                       padj = 1.5
                     )
                     
                     
                     
                   }) }else {
                     showModal(modalDialog(
                       title = "Error",
                       div(style = "color: red;", "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more methods can be select.")
                     ))
                     
                    
                     
                   }
                   }, error = function(e) {
                     print("")
                   })
                   
                   tryCatch({
                     if  (exists("iG") ) {
                       
                   output$GenesList = renderTable({
                     s1 = lapply(input$genes,
                                 function(i)
                                   iG[1:i,])
                     s1 = as.data.frame(s1)
                     rownames(s1) = make.names(row.names(iG)[1:nrow(s1)], unique = TRUE)
                     
                     
                     colnames(s1) = paste("# of Genes which operate as Biomarkers: ", nrow(iG))
                     s1[, 1] = rownames(s1)
                     
                     s1
                   }
                   , rownames = FALSE)} else {
                     showModal(modalDialog(
                       title = "Message",
                       "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more can be select."
                     ))
                   }
                   }, error = function(e) {
                     print("")
                   })
                   
                   
                   tryCatch({
                     if (exists("FilterData") ) {
                       
                       iG <- FilterData$ig
                       NewData = FilterData$newdata
                     
                       
                   output$HeatmapList = renderTable({
                     if (input$HeatMap1 == TRUE) {
                       complexHeatMapFun(NewData,iG,Plot=FALSE)
                      
                    
                     }
                     
                   })
                   
                  
                   output$KnnClassifier = renderPlot({
                     KnnClassifier(data = NewData,iG, Labels = NewData[, ncol(NewData)])
                     
                   })
                   output$HeatMap <-
                     renderPlot(execOnResize = FALSE, {
                       if (input$HeatMap1 == TRUE) {
                         complexHeatMapFun(NewData,iG,Plot=TRUE)
                       }
                       
                     })
                   
                   }else {
                       showModal(modalDialog(
                         title = "Message",
                         "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more can be select."
                       ))
                     }
                   }, error = function(e) {
                     print("")
                   })
                   
                   
                   
                   output$downloadData <- downloadHandler(
                     filename = function() {
                       paste("FilterData-", Sys.Date(), ".csv", sep = "")
                     },
                     content = function(file) {
                       write.csv(NewData, file)
                     }
                   )
                   {
                     incProgress(10 / 10)
                     Sys.sleep(0.45)
                   }
                   
                 })
    
    showModal(
      modalDialog(
        title = "The filtered data can be downloaded!",
        paste0("Click on the Download Button"),
        easyClose = TRUE,
        footer = NULL
      )
    )
    
  }
  else {
    if (input$GENEid == "EnsemblGenes") {
      CSV_file <- read.csv(input$file1$datapath)
     # colnames(CSV_file) <- tools::toTitleCase(tolower(colnames(CSV_file)))
      
      CSV_file1 <- LexikonFun(CSV_file, input$organismus, input$GENEid)
      rownames(CSV_file1)=rownames(CSV_file)
      validate(need(
        ncol(CSV_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
      
    }
    if (input$GENEid == "ENTREZID") {
      CSV_file <- read.csv(input$file1$datapath)
      #colnames(CSV_file) <- tools::toTitleCase(tolower(colnames(CSV_file)))
      CSV_file1 <- LexikonFun(CSV_file, input$organismus, input$GENEid)
      rownames(CSV_file1)=rownames(CSV_file)
      validate(need(
        ncol(CSV_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
      
    }
    if (input$GENEid == "SYMBOL") {
      CSV_file <- read.csv(input$file1$datapath)
    #  colnames(CSV_file) <- tools::toTitleCase(tolower(colnames(CSV_file)))
      CSV_file1 <- LexikonFun(CSV_file, input$organismus, input$GENEid)
      rownames(CSV_file1)=rownames(CSV_file)
      validate(need(
        ncol(CSV_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
    }
    withProgress(message = 'Please wait........', value = 0, {
      # Number of times we'll go through the loop
      {
        incProgress(4 / 10)
        Sys.sleep(0.10)
      }
      
      
      
      
      
      req(input$click) #to prevent print at first lauch
      if (input$P_method == "Empty" &
          input$VariableM == "NoMethod"  &
          input$ML_Method == "Empty") {
        print("EnsemleMethod")
        FilterData <- EnsemleMethod(obj = CSV_file1[,-ncol(CSV_file1)],
                                    EnseLabels = CSV_file1[, ncol(CSV_file1)])
      }
      
      
      else if (input$P_method == "Empty" &
               input$VariableM == "NoMethod"   &
               input$ensembleVar=="NoMethod" & 
               input$ensemblePvalue=="NoMethod") {
        print("ML_Method")
        FilterData <- SelectionFilter1(
          data = CSV_file1[,-ncol(CSV_file1)],
          Labels = CSV_file1[, ncol(CSV_file1)],
          MLmethod = input$ML_Method

        )
      }
      
      else if (input$VariableM == "SCMarker" &  input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
               & input$ensembleWrapper=="NoMethod" & 
               input$ML_Method=="Empty") {
        FilterData <- SCMarkerfun(
          CSV_file1,
          Labels = CSV_file1[, ncol(CSV_file1)],
          GeneSK = input$geneK,
          CellSK = input$cellK
        )
        
      }
      
      else if (input$VariableM == "SelfE" &  input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
               & input$ensembleWrapper=="NoMethod" & 
               input$ML_Method=="Empty") {
        FilterData <- SelfEGenes(CSV_file1,
                                 Labels = CSV_file1[, ncol(CSV_file1)],
                                 n = input$n)
        
        
      }
      else if (input$VariableM == "DUBStepR" & input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
               & input$ensembleWrapper=="NoMethod" & 
               input$ML_Method=="Empty") {
        
        FilterData <- DUBStepRfun(CSV_file1,
                                  Labels = CSV_file1[, ncol(CSV_file1)])
        
        
      }
      
      else if (input$VariableM == "ScPNMF" & input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
               & input$ensembleWrapper=="NoMethod" & 
               input$ML_Method=="Empty") {
        FilterData <- scPNMFfun(CSV_file1,
                                Labels = CSV_file1[, ncol(CSV_file1)],
                                DM <- input$distMethod)
        
        
      }
      else if (input$VariableM == "M3Drop"  & input$P_method=="Empty" & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
               & input$ensembleWrapper=="NoMethod" & 
               input$ML_Method=="Empty") {
        
        FilterData <- M3Dropfun(CSV_file1,
                                Labels = CSV_file1[, ncol(CSV_file1)])
      }
      else if  (input$VariableM == "NoMethod"  & input$ensembleVar=="NoMethod" & input$ensemblePvalue=="NoMethod"
                & input$ensembleWrapper=="NoMethod" & input$P_method!="Empty"& 
                input$ML_Method=="Empty"){
        FilterData <- SelectionFilter( 
          data = CSV_file1[,-ncol(CSV_file1)],
          Labels = CSV_file1[, ncol(CSV_file1)],
          PvalueNum = input$PvalueNum, logfc= input$logfc
        )
      }else {
        shinyalert(title = "Message",type = "error",
                   text = "Please choose the appropriate options. It seems that either data was not uploaded or the uploaded data is incorrect, or the necessary conditions were not met.", closeOnClickOutside = TRUE) 
      }
      { 
        incProgress(7 / 10)
        Sys.sleep(0.20)
      }
      
      tryCatch({
        iG <- FilterData$ig
        if (exists("iG")) {
          output$TheBarPlot <- renderPlot({
            
            
            
            dfbar = as.data.frame(head(iG, input$genes) )
            ColorFun <-
              colorRampPalette(c("#CCCCCC" , "#104E8B"))
            ColorPaleta <- ColorFun(n = nrow(x = dfbar))
            
            dfbar$Color <-
              as.character(x = cut(
                x = rank(x = dfbar[, 1])  # used to assign order in the event of ties
                ,
                breaks = nrow(x = dfbar)  # same as the 'n' supplied in ColorFun
                ,
                labels = ColorPaleta  # label the groups with the color in ColorPaleta
              ))
            
            par(mar = c(7, 4.2, 4.1, 3))
            barplot(
              height = dfbar[, 1],
              names.arg = rownames(dfbar),
              las = 2,
              col = dfbar$Color,
              border = NA,
              main = "Most Important Genes which Act as Potential Biomarkers for the given case-study",
              cex.main = 1.2,
              #xlab = "Genes ID",
              cex.names = 0.5
            )
            mtext(
              "Genes ID",
              side = 1,
              line = 3,
              cex = 1.2,
              font = 2,
              col = "black",
              family = "Calibri Light",
              padj = 1.5
            )
            
            
            
          }) }else {
            showModal(modalDialog(
              title = "Error",
              div(style = "color: red;", "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more methods can be select.")
            ))
            
            
            
          }
      }, error = function(e) {
        print("")
      })
      
      tryCatch({
        if  (exists("iG") ) {
          
          output$GenesList = renderTable({
            s1 = lapply(input$genes,
                        function(i)
                          iG[1:i,])
            s1 = as.data.frame(s1)
            rownames(s1) = make.names(row.names(iG)[1:nrow(s1)], unique = TRUE)
            
            
            colnames(s1) = paste("# of Genes which operate as Biomarkers: ", nrow(iG))
            s1[, 1] = rownames(s1)
            
            s1
          }
          , rownames = FALSE)} else {
            showModal(modalDialog(
              title = "Message",
              "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more can be select."
            ))
          }
      }, error = function(e) {
        print("sdsd")
      })
      
      
      tryCatch({
        if (exists("FilterData") ) {
          iG <- FilterData$ig
          NewData = FilterData$newdata
          
          output$HeatmapList = renderTable({
            if (input$HeatMap1 == TRUE) {
              complexHeatMapFun(NewData,iG,Plot=FALSE)
              
              
            }
            
          })
          
          output$KnnClassifier = renderPlot({
            KnnClassifier(data = NewData,iG, Labels = NewData[, ncol(NewData)])
            
          })
          output$HeatMap <-
            renderPlot(execOnResize = FALSE, {
              if (input$HeatMap1 == TRUE) {
                complexHeatMapFun(NewData,iG,Plot=TRUE)
                
              }
            })
          
        }else {
          showModal(modalDialog(
            title = "Message",
            "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more can be select."
          ))
        }
      }, error = function(e) {
        print("")
      })
      
      output$downloadData <- downloadHandler(
        filename = function() {
          paste("FilterData-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(NewData, file)
        }
      )
      {
        incProgress(10 / 10)
        Sys.sleep(0.30)
      }
      
    })
    showModal(
      modalDialog(
        title = "The filtered data can be downloaded!",
        paste0("Click on the Download Button"),
        easyClose = TRUE,
        footer = NULL
      )
    )
  }
  tryCatch({
    if  (exists("FilterData")) {
      FilterData=FilterData
    } else {
      FilterData="The are not available data!!!"
      showModal(modalDialog(
        title = "Message",
        "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more methods can be select.",
        easyClose = TRUE,footer = modalButton("OK")
      ))
    }
  }, error = function(e) {
   
  })
  
  return(FilterData)
  
}
