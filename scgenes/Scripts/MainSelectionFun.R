#Here is the main function that performs the feature selection method and visualizes the isolated genes.

MethodData = function() {
  has_csv <- !is.null(input$file1) && length(input$file1) > 0L
  has_rds <- !is.null(input$rdsFile) && length(input$rdsFile) > 0L
  validate(need(
    xor(has_csv, has_rds),
    "Upload exactly one dataset: either CSV or RDS, but not both."
  ))

  #check if the input data are rds. and do the all process ........
  if (!has_csv) {
    RDS_file <- read_uploaded_rds(input$rdsFile)

    
    
    if (input$GENEid == "EnsemblGenes") {
    
      # Convert to title case
    #  colnames(RDS_file) <- tools::toTitleCase(tolower(colnames(RDS_file)))
      RDS_file1 <- LexikonFun(RDS_file, input$organismus, input$GENEid)
  
       
      validate(need(
        ncol(RDS_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
      
    }
    if (input$GENEid == "ENTREZID") {
    
     # colnames(RDS_file) <- tools::toTitleCase(tolower(colnames(RDS_file)))
      
      RDS_file1 <- LexikonFun(RDS_file, input$organismus, input$GENEid)
      validate(need(
        ncol(RDS_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
      
    }
    if (input$GENEid == "SYMBOL") {
      
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
                       input$ML_Method == "Empty" &&
                       (is.null(input$SHAP_Method) || input$SHAP_Method == "Empty")) {
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
                   }
                   else if (input$VariableM == "NoMethod" & input$P_method == "Empty" & 
                           input$ensembleVar == "NoMethod" & input$ensemblePvalue == "NoMethod" &
                           input$ensembleWrapper == "NoMethod" & 
                           ((!is.null(input$SHAP_Method) && input$SHAP_Method != "Empty" && 
                            (input$SHAP_Method == "shap_rf" | input$SHAP_Method == "shap_xgb")) ||
                            (input$ML_Method == "shap_rf" | input$ML_Method == "shap_xgb"))) {
                     
                     # Handle SHAP Values methods - check both ML_Method and SHAP_Method
                     # Determine which method to use
                     shap_method <- NULL
                     if (!is.null(input$SHAP_Method) && input$SHAP_Method != "Empty" && 
                         (input$SHAP_Method == "shap_rf" || input$SHAP_Method == "shap_xgb")) {
                       shap_method <- input$SHAP_Method
                     } else if (input$ML_Method == "shap_rf" || input$ML_Method == "shap_xgb") {
                       shap_method <- input$ML_Method
                     }
                     
                     # Get importance limit
                     importance_limit <- input$importanceLimit
                     if (!is.null(input$SHAP_importanceLimit) && input$SHAP_importanceLimit != input$importanceLimit) {
                       importance_limit <- input$SHAP_importanceLimit
                     }
                     
                     # Call appropriate SHAP method
                     if (!is.null(shap_method)) {
                       if (shap_method == "shap_rf") {
                         FilterData <- ShapValuesFilter(
                           data = RDS_file1[,-ncol(RDS_file1)],
                           Labels = RDS_file1[, ncol(RDS_file1)],
                           MLmethod = "rf",
                           importanceLimit = importance_limit
                         )
                       } else if (shap_method == "shap_xgb") {
                         FilterData <- ShapValuesFilter(
                           data = RDS_file1[,-ncol(RDS_file1)],
                           Labels = RDS_file1[, ncol(RDS_file1)],
                           MLmethod = "xgbTree",
                           importanceLimit = importance_limit
                         )
                       }
                     }
                     
                     # Check if SHAP analysis was successful
                     if (!exists("FilterData") || is.null(FilterData) || 
                         !is.list(FilterData) || is.null(FilterData$ig)) {
                       # If SHAP analysis failed, show error and set FilterData to NULL
                       shinyalert(title = "SHAP Analysis Error", type = "error",
                                 text = "SHAP values analysis failed. Please check your data and try again.", 
                                 closeOnClickOutside = TRUE)
                       FilterData <- NULL  # Set to NULL so it's handled later
                     } else {
                       # Ensure newdata is also set (some functions might use NewData with capital N)
                       if (!is.null(FilterData$newdata)) {
                         NewData <- FilterData$newdata
                       }
                     }
                   }
                   else if (input$VariableM == "NoMethod" & input$P_method == "Empty" & 
                           input$ensembleVar == "NoMethod" & input$ensemblePvalue == "NoMethod" &
                           input$ensembleWrapper == "NoMethod" & 
                           input$ML_Method != "Empty" & input$ML_Method != "shap_rf" & input$ML_Method != "shap_xgb" &&
                           (is.null(input$SHAP_Method) || input$SHAP_Method == "Empty")) {
                     FilterData <- SelectionFilter1(
                       data = RDS_file1[,-ncol(RDS_file1)],
                       Labels = RDS_file1[, ncol(RDS_file1)],
                       MLmethod = input$ML_Method
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
                    # Ensure FilterData exists before accessing it
                    if (!exists("FilterData") || is.null(FilterData) || 
                        !is.list(FilterData) || is.null(FilterData$ig)) {
                      shinyalert(title = "Analysis Error", type = "error",
                                text = "Analysis failed. Please check your method selections and try again.", 
                                closeOnClickOutside = TRUE)
                      return()
                    }
                    
                    iG <- FilterData$ig
                    newdata <- FilterData$newdata
                    if (exists("iG")) {
                   output$TheBarPlot <- renderGirafe({
                     dfbar <- as.data.frame(head(iG, input$genes))
                     dfbar$gene <- factor(rownames(dfbar), levels = rev(rownames(dfbar)))
                     colnames(dfbar)[1] <- "score"

                     p <- ggplot2::ggplot(dfbar, ggplot2::aes(x = gene, y = score)) +
                       ggiraph::geom_bar_interactive(
                         ggplot2::aes(
                           tooltip = paste0("<b>", gene, "</b><br>Score: ", round(score, 4)),
                           data_id = gene,
                           fill = score
                         ),
                         stat = "identity", width = 0.7
                       ) +
                       ggplot2::coord_flip() +
                       ggplot2::scale_fill_gradient(low = "#cbd5e1", high = "#1e40af", guide = "none") +
                       ggplot2::theme_minimal(base_size = 13) +
                       ggplot2::theme(
                         axis.text.y = ggplot2::element_text(size = 10, color = "#1e293b"),
                         axis.text.x = ggplot2::element_text(size = 10),
                         panel.grid.major.y = ggplot2::element_blank(),
                         panel.grid.minor = ggplot2::element_blank(),
                         plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
                       ) +
                       ggplot2::labs(x = NULL, y = "Importance Score",
                            title = "Potential Biomarker Genes")

                     ggiraph::girafe(
                       ggobj = p,
                       width_svg = 10,
                       height_svg = max(5, nrow(dfbar) * 0.22),
                       options = list(
                         ggiraph::opts_hover(css = "fill:#3b82f6;stroke:#1e40af;cursor:pointer;"),
                         ggiraph::opts_tooltip(css = "background:rgba(15,23,42,0.95);color:white;padding:8px 12px;border-radius:8px;font-size:13px;"),
                         ggiraph::opts_toolbar(saveaspng = TRUE)
                       )
                     )
                   }) }else {
                     showModal(modalDialog(
                       title = "Error",
                       div(style = "color: red;", "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more methods can be select.")
                     ))
                     
                    
                     
                   }
                   }, error = function(e) {
                     showNotification(paste("Analysis error:", e$message), type = "error", duration = 15)
                   })
                   
                   tryCatch({
                     if  (exists("iG") ) {
                       
                   output$GenesList = DT::renderDataTable({
                     s1 = lapply(input$genes,
                                 function(i)
                                   iG[1:i,])
                     s1 = as.data.frame(s1)
                     rownames(s1) = make.names(row.names(iG)[1:nrow(s1)], unique = TRUE)

                     colnames(s1) = paste("# of Genes which operate as Biomarkers: ", nrow(iG))
                     s1[, 1] = rownames(s1)

                     DT::datatable(s1, rownames = FALSE,
                       options = list(pageLength = 20, scrollY = "400px", dom = 'ftip'),
                       class = 'cell-border stripe')
                   })} else {
                     showModal(modalDialog(
                       title = "Message",
                       "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more can be select."
                     ))
                   }
                   }, error = function(e) {
                     showNotification(paste("Analysis error:", e$message), type = "error", duration = 15)
                   })
                   
                   
                   tryCatch({
                     if (exists("FilterData") ) {
                       
                       iG <- FilterData$ig
                       NewData = FilterData$newdata
                     
                       
                   output$HeatmapList = DT::renderDataTable({
                     if (input$HeatMap1 == TRUE) {
                       S <- complexHeatMapFun(NewData,iG,Plot=FALSE)
                       DT::datatable(S, rownames = FALSE,
                         options = list(pageLength = 20, scrollY = "400px", dom = 'ftip'),
                         class = 'cell-border stripe')
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
                     showNotification(paste("Analysis error:", e$message), type = "error", duration = 15)
                   })
                   
                   
                   
                   output$downloadData <- downloadHandler(
                     filename = function() {
                       paste0("scGenesFinder-results-", Sys.Date(), ".zip")
                     },
                     content = function(file) {
                       dpi <- suppressWarnings(as.numeric(input$downloadDpi))
                       if (is.na(dpi) || !dpi %in% c(300, 600)) dpi <- 300

                       export_dir <- tempfile("scgenes-download-")
                       dir.create(export_dir)
                       on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

                       # One PNG device per plot, closed immediately after that plot - so an
                       # error in one plot cannot leave a device open or corrupt another's file.
                       render_png <- function(png_name, draw) {
                         path <- file.path(export_dir, png_name)
                         grDevices::png(path, width = 11, height = 7, units = "in", res = dpi)
                         tryCatch({
                           draw()
                           TRUE
                         }, error = function(e) {
                           showNotification(paste0(png_name, " was skipped: ", conditionMessage(e)), type = "warning")
                           FALSE
                         }, finally = grDevices::dev.off())
                       }

                       parts <- "FilterData.csv"
                       write.csv(NewData, file.path(export_dir, "FilterData.csv"))

                       if (isTRUE(input$HeatMap1) &&
                           render_png("ExpressionHeatmap.png", function() complexHeatMapFun(NewData, iG, Plot = TRUE))) {
                         parts <- c(parts, "ExpressionHeatmap.png")
                       }

                       if (render_png("KnnClassification.png",
                                      function() KnnClassifier(data = NewData, iG, Labels = NewData[, ncol(NewData)]))) {
                         parts <- c(parts, "KnnClassification.png")
                       }

                       zip::zip(file, files = parts, root = export_dir)
                     },
                     contentType = "application/zip"
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
    CSV_file <- read_uploaded_csv(
      input$file1,
      header = if (is.null(input$header)) TRUE else as.logical(input$header),
      sep = if (is.null(input$sep)) "," else input$sep,
      quote = if (is.null(input$quote)) "\"" else input$quote
    )
    if (input$GENEid == "EnsemblGenes") {
     # colnames(CSV_file) <- tools::toTitleCase(tolower(colnames(CSV_file)))
      
      CSV_file1 <- LexikonFun(CSV_file, input$organismus, input$GENEid)
      rownames(CSV_file1)=rownames(CSV_file)
      validate(need(
        ncol(CSV_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
      
    }
    if (input$GENEid == "ENTREZID") {
      #colnames(CSV_file) <- tools::toTitleCase(tolower(colnames(CSV_file)))
      CSV_file1 <- LexikonFun(CSV_file, input$organismus, input$GENEid)
      rownames(CSV_file1)=rownames(CSV_file)
      validate(need(
        ncol(CSV_file1) != 0,
        "The Genes Id or the Organismus is not correct"
      ))
      
    }
    if (input$GENEid == "SYMBOL") {
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
          input$ML_Method == "Empty" &&
          (is.null(input$SHAP_Method) || input$SHAP_Method == "Empty")) {
        print("EnsemleMethod")
        FilterData <- EnsemleMethod(obj = CSV_file1[,-ncol(CSV_file1)],
                                    EnseLabels = CSV_file1[, ncol(CSV_file1)])
      }
      
      
      else if (input$P_method == "Empty" &
               input$VariableM == "NoMethod"   &
               input$ensembleVar=="NoMethod" & 
               input$ensemblePvalue=="NoMethod" &
               input$ensembleWrapper=="NoMethod" & 
               ((!is.null(input$SHAP_Method) && input$SHAP_Method != "Empty" && 
                (input$SHAP_Method == "shap_rf" | input$SHAP_Method == "shap_xgb")) ||
                (input$ML_Method == "shap_rf" | input$ML_Method == "shap_xgb"))) {
        
        # Handle SHAP Values methods for CSV data - check both ML_Method and SHAP_Method
        shap_method <- NULL
        if (!is.null(input$SHAP_Method) && input$SHAP_Method != "Empty" && 
            (input$SHAP_Method == "shap_rf" || input$SHAP_Method == "shap_xgb")) {
          shap_method <- input$SHAP_Method
        } else if (input$ML_Method == "shap_rf" || input$ML_Method == "shap_xgb") {
          shap_method <- input$ML_Method
        }
        
        # Get importance limit
        importance_limit <- input$importanceLimit
        if (!is.null(input$SHAP_importanceLimit) && input$SHAP_importanceLimit != input$importanceLimit) {
          importance_limit <- input$SHAP_importanceLimit
        }
        
        # Call appropriate SHAP method
        if (!is.null(shap_method)) {
          if (shap_method == "shap_rf") {
            FilterData <- ShapValuesFilter(
              data = CSV_file1[,-ncol(CSV_file1)],
              Labels = CSV_file1[, ncol(CSV_file1)],
              MLmethod = "rf",
              importanceLimit = importance_limit
            )
          } else if (shap_method == "shap_xgb") {
            FilterData <- ShapValuesFilter(
              data = CSV_file1[,-ncol(CSV_file1)],
              Labels = CSV_file1[, ncol(CSV_file1)],
              MLmethod = "xgbTree",
              importanceLimit = importance_limit
            )
          }
        }
        
        # Check if SHAP analysis was successful
        if (!exists("FilterData") || is.null(FilterData) || 
            !is.list(FilterData) || is.null(FilterData$ig)) {
          # If SHAP analysis failed, show error and set FilterData to NULL
          shinyalert(title = "SHAP Analysis Error", type = "error",
                    text = "SHAP values analysis failed. Please check your data and try again.", 
                    closeOnClickOutside = TRUE)
          FilterData <- NULL  # Set to NULL so it's handled later
        } else {
          # Ensure newdata is also set (some functions might use NewData with capital N)
          if (!is.null(FilterData$newdata)) {
            NewData <- FilterData$newdata
          }
        }
      }
      else if (input$P_method == "Empty" &
               input$VariableM == "NoMethod"   &
               input$ensembleVar=="NoMethod" & 
               input$ensemblePvalue=="NoMethod" &
               input$ensembleWrapper=="NoMethod" & 
               input$ML_Method != "Empty" & input$ML_Method != "shap_rf" & input$ML_Method != "shap_xgb" &&
               (is.null(input$SHAP_Method) || input$SHAP_Method == "Empty")) {
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
        # Ensure FilterData exists before accessing it
        if (!exists("FilterData") || is.null(FilterData) || 
            !is.list(FilterData) || is.null(FilterData$ig)) {
          shinyalert(title = "Analysis Error", type = "error",
                    text = "Analysis failed. Please check your method selections and try again.", 
                    closeOnClickOutside = TRUE)
          return()
        }
        
        iG <- FilterData$ig
        newdata <- FilterData$newdata
        if (exists("iG")) {
          output$TheBarPlot <- renderGirafe({
            dfbar <- as.data.frame(head(iG, input$genes))
            dfbar$gene <- factor(rownames(dfbar), levels = rev(rownames(dfbar)))
            colnames(dfbar)[1] <- "score"

            p <- ggplot2::ggplot(dfbar, ggplot2::aes(x = gene, y = score)) +
              ggiraph::geom_bar_interactive(
                ggplot2::aes(
                  tooltip = paste0("<b>", gene, "</b><br>Score: ", round(score, 4)),
                  data_id = gene,
                  fill = score
                ),
                stat = "identity", width = 0.7
              ) +
              ggplot2::coord_flip() +
              ggplot2::scale_fill_gradient(low = "#cbd5e1", high = "#1e40af", guide = "none") +
              ggplot2::theme_minimal(base_size = 13) +
              ggplot2::theme(
                axis.text.y = ggplot2::element_text(size = 10, color = "#1e293b"),
                axis.text.x = ggplot2::element_text(size = 10),
                panel.grid.major.y = ggplot2::element_blank(),
                panel.grid.minor = ggplot2::element_blank(),
                plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5)
              ) +
              ggplot2::labs(x = NULL, y = "Importance Score",
                   title = "Potential Biomarker Genes")

            ggiraph::girafe(
              ggobj = p,
              width_svg = 10,
              height_svg = max(5, nrow(dfbar) * 0.22),
              options = list(
                ggiraph::opts_hover(css = "fill:#3b82f6;stroke:#1e40af;cursor:pointer;"),
                ggiraph::opts_tooltip(css = "background:rgba(15,23,42,0.95);color:white;padding:8px 12px;border-radius:8px;font-size:13px;"),
                ggiraph::opts_toolbar(saveaspng = TRUE)
              )
            )
          }) }else {
            showModal(modalDialog(
              title = "Error",
              div(style = "color: red;", "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more methods can be select.")
            ))
            
            
            
          }
      }, error = function(e) {
        showNotification(paste("Analysis error:", e$message), type = "error", duration = 15)
      })
      
      tryCatch({
        if  (exists("iG") ) {
          
          output$GenesList = DT::renderDataTable({
            s1 = lapply(input$genes,
                        function(i)
                          iG[1:i,])
            s1 = as.data.frame(s1)
            rownames(s1) = make.names(row.names(iG)[1:nrow(s1)], unique = TRUE)

            colnames(s1) = paste("# of Genes which operate as Biomarkers: ", nrow(iG))
            s1[, 1] = rownames(s1)

            DT::datatable(s1, rownames = FALSE,
              options = list(pageLength = 20, scrollY = "400px", dom = 'ftip'),
              class = 'cell-border stripe')
          })} else {
            showModal(modalDialog(
              title = "Message",
              "The analysis could not be executed; something is wrong with your selection. Make sure that the data you uploaded is in the correct format and that only one methoth from the Gene Selection field is selected; only in the Ensemble Aproach tab more can be select."
            ))
          }
      }, error = function(e) {
        showNotification(paste("Analysis error:", e$message), type = "error", duration = 15)
      })
      
      
      tryCatch({
        if (exists("FilterData") ) {
          iG <- FilterData$ig
          NewData = FilterData$newdata
          
          output$HeatmapList = DT::renderDataTable({
            if (input$HeatMap1 == TRUE) {
              S <- complexHeatMapFun(NewData,iG,Plot=FALSE)
              DT::datatable(S, rownames = FALSE,
                options = list(pageLength = 20, scrollY = "400px", dom = 'ftip'),
                class = 'cell-border stripe')
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
        showNotification(paste("Analysis error:", e$message), type = "error", duration = 15)
      })
      
      output$downloadData <- downloadHandler(
        filename = function() {
          paste0("scGenesFinder-results-", Sys.Date(), ".zip")
        },
        content = function(file) {
          dpi <- suppressWarnings(as.numeric(input$downloadDpi))
          if (is.na(dpi) || !dpi %in% c(300, 600)) dpi <- 300

          export_dir <- tempfile("scgenes-download-")
          dir.create(export_dir)
          on.exit(unlink(export_dir, recursive = TRUE, force = TRUE), add = TRUE)

          # One PNG device per plot, closed immediately after that plot - so an
          # error in one plot cannot leave a device open or corrupt another's file.
          render_png <- function(png_name, draw) {
            path <- file.path(export_dir, png_name)
            grDevices::png(path, width = 11, height = 7, units = "in", res = dpi)
            tryCatch({
              draw()
              TRUE
            }, error = function(e) {
              showNotification(paste0(png_name, " was skipped: ", conditionMessage(e)), type = "warning")
              FALSE
            }, finally = grDevices::dev.off())
          }

          parts <- "FilterData.csv"
          write.csv(NewData, file.path(export_dir, "FilterData.csv"))

          if (isTRUE(input$HeatMap1) &&
              render_png("ExpressionHeatmap.png", function() complexHeatMapFun(NewData, iG, Plot = TRUE))) {
            parts <- c(parts, "ExpressionHeatmap.png")
          }

          if (render_png("KnnClassification.png",
                         function() KnnClassifier(data = NewData, iG, Labels = NewData[, ncol(NewData)]))) {
            parts <- c(parts, "KnnClassification.png")
          }

          zip::zip(file, files = parts, root = export_dir)
        },
        contentType = "application/zip"
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
