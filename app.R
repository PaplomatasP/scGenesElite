source("ui.R", local = TRUE)
# Define server logic to read selected file ----
server <- function(input, output) { 
  #  shinyjs::runjs("$('navbarPage > *').css('zoom', '10%');")
  
  tags$head(tags$script(src = "https://code.jquery.com/jquery-3.6.0.min.js"))
  
  
  # shinyjs::runjs("$(window).off('resize');")
  #loAd all necessary scripts
  source("./Scripts/VarReductionFilter.R", local = TRUE)
  source('./Scripts/SelectionFilter.R', local = TRUE)
  source('./Scripts/StatisticalPvalueFilter.R', local = TRUE)
  source('./Scripts/ImportanceFilter.R', local = TRUE)
  source('./Scripts/Lexikon.R', local = TRUE)
  source('./Scripts/VariableGenesMethods.R', local = TRUE)
  source('./Scripts/pathway.R', local = TRUE)
  source('./Scripts/EnsemleMethod.R', local = TRUE)
  source('./Scripts/HeatmapGraph.R', local = TRUE)
  source('./Scripts/KnnClassifier.R', local = TRUE)
  source('./Scripts/MainSelectionFun.R', local = TRUE)
  #Load the KEGG pahtway dataframe
  Hpaths = readRDS("./data/KEGGpaths.rds")
  
  # Download the example data from github repository
  output$Example <- downloadHandler(
    filename = "ExampleData.csv",
    content = function(file) {
      withProgress(message = 'Data are downloading',
                   detail = '        This may take a while...',
                   value = 0,
                   {         
                     url = "https://raw.githubusercontent.com/PaplomatasP/scGenesElite/Master/scgenes/data/ExampleData.csv"
                     {
                       incProgress(4 / 10)
                       Sys.sleep(0.25)
                     }
                     
                     ExampleData = readr::read_csv(url)
                     ExampleData = data.frame(ExampleData)
                     rownames(ExampleData) = ExampleData[, 1]
                     ExampleData = ExampleData[,-c(1)]
                     write.csv(ExampleData, file, row.names = TRUE)
                     {
                       incProgress(10)
                       Sys.sleep(0.45)
                     }
                   })
    }
    
  )
  
  #Read the import .rds dataset and visualize it.
  output$Rvalue <- renderTable({
    req(input$rdsFile)
    
    RDS_file <- input$rdsFile
    if (is.null(RDS_file)) {
      return()
    } else {
      rd <- readRDS(RDS_file$datapath)
      return(rd[1:10, (ncol(rd) - 10):ncol(rd)])
    }
    
  })
  
  #Read the import .csv dataset and visualize it.
  output$contents <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # as uploads a file, head of that data file by default,
    # or all rows if selected, will be shown.
    
    req(input$file1)
    df  = read.csv(
      input$file1$datapath,
      header = input$header,
      sep = input$sep,
      quote = input$quote
    )
    CSV_file <- as.data.frame(df)
    if (input$disp == "head") {
      return(CSV_file[1:10, (ncol(CSV_file) - 10):ncol(CSV_file)])
    }
    else {
      return(df)
    }
    
    
  })
  
  
  
  
  
  
  # active the enrichments' Analysis process after click the button
  #Here we load all available ontologies and do a preprocessing in order to display them in bar chart form.
  output$Enrichment = eventReactive(input$click1,  {
    FilterData <- MethodData()
    iG <- FilterData$ig
    rownames(iG) <- tools::toTitleCase(tolower(rownames(iG)))
    if (input$all) {
      withProgress(message = 'Please wait........', value = 0, {
        {
          incProgress(3 / 10)
          Sys.sleep(0.10)
        }
        if (websiteLive) {
          dbs <-
            c(
              "KEGG_2021_Human",
              "WikiPathway_2021_Human",
              "BioPlanet_2019",
              "BioCarta_2016",
              "Reactome_2016",
              "MSigDB_Hallmark_2020",
              "GO_Biological_Process_2021",
              "GO_Molecular_Function_2021",
              "GO_Cellular_Component_2021",
              "MGI_Mammalian_Phenotype_Level_4_2021",
              "Human_Phenoty pe_Ontology",
              "Jensen_DISEASES",
              "DisGeNET",
              "DSigDB",
              "DrugMatrix",
              "OMIM_Disease", 
              "HDSigDB_Human_2021",
              "COVID-19_Related_Gene_Sets_2021"
            )
          
          enriched <- enrichr(rownames(iG)[1:input$genes1], dbs)
          count = 0
          n = list()
          for (i in seq(1:18)) {
            if (nrow(enriched[[i]]) == 0) {
              n[[i]] <- i
              count = count + 1
            }
            else if (count == 0) {
              n[[i]] = 0
            }
            
          }
          
          
          nList = unlist(n)
          nList = nList[nList > 0]
          
          
          if (length(nList) >= 1) {
            for (i in (1:length(nList))) {
              nNew = as.numeric(nList[i])
              enriched[[nNew]] <-
                as.data.frame(matrix(rbinom(10 * 10, 1, 0), ncol = 9))
              colnames(enriched[[nNew]]) = c(
                "Term",
                "Overlap",
                "P.value",
                "Adjusted.P.value",
                "Old.P.value",
                "Old.Adjusted.P.value",
                "Odds.Ratio",
                "Combined.Score",
                "Genes"
              )
            }
          }
          
          {
            incProgress(5 / 10)
            Sys.sleep(0.10)
          }
          BBP <- list()
          for (i in seq_along(enriched)) {
            print(i)
            enriched[[i]] <-
              data.frame(enriched[[i]][order(enriched[[i]]$Combined.Score, decreasing = TRUE), drop = FALSE, ])
            
            
            
            
            allPlotData = as.data.frame(enriched[[i]][1:10,])
            allPlotData <-
              allPlotData[order(allPlotData$Combined.Score, decreasing = FALSE), drop = FALSE, ]
            
            
            BBP[[i]] <-
              ggplot(data = allPlotData, aes(
                x = reorder(Term, +Combined.Score),
                y = Combined.Score
              )) +
              ggtitle(dbs[i]) +
              geom_bar(stat = "identity",
                       position = "dodge",
                       aes(fill = Combined.Score)) + coord_flip() + theme_gray() +
              theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                legend.position = "none"
              ) +
              geom_text(
                aes(label = Term),
                hjust = -0.1,
                size = 3,
                position = position_stack(vjust = 0) ,
                inherit.aes = TRUE
              ) +
              scale_fill_distiller(name = "Value",
                                   palette = "Reds",
                                   direction = 1)
            
            output$BioBarPlot <- renderPlot({
              mp <- marrangeGrob(BBP, nrow = 6, ncol = 3)
              mp
              
            })
            
            {
              incProgress(10 / 10)
              Sys.sleep(0.50)
            }
            
            plotInput <- reactive({
            
              mp <- marrangeGrob(BBP, nrow = 6, ncol = 3)
              
             
            })
            
            #Download button for the images and dataset from the enrichment results.
            # output$GP <- downloadHandler(
            #   filename = function() { paste("Enrichment_Analysis", '.png', sep='') },
            #   content = function(file) {
            #     filePath <- file.path(tempdir(), "Enrichment_Analysis.png")
            #     device <- function(..., width, height) grDevices::png(filePath, width = 18, height = 19, res = 300, units = "in")
            #     ggsave(filePath, plot = plotInput(), device = device)
            #     file.copy(filePath, file)
            #   }
            # )
            
            
        #  })
            
            
            #Pre-processing of download objects.
            to_download <-
              reactiveValues(
                KEGG_2021_Human = as.data.frame(enriched[[1]][1:10,][c(1, 2, 3, 9)]),
                WikiPathway_2021_Human = as.data.frame(enriched[[2]][1:10,][c(1, 2, 3, 9)]),
                BioPlanet_2019 = as.data.frame(enriched[[3]][1:10,][c(1, 2, 3, 9)]),
                BioCarta_2016 = as.data.frame(enriched[[4]][1:10,][c(1, 2, 3, 9)]),
                Reactome_2016 = enriched[[5]][1:10,][c(1, 2, 3, 9)],
                MSigDB_Hallmark_2020 = enriched[[6]][1:10,][c(1, 2, 3, 9)],
                GO_Biological_Process_2021 = enriched[[7]][1:10,][c(1, 2, 3, 9)],
                GO_Molecular_Function_2021 = enriched[[8]][1:10,][c(1, 2, 3, 9)],
                GO_Cellular_Component_2021 = enriched[[9]][1:10,][c(1, 2, 3, 9)],
                MGI_Mammalian_Phenotype_Level_4_2021 = enriched[[10]][1:10,][c(1, 2, 3, 9)],
                Human_Phenotype_Ontology = enriched[[11]][1:10,][c(1, 2, 3, 9)],
                Jensen_DISEASES = enriched[[12]][1:10,][c(1, 2, 3, 9)],
                DisGeNET = enriched[[13]][1:10,][c(1, 2, 3, 9)],
                DSigDB = enriched[[14]][1:10,][c(1, 2, 3, 9)],
                DrugMatrix = enriched[[15]][1:10,][c(1, 2, 3, 9)],
                OMIM_Disease = enriched[[16]][1:10,][c(1, 2, 3, 9)],
                HDSigDB_Human_2021 = enriched[[17]][1:10,][c(1, 2, 3, 9)],
                COVID19_Related_Gene_Sets_2021 = enriched[[18]][1:10,][c(1, 2, 3, 9)]
              )
            
        #     output$GT <- downloadHandler(
        #       filename = function() {
        #         paste("Enrichment_Analysis_",
        #               Sys.Date(),
        #               ".zip",
        #               sep = "")
        #       },
        #       content = function(file) {
        #         temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
        #         dir.create(temp_directory)
        #         
        #         reactiveValuesToList(to_download) %>%
        #           imap(function(x, y) {
        #             if (!is.null(x)) {
        #               file_name <- glue("{y}_data.csv")
        #               readr::write_csv(x, file.path(temp_directory, file_name))
        #             }
        #           })
        #         
        #         
        #         zip::zip(
        #           zipfile = file,
        #           files = dir(temp_directory),
        #           root = temp_directory
        #         )
        #         
        #         
        #         
        #       },
        #       contentType = "application/zip"
        #       
        #     )
        #     
        #     
          }
        }
        
        
      })
      
      
      
      
    }
    
    #Preprocessing to load only the requested ontologies and  display them in bar chart format
    else{
      
      withProgress(message = 'Please wait........', value = 0
                   , {
                     {
                       incProgress(3 / 10)
                       Sys.sleep(0.10)
                     }
                     dbs1 <- c(input$BP, input$BO, input$DD)
                     oo = which(dbs1 != "-")
                     dbs = dbs1[c(oo)]
                     
                     if (websiteLive) {
                     
                       enriched <- enrichr(rownames(iG)[1:input$genes1], dbs)
                       
                       
                       count = 0
                       n = list()
                       for (i in seq_along(enriched)) {
                         if (nrow(enriched[[i]]) == 0) {
                           n[[i]] <- i
                           count = count + 1
                         }
                         else if (count == 0) {
                           n = 0
                         }
                         
                       }
                       
                       
                       nList = unlist(n)
                       nList = nList[nList > 0]
                       
                       
                       if (length(nList) >= 1) {
                         for (i in (1:length(nList))) {
                           nNew = as.numeric(nList[i])
                           enriched[[nNew]] <-
                             as.data.frame(matrix(rbinom(10 * 10, 1, 0), ncol = 9))
                           colnames(enriched[[nNew]]) = c(
                             "Term",
                             "Overlap",
                             "P.value",
                             "Adjusted.P.value",
                             "Old.P.value",
                             "Old.Adjusted.P.value",
                             "Odds.Ratio",
                             "Combined.Score",
                             "Genes"
                           )
                         }
                       }
                       
                       {
                         incProgress(3 / 10)
                         Sys.sleep(0.10)
                       }
                       
                       BBP <- list()
                       for (i in seq_along(enriched)) {
                         enriched[[i]] <-
                           data.frame(enriched[[i]][order(enriched[[i]]$Combined.Score, decreasing = TRUE), drop = FALSE, ])
                         
                         allPlotData = as.data.frame(enriched[[i]][1:10,])
                         allPlotData <-
                           allPlotData[order(allPlotData$Combined.Score, decreasing = FALSE), drop = FALSE, ]
                         
                         {
                           incProgress(10 / 10)
                           Sys.sleep(0.10)
                         }
                         BBP[[i]] <-
                           ggplot(data = allPlotData, aes(
                             x = reorder(Term, +Combined.Score),
                             y = Combined.Score
                           )) +
                           ggtitle(dbs[i]) +
                           geom_bar(stat = "identity",
                                    position = "dodge",
                                    aes(fill = Combined.Score)) + coord_flip() + theme_gray() +
                           theme(
                             axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks.x = element_blank(),
                             axis.text.y = element_blank(),
                             axis.title.y = element_blank(),
                             legend.position = "none"
                           ) +
                           geom_text(
                             aes(label = Term),
                             hjust = -0.1,
                             size = 3,
                             position = position_stack(vjust = 0) ,
                             inherit.aes = TRUE
                           ) +
                           scale_fill_distiller(name = "Value",
                                                palette = "Reds",
                                                direction = 1)
                         
                         output$BioBarPlot <- renderPlot({
                           mp <- marrangeGrob(BBP, nrow = 2, ncol = 2)
                           mp
                           
                         })
                         plotInput <- reactive({

                           mp <- marrangeGrob(BBP, nrow = 3, ncol = 1)
                        
                          
                         })
                         
                         #Download button for the images and dataset from the enrichment results.
                         
                         # output$GP <- downloadHandler(
                         #   filename = function() { paste("Enrichment_Analysis", '.png', sep='') },
                         #   content = function(file) {
                         #     filePath <- file.path(tempdir(), "Enrichment_Analysis.png")
                         #     device <- function(..., width, height) grDevices::png(filePath, width = 18, height = 19, res = 300, units = "in")
                         #     ggsave(filePath, plot = plotInput(), device = device)
                         #     file.copy(filePath, file)
                         #   }
                         # )
                         
                         
                         to_download <- reactiveValues(
                           Enrichment1 = if (length(enriched) >= 1) {
                             Pathway_Enrichment = as.data.frame(enriched[[1]][1:10,][c(1, 2, 3, 9)])
                           },
                           
                           Enrichment2 = if (length(enriched) >= 2) {
                             Ontologies_Enrichment = as.data.frame(enriched[[2]][1:10,][c(1, 2, 3, 9)])
                           },
                           
                           Enrichment3 =  if (length(enriched) >= 3) {
                             Dieseses_Drugs_Enrichment = as.data.frame(enriched[[3]][1:10,][c(1, 2, 3, 9)])
                           },
                           
                           
                           
                         )
                         
                         
                         # output$GT <- downloadHandler(
                         #   filename = function() {
                         #     paste("Enrichment_Analysis_",
                         #           Sys.Date(),
                         #           ".zip",
                         #           sep = "")
                         #   },
                         #   content = function(file) {
                         #     temp_directory <- file.path(tempdir(), as.integer(Sys.time()))
                         #     dir.create(temp_directory)
                         #     
                         #     reactiveValuesToList(to_download) %>%
                         #       imap(function(x, y) {
                         #         if (!is.null(x)) {
                         #           file_name <- glue("{y}_data.csv")
                         #           readr::write_csv(x, file.path(temp_directory, file_name))
                         #         }
                         #       })
                         #     
                         #     
                         #     zip::zip(
                         #       zipfile = file,
                         #       files = dir(temp_directory),
                         #       root = temp_directory
                         #     )
                         #     
                         #     
                         #     
                         #   },
                         #   contentType = "application/zip"
                         #   
                         # )
                         
                         
                       }
                       
                       
                       
                       print("finish")
                     }
                   })
    }
    
    
    
    
  })
  
  #load the pathway dataframe and  display it.
  PathsData <- data.frame(iD = rownames(Hpaths),
                          Paths = Hpaths[, 1])
  
  
  output$KEGG <- renderDT({
    
    datatable(PathsData, options = list(
      scrollY = "200px",
      scrollCollapse = TRUE,
      pagination = FALSE,
      pageLength = nrow(PathsData)
    ))
  })
  
  
  # connect the go button to the active function by clicking .
  text1 = eventReactive(input$click2, {
    return(input$click2)
  })
  #doing the pre-process check for the Organismus and plot the pathaway.
  output$KEGGmap <- renderPlot({
    text1()
    if (input$organismus == "Human") {
      Species = "hsa"
      
    }
    else{
      Species = "mmu"
    }
    withProgress(message = 'Please wait........', value = 0, {
      {
        incProgress(8 / 10)
        Sys.sleep(0.10)
      }
      FilterData <- MethodData()
      iG <- FilterData$ig
      # A dictionary for translating gene ID to a gene symbol before executing the KEGG map process.
      ScaledData <- iGlexikon(iG, input$GENEid,input$organismus)
    
      plot_pathview(gene.data =ScaledData,
                    pathway.id = input$inText,
                    species = Species, out.suffix = "", save_image = FALSE,
                    keys.align = "y", kegg.native = T, match.data = T, multi.state = T,
                    same.layer = F)
     
      
      

      
      {
        incProgress(10 / 10)
        Sys.sleep(0.10)
      }
      print("Fire")
      
    })
    
  })
  withProgress(message = 'Please wait........', value = 0, {
    {
      incProgress(2 / 10)
      Sys.sleep(0.10)
    }
    # connect the "Plotting Graphs" button and active  by clicking .
    text3 = eventReactive(input$run_button1, {
      return(input$run_button1)
      
    })
    
    
    
    
    {
      incProgress(7 / 10)
      # connect the " lysis" button and active  by clicking .
      output$text = eventReactive(input$click, {
        FilterData <- MethodData()
        iG <- FilterData$ig
        NewData = FilterData$newdata 
        
        #Plot the Similarity Graph
        output$graph <- renderVisNetwork({
          text3()
          
          if (input$graph1) {
            FilterData <- MethodData()
            NewData = FilterData$newdata 
            iG <- FilterData$ig
            
            isolate({
           
              GraphsFun(NewData,iG)
            })
            
            
          }
        })
        
      })
      Sys.sleep(0.10)
    }
    
    {
      incProgress(9 / 10)
      #Plot the PPI network Analysis Graph
      output$PPInetwork <- renderPlot(execOnResize = FALSE, {
        text3()
        if (input$PPInetwork1) {
          FilterData <- MethodData()
          iG <- FilterData$ig
         
          print("i am here to PPInetwork1")
          isolate({
           
            PPInetwork(iG)
          })
          
          
        }
      })
    }
  })
}

#  Run the app ----
shinyApp(ui, server)
