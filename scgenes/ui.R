#Load the necessary libraries
libs <-
  c(
    "shiny","shinyalert",
    "shinythemes",
    "enrichR",
    "ggplot2",
    "gridExtra",
    "glue",
    "tidyverse",
    "shinyWidgets",
    "shinydashboard",
    "twoddpcr",
    "SCMarker",
    "DT",
    "Seurat",
    "pathview",
    "grid",
    "png",
    "ggiraph",
    "AnnotationDbi",
    "AnnotationFilter",
    "Biobase",
    "BiocFileCache",
    "BiocGenerics",
    "BiocParallel",
    "BiocStyle",
    "BiocManager",
    "fastAdaboost",
    "votesys",
    "shinycustomloader",
    "M3Drop",
    "ComplexHeatmap",
    "igraph",
    "visNetwork",
    "SingleR",
    "shinyjs",
    "STRINGdb"
  )
lapply(libs, require, character.only = TRUE)



options(repos = BiocManager::repositories())


websiteLive <- TRUE
#Increase the size of the acceptable import dataset
options(shiny.maxRequestSize = 256 * 2048 ^ 2)





# Define UI for data upload app ----

ui <- #fluidPage(div(class = "tab-content",
  navbarPage(
  
    theme = shinythemes::shinytheme("cosmo"),
    "",
    tags$head(
      tags$style(
        "
        h1 {
          text-align: center;
          font-size: 28px;
          font-weight: bold;
        }
        p {
          font-family: system-ui;
          font-size: 15px;
        }
        
         h3 {
          font-family: system-ui;
          text-align: center;
          font-size: 22px;
          font-weight: bold;
         }
        
         .footer {
          background-color: #333;
          color: #fff;
          text-align: center;
          padding: 18px;
         }
        "
        # .main-panel {
        #  padding: 50px;
        #  background-color: #f2f2f2;
        #      
        # }
        # .header {
        #   background-color: #333;
        #   color: #fff;
        #   text-align: center;
        #   padding: 19px;
        # }
        # 
        # }
        # .sidebar-layout {
        #   padding: 20px;
        # }
        # .overview-image {
        #   display: block;
        #   margin: 0 auto;
        # }
        
   
      )
    ),
    # App title ----
    tabPanel (
     
      
        title = tags$img(src='scGenesElite.jpg.png', width = '200px', height = '80px'),
    
      
      sidebarLayout(
        position = "right",
        #sidebarPanel(
     
        div(tags$img(
          src = "overview.jpg",
          width = '520px',
           height = '600px'
        ), style = "position:relative; top:0px;"),
        mainPanel(
          br(),
          br(),
          strong(
            h1(
              "A tool-guide for the identification of more robust genes-level drug targets and biomarkers for complex diseases."
            )
          ),
      
          tags$hr(),
          br(),
          br(),
          p(
            "The app offers a flexible platform for the identification of dominant genes in a single-cell RNA-sequencing (scRNA-seq) dataset which operate as disease biomarkers.
          It includes three different types of gene selection methods, exploiting the statistical aspect, the machine-leaning aspect along with state-of-the-art feature selection methods tailored for scRNA-seq data.
          The feature selection operation modes include 20 different methodologies covering a broad range of such approaches. The extracted gene list is further examined for enrichment in various biological and pharmacological features
          including (i) pathway terms, (ii) gene ontology (GO) terms of molecular function, biological processes and cellular components, (iii) disease terms, (iv) drug substances based on the EnrichR tool.
          Snapshots of KEGG pathway maps enhance the investigation of the exported genes as biomarkers and provide insight into the functional and structural characteristics of the biological system under study.
          The final tab offers the user the opportunity to undertake a protein-protein interaction (PPI) network analysis and similarity graph analysis. The PPI network analysis is aimed at determining the functional 
            relationships between proteins by exploring the interactions between the proteins in a biological system. This information can facilitate a more comprehensive appreciation of the biological system's overall functioning and provide valuable insights into the mechanisms of diseases and potential drug targets.
            The similarity graph analysis,
            on the other hand, enables the identification of molecular modules within genetic networks through the assessment of the similarity between gene interaction profiles within a cell."
           
          ),
          
          
          br(),
          br(),
          br(),
      
          div(
            class = "footer",
            p(
              a(href = "mailto:p.paplomatas@hotmail.com", "Feel free to contact us for any issue or question at p.paplomatas@hotmail.com", 
                style = "color:white; text-align:center; font-size: 22px; "),
            )
          ),
          br(),
          
        )
      )
    ),
    tabPanel (
      tags$head(tags$style(
        HTML("
      .shiny-output-error-validation {
        color: red;
      }
    ")
      )),
      title = tags$h3("Data Upload"),
      
      #Info button
      dropMenu(
        dropdownButton(
          "Info",
          status = 'info',
          size = "xs",
          icon = icon('info')
        ),
        h3(strong('Information')),
        br(),
       
        p(
          'The Shiny app requires single-cell RNA-sequencing data in the form of read counts with annotations.
      The dataset needs to be formatted in a matrix with size NxE, with N cell samples and E-1 gene expressions (i.e. X(i,j) is the expression value of gene j for the cell i).
      The app requires the annotation last column while it should have two classes (e.g., health-disease). The app can also be run with multiple classes.
      In this case, the exported potential biomarkers indicate the separability among all classes offering the interpretation that the reader considers. It operates more effectively if data are normalized, and the initial condition (cell classes) includes a binary case study (control - state) are computed.
         The app can also be run with multiple classes.'
        ),
        
        
        
        
        
        placement = "bottom-start",
        arrow = TRUE,
        theme = "material",
        maxWidth = 1000
      ),
      
      # Sidebar layout with input and output definitions ----
      sidebarLayout(
        sidebarPanel(
          downloadButton("Example", "Example DataSet" , style = "width:100%"),
          
          radioButtons(
            "organismus",
            "Organismus:",
            
            choices = c("H.sapiens" = "Human",
                        "M. musculus" = "Mouse"),
            selected = "Human",
            inline = TRUE
          ),
          
          radioButtons(
            "GENEid",
            "Genes ID",
            choices = c(
              "Genes Symbol" = "SYMBOL",
              "Ensembl ID" = "EnsemblGenes",
              "Entrez ID " = "ENTREZID"
            ),
            selected = "SYMBOL",
            inline = TRUE
          ),
          
          # Sidebar panel for inputs ----
          fileInput("rdsFile", "Choose RDS File", multiple = TRUE),
          
          # Input: Select a file ----
          tags$hr(),
          fileInput(
            "file1",
            "Choose CSV File",
            multiple = TRUE,
            accept = c("text/csv",
                       "text/comma-separated-values,text/plain",
                       ".csv")
            
          ),
          
          # Horizontal line ----
          tags$hr(),
          
          # Input: Checkbox if file has header ----
          checkboxInput("header", "Header", TRUE),
          
          # Input: Select separator ----
          radioButtons(
            "sep",
            "Separator",
            choices = c(
              Comma = ",",
              Semicolon = ";",
              Tab = "\t"
            ),
            selected = ","
          ),
          
          # Input: Select quotes ----
          radioButtons(
            "quote",
            "Quote",
            choices = c(
              None = "",
              "Double Quote" = '"',
              "Single Quote" = "'"
            ),
            selected = '"'
          ),
          
          # Horizontal line ----
          tags$hr(),
          
          # Input: Select number of rows to display ----
          radioButtons(
            "disp",
            "Display",
            choices = c(Head = "head",
                        All = "all"),
            selected = "head"
          )
          
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(# Output: Data file ----
                  tableOutput("contents"),
                  tableOutput("Rvalue"))
        
      )
    ),
    #Second tab
    tabPanel (
      title = tags$h3("Run Analysis"),
      # dashboardSidebar layout with input and output definitions ----
      dashboardSidebar(width = 5000,
                       column(
                         1,
                         
                         #Info Button
                         dropMenu(
                           dropdownButton(
                             "Info",
                             status = 'info',
                             size = "xs",
                             icon = icon('info')
                           ),
                           h3(strong('Variance Filter Information')),
                           br(),
                           h5(
                             div(style = "color:blue",
                                 strong("1.Remove Low Variance:"), ),
                             "
Using the nearZeroVar function from R package:Caret identifies predictors with one unique value (zero variance predictors) or predictors with both of the following characteristics: they have a small number of unique values compared to the number of samples and a high frequency of the most frequent value. A threshold option is provided as a cutoff for the percentage of distinct values in relation to the total number of samples in the dataset.
" ,
                             br(),
                             div(style = "color:blue",
                                 strong("2.Keep High Variable Genes :"), ),
                             "
The FindVariableFeatures function of the seraut package is utilized to reduce the number of genes based on their <<high variability>> within the dataset. The user can specify a threshold, which corresponds to the exact number of genes to be retained for further analysis. This function is a useful tool for identifying the most informative genes that are likely to have the greatest impact on the outcome of the analysis.
"
                           ),
                           
                           placement = "bottom-start",
                           arrow = TRUE,
                           theme = "material",
                           maxWidth = 1000
                         )
                       )),
      #------------------------------------------------------------splitLayout
      
      splitLayout(
        cellWidths = 950,
        style = "width:auto; max-height: 50vh; overflow-x: hidden; overflow-y: auto;" ,
        cellArgs = list(style = "padding: 6px"),
        
        
        inputPanel(
          radioButtons(
            "VarFilter",
            "Data Filtering:",
            
            choices = c(
              "Remove Low Variance" = "Strict_Filter",
              "Keep High Variable Genes" = "Lightweight_Filter",
              "No Filter" = "Unselect"
            ),
            selected = "Strict_Filter",
            inline = TRUE,
          ),
          #change the size
          tags$head(tags$style(
            HTML(
              '
             #uniqueCut{height: 20px}
             #uniqueCut{width:  100px}
             #nfeatures{height: 20px}
             #nfeatures{width:  100px}
             '
            )
          )),
          # numericInput function
          numericInput("uniqueCut", "Low Variance cutoff", 15, 0, 100),
          numericInput("nfeatures", "High Variable nfeatures", 2000, 0, 10000)
        ),
        
      ),
      # sidebarLayout layout with input and output definitions ----
      sidebarLayout(
        column(
          1,
          offset = 4,
          actionButton("click", "Run Analysis",
                      
                       style = "flex-start"),# #position: absolute; right: -130px;  bottom: -85px
        ),
        column(
          2 ,
          offset = 2,
          downloadLink(
            "downloadData",
            "Download",
            
            style = "color: #fefbd8; background-color: rgb(230, 50  , 50);
                              position: absolute;   bottom: -60px",
          )
        ),
        
        
        
      ),
      #------------------------------------------------------------dashboardSidebar
      dashboardSidebar(
        width = 450,
        column(
          width = 1,
          dropMenu(
            dropdownButton(
              "Info",
              status = 'info',
              size = "xs",
              icon = icon('info')
            ),
            
            h3(strong('Information')),
            br(),
            h5(
              "In this suite of analysis tools, a wide range of strategies are provided for the analysis of the data. These include methods for identifying variable genes based on single cell RNA sequencing, statistical approaches based on p-value threshold, and wrapper and tree-based machine learning models. A unique feature of this suite is the ability for the user to create their own ensemble method by combining one or more of these methods. Additionally, various threshold options are available for each method (for more information, consult the tutorial provided). It should be noted that each method employs a distinct approach and may require varying amounts of computational resources, depending on the size of the data. Therefore, it is recommended that the user employs preprocessing filters to reduce the data size and improve the reliability of the results. The normalization button allows the user to normalize their data prior to analysis. On the Ensemble tab, the options selected from the previous tabs for the various methods are applied. The genes identified through the analysis are depicted in a visual format utilizing both a barplot and a heatmap. The heatmap specifically employs cell type or State labeling as a means of organizing and presenting the data. Furthermore, a classification K-nearest neighbors (K-nn) model is executed to assess the ability of the model to accurately classify based solely on the isolated genes. The results of the K-NN analysis are presented in a comprehensive manner through the utilization of a confusion matrix and a data table."
            )
            ,
            
            
            placement = "bottom-start",
            #"light-border"
            arrow = TRUE,
            theme = "material",
            maxWidth = 1000
          )
        ),
        
        
        strong("Genes Selection Methods", style =
                 "color:black")
      ),
      
      br(),
      #TabsetS
      tabsetPanel(
        type = "tabs",
        
    
        # Variable Genes Tab
        tabPanel(
          "Variable Genes",
          
          # ---------------sidebarPanel
          sidebarPanel(
            radioButtons (
              "Norm",
              "Normalization:",
              c("Normalization" = "Normal",
                "No-Normalization" = "No_Normal"),
              selected = "Normal",
              inline = TRUE
              
            ),
            tags$hr(),
            fluidRow(
              column(
                12,
                radioButtons(
                  'VariableM',
                  label    = "Methods:",
                  choices = list(
                    "SCMarker"       = "SCMarker",
                    "DUBStepR"       = "DUBStepR",
                    "ScPNMF"         =  "ScPNMF",
                    "SelfE"          = "SelfE" ,
                    "M3Drop"          = "M3Drop" ,
                    "None" = "NoMethod"
                  ),
                  selected = "NoMethod"
                  
                ),
              ),
              strong(div(style = "color:black",
                         "Parameters Selection:", ), ),
              
              #Give the right size
              tags$head(tags$style(
                HTML(
                  '
        #cellK{height: 20px}
        #cellK{width:  80px}
        #geneK{height: 20px}
        #geneK{width:  80px}
        #n{height: 20px}
        #n{width:  80px}
        #k{height: 20px}
        #k{width:  80px}
        #np{height: 20px}
        #np{width:  80px}
        #distMethod{height: 20px}
        #distMethod{width:  80px}
        #gM{height: 20px}
        #gM{width:  80px}
        #M3dropthreshold{height: 20px}
        #M3dropthreshold{width:  90px}
        #M3Method{height: 20px}
        #M3Method{width:  80px}



    '
                )
              )),
              
              column(4,
                     div(style = "color:black",
                         "SCMarker:", ), ),
              column(4,
                     numericInput(
                       "geneK",
                       "geneK",
                       value = 20,
                       min = 0,
                       step = 1
                     )),
              column(4,
                     numericInput(
                       "cellK",
                       "cellK",
                       value = 20,
                       min = 0,
                       step = 1
                     )),
              
              column(4,
                     div(style = "color:black",
                         "DUBStepR:", ), ),
              column(4,
                     numericInput(
                       "k",
                       "k",
                       value = 10,
                       min = 1,
                       step = 1
                     )),
              
              column(4,
                     numericInput(
                       "np",
                       "num.pcs",
                       value = 10,
                       min = 1,
                       step = 1
                     )),
              column(4,
                     div(style = "color:black",
                         "ScPNMF:", ), ),
              column(4,
                     numericInput(
                       "gM",
                       "M Genes #",
                       value = 300,
                       min = 1,
                       step = 1
                     )),
              column(4,
                     selectInput(
                       "distMethod",
                       "distMethod",
                       c(
                         "EucDist" = "EucDist",
                         "KL" = "KL",
                         "DPNMF" = "DPNMF"
                       )
                     )),
              column(width = 4,
                     div(style = "color:black",
                         "M3Drop:", ), ),
              column(
                4,
                numericInput(
                  "M3dropthreshold",
                  "mt_threshold",
                  value = 0.05,
                  min = 0.001,
                  max = 1,
                  step = 0.01
                )
              ),
              column(4,
                     selectInput(
                       "M3Method", "mt_method",
                       c("fdr" = "fdr",
                         "bon" = "bon")
                     )),
              column(4,
                     div(style = "color:black",
                         "SelfE:", ), ),
              column(4,
                     numericInput(
                       "n",
                       "# Features",
                       value = 150,
                       min = 1,
                       step = 1
                     )),
              
            )  ,
            
            
            
          ),
        ),
        
        
        
        # ---------------P-value tab
        tabPanel(
          "P-value",
          
          # ---------------sidebarPanel
          sidebarPanel(
            position = "right",
            width = 4,
            
            
            
            radioButtons (
              "P_method",
              "P-Value Methods:",
              c(
                "Wilcoxon rank sum test" = "Seurat_method",
                "Beta-Poisson generalized linear model" = "BPSC_metchod",
                "Wald test" = "MAST_method",
                "Likelihood Ratio Test" = "DESeq2_method",
                "None" = "Empty"
              ),
              selected = "Empty",
            ),
            
            
            numericInput(
              "PvalueNum",
              "P-Value Threshold",
              value = 0.05,
              min = 0.01,
              max = 0.99,
              step = 0.01
            ),
            
            
          ),
          
        ),
        # ---------------Wrapper Based ML tab
        tabPanel(
          "Wrapper Based ML",
          
          
          # ---------------sidebarPanel
          sidebarPanel(
            position = "right",
            width = 4,
            
            
            
            
            
            selectInput(
              "Wrapper_ML_Method",
              "Wrapper Machine Learning Methods:",
              c(
                "Linear Discriminant Analysis" = "lda",
                "Lasso and Elastic-Net Regularized Generalized Linear Models" = "glmnet",
                "K-nearest neighbors algorithm" = "knn",
                "Support Vector Machine-Radial" = "svmRadial",
                "None" = "Empty"
                
              ),
              selected = "Empty",
            ),
            numericInput(
              "importanceLimit",
              "Significant Threshold",
              value = 30,
              min = 0,
              max = 100,
              step = 1
            ),
            
            
            
            
          ) ,
        ),
        # ---------------Tree Based ML Tab
        tabPanel(
          "Tree Based ML",
          
          
          # ---------------sidebarPanel
          sidebarPanel(
            position = "right",
            width = 4,
            
            
            
            selectInput(
              "ML_Method",
              "Tree Based Machine Learning Methods:",
              c(
                "Random Forest Algorithm" = "rf",
                "eXtreme Gradient Boosting" = "xgbTree",
                "Bagged CART" = "treebag",
                "Recursive Partitioning and Regression Trees" = "rpart",
                "C5.0 Decision Trees and Rule-Based Models" = "C5.0",
                "None" = "Empty"
                
              ),
              selected = "Empty",
            ),
            numericInput(
              "importanceLimit",
              "Significant Threshold",
              value = 30,
              min = 0,
              max = 100,
              step = 1
            ),
            
            
            
            
          ) ,
        ),
        # ---------------Ensemble Approach  Tab
        tabPanel(
          "Ensemble Approach",
          
          
          # ---------------sidebarPanel
          sidebarPanel(
            position = "right",
            width = 4,
            
            radioButtons(
              'ensembleVar',
              label    = "Variable Genes Methods:",
              choices = list(
                "SCMarker"       = "SCMarker",
                "DUBStepR"       = "DUBStepR",
                "ScPNMF"         =  "ScPNMF",
                "M3Drop"          = "M3Drop" ,
                "SelfE"          = "SelfE" ,
                
                "None" = "NoMethod"
              ),
              selected = "NoMethod"
              
            ),
            
            radioButtons(
              'ensemblePvalue',
              label    = "P-value Methods:",
              choices = list(
                "Wilcoxon rank sum test" = "Seurat_method",
                "Beta-Poisson generalized linear model" = "BPSC_metchod",
                "Wald test" = "MAST_method",
                "Likelihood Ratio Test" = "DESeq2_method",
                "None" = "NoMethod"
              ),
              selected = "NoMethod"
              
            ),
            
            radioButtons(
              'ensembleWrapper',
              label    = "Wrapper Based Methods:",
              choices = list(
                "Random Forest Algorithm" = "rf",
                "eXtreme Gradient Boosting" = "xgbTree",
                "Bagged CART" = "treebag",
                "Recursive Partitioning and Regression Trees" = "rpart",
                "C5.0 Decision Trees and Rule-Based Models" = "C5.0",
                "None" = "NoMethod"
              ),
              selected = "NoMethod"
              
            ),
            radioButtons(
              'ensembleMLBased',
              label    = "ML Based Methods:",
              choices = list(
                "Linear Discriminant Analysis" = "lda",
                "Lasso and Elastic-Net Regularized Generalized Linear Models" = "glmnet",
                "K-nearest neighbors algorithm" = "knn",
                "Support Vector Machine-Radial" = "svmRadial",
                "None" = "NoMethod"
              ),
              selected = "NoMethod"
              
            ),
            
            
            
            
            
          ) ,
        ),
        
        
      ),
      
      
      # -----------   Main panel for displaying outputs
      mainPanel(tableOutput("text")),
      
      
      
      
      
      
      # -----------   sidebarPanel
      sidebarPanel(
        position = "right",
        width = 8,
        
        sidebarLayout(
          column(
            12,
            
            sliderInput(
              "genes",
              "Number of genes:",
              
              min = 2,
              max = 500,
              value = 100
            ),
          ),
         
          column(
          
            12,
            p("⚠ Don't forget to choose the right Organismus !!! "),
            checkboxInput("HeatMap1", "Heatmap🔥", value = FALSE),
            
            column(
              width = 2,
             
              radioButtons(
                "clustering",
                "Heatmap Clustering:",
                
                choices = c("Between Groups" = "cluster_between_groups",
                            "None" = "None_Clustering"),
                selected = "None_Clustering",
                inline = TRUE
              )
            ),
            column(
              width = 2,
              radioButtons(
                "Split",
                "Heatmap Splitting:",
                
                choices = c("Cell Type" = "CellType",
                            "State" = "labels"),
                selected = "CellType",
                inline = TRUE
              )
            ),
            
          ),
          
        ),
        column(12,
               mainPanel(
                 tableOutput("error"),
                 
                 splitLayout(
                   cellWidths = 880,
                   tags$div(
                     plotOutput("TheBarPlot"),
                     width = "100%",
                     height = "400",
                     style = "margin-top: 40px;"
                   ),
                   # margin-left: 100px;
                   tags$div(tableOutput("GenesList"), style = "margin-top: 18px; height: 40vh; width: 48vh;  overflow-x: hide;
                               overflow-y: scroll; "),
                   
                  
                 ),
                 
               ), ),
        column(12,
               mainPanel(

                 splitLayout(
                   cellWidths = 950,
                   tags$div(
                     plotOutput("KnnClassifier"),
                     width = "100%",
                     height = "400",
                     style = "margin-top: 20px;"
                   ),
                   
                   
                   
                 ),
                 
               ), ),
        
        column(12,
               mainPanel(
                 
                 splitLayout(
                   cellWidths = 1000,
                   tags$div(
                     plotOutput("HeatMap"),
                     width = "100%",
                     height = "400",
                     style = "margin-top: 40px;"
                   ),
                   # margin-left: 100px;
                   tags$div(tableOutput("HeatmapList"), style = "margin-top: 18px; height: 45vh; width: 30vh;  overflow-x: hide;
                               overflow-y: scroll; "),
                   
                   
                 ),
                 
               ), ),
        
        
      ),
      
      
      
    ),
    
    # -----------   Enrichment Analysis Tab
    tabPanel (
      title = tags$h3("Enrichment Analysis"),
      
      #Info Button
      dropMenu(
        dropdownButton(
          "Info",
          status = 'info',
          size = "xs",
          icon = icon('info'),
          width = 500
        ),
        h3(strong('Information')),
        br(),
        h5(
          'In this section, the user is provided with the capability to investigate the potential biomarkers that have been identified through previous analysis. The user has the option to select either the entire set of isolated genes or a specific subset of genes. For example, if the user selects 50 genes, the top 50 genes with the highest scores will be selected for further analysis. Additionally, the user has the option to select from 18 different Pathways Datasets, which are divided into three distinct categories: Biological Pathway, Biological Ontologies, and Diseases-Drugs. The user can choose to analyze a specific pathway, a combination of three ontology terms, or all available ontology terms. The Enrichr database is utilized to perform this analysis.  '
        ),
        
        placement = "bottom-start",
        arrow = TRUE,
        theme = "material",
        maxWidth = 1000
      ),
      
      fluidRow(
        # -----------   sidebarPanel
        sidebarPanel(
          position = "right",
          width = 12,
          
          sliderInput(
            "genes1",
            "Select the number of genes for the enrichment analysis:",
            
            min = 10,
            max = 9000,
            value = 50,
            step = 1
          ),
          column(
            3,
            (""),
            checkboxInput("all", "All Available Ontology Terms", value = FALSE)
          ),
          
          
          
          column(3,
                 selectInput(
                   "BP",
                   "Biological Pathway",
                   c(
                     "-",
                     "KEGG 2021 Human" = "KEGG_2021_Human",
                     "WikiPathway 2021 Human" = "WikiPathway_2021_Human",
                     "BioPlanet 2019" = "BioPlanet_2019",
                     "BioCarta 2016" = "BioCarta_2016",
                     "MSigDB Hallmark 2020" = "MSigDB_Hallmark_2020",
                     "Reactome 2016" = "Reactome_2016"
                     
                   )
                 )),
          column(3,
                 selectInput(
                   "BO",
                   "Biological Ontologies",
                   c(
                     "-",
                     "GO Biological Process 2021" = "GO_Biological_Process_2021",
                     "GO Molecular Function 2021" = "GO_Molecular_Function_2021",
                     "GO Cellular Component 2021" = "GO_Cellular_Component_2021",
                     "MGI Mammalian Phenotype Level 4 2021" = "MGI_Mammalian_Phenotype_Level_4_2021",
                     "Human Phenotype Ontology" = "Human_Phenotype_Ontology",
                     "Jensen_DISEASES" = "Jensen_DISEASES"
                     
                     
                   )
                 )),
          column(3,
                 selectInput(
                   "DD",
                   "Diseases-Drugs",
                   c(
                     "-",
                     "DisGeNET" = "DisGeNET",
                     "DSigDB" = "DSigDB",
                     "DrugMatrix" = "DrugMatrix",
                     "OMIM Disease" = "OMIM_Disease",
                     "HDSigDB Human 2021" = "HDSigDB_Human_2021",
                     "COVID-19 Related Gene_Sets 2021" = "COVID-19_Related_Gene_Sets_2021"
                     
                     
                   )
                 )),
          column(3,
                 actionButton("click1", "Enrichment Analysis ")),
           column(
             7,
          # 
          #   downloadButton("GP", "Export Plots"),
          #   downloadButton("GT", "Export Tables"),
          ),
          br(),
          fluidRow(sidebarPanel(mainPanel(
            tags$div( plotOutput("BioBarPlot", width = "280%", height = "1300"),
                      style = "margin-top: 80px;  margin-left: -480px;"),
            
            tableOutput('Enrichment'),
            
          ), ))
        )
      )
    )
    
    ,
    # -----------   KEGG Tab
    tabPanel(
      tags$h3("KEGG Maps"),
      value = 1,
      dropMenu(
        dropdownButton(
          "Info",
          status = 'info',
          size = "xs",
          icon = icon('info')
        ),
        h3(strong('Information')),
        br(),
        h5(
          'In this tab, the user is presented with the capability to visualize the KEGG Pathway Maps, which provide a holistic understanding of the functional and structural aspects of the biological system under investigation. This is achieved by mapping the isolated genes onto the KEGG Pathway Maps. The user is provided with a data frame containing the pathway identifiers, and can simply copy the desired identifier and press the "go" button to initiate the visualization process. It is important to note that the correct organism must be selected in the first tab prior to utilizing this feature, as the KEGG Pathway Maps are organism-specific.'
        ),
        
        
        
        
        
        placement = "bottom-start",
        arrow = TRUE,
        theme = "material",
        maxWidth = 1000
      ),
      # -----------   sidebarLayout
      sidebarLayout(
        # -------------- sidebarPanel
        sidebarPanel(
          h5("Give ID of KEGG pathways. ⚠ Don't forget to choose the right Organismus !!!"),
          textInput("inText", "Pathway iD"),
          column(
            1,
            offset = 2,
            actionButton("click2", "  Go  ",
                         style = "position: absolute; right: -130px; bottom: -72px")
          ),
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          dataTableOutput("KEGG"),
          textOutput ("text1"),
          
          sidebarPanel(
            position = "left",
            width = 12,
            
            
            plotOutput("KEGGmap", width = "100%",height = 1000) ,#, width = "100%", height = 500
            
            
          )
        )
        
        
      )
    ),
    
    # ----- Graphs Network Tab
    tabPanel(
      title = tags$h3("Graph Analysis"),
      
      dropMenu(
        dropdownButton(
          "Info",
          status = 'info',
          size = "xs",
          icon = icon('info')
        ),
        h3(strong('Information')),
        br(),
        h5(
          'On the final tab, the user is provided with the capability to perform protein-protein interaction (PPI) network analysis and similarity graph analysis. The user can specify a threshold for the PPI network, which ensures that only interactions with a combined score greater than the specified threshold are included in the network. The PPI network is based on the STRINGdb database, which is a widely used database for protein-protein interactions. The similarity graph feature allows the user to discover molecular modules in genetic networks by measuring the similarity between the profiles of gene interactions in a cell. A pearson correlation threshold is available for this analysis, allowing the user to adjust the level of similarity required for two genes to be considered as similar.'
        ),
        
        
        
        
        
        placement = "bottom-start",
        arrow = TRUE,
        theme = "material",
        maxWidth = 1000
      ),
      
      sidebarPanel(
        width = 4,
        
          column(
            12,
            p(
              "Select the number of genes to analyze based on the scoring priority of the algorithm used in the analysis.⚠ Don't forget to choose the right Organismus !!! " ,
            sliderInput(
              "Genes",
              "Number of  Genes:",
              
              min = 2,
              max = 3000,
              value = 100
            ) ) ),
         
        checkboxInput(
          "PPInetwork1",
          "Protein–protein interaction (PPI) network analysis ️🕸",
          value = FALSE
        ),
        
        numericInput(
          "Score_Threshold_PPI",
          "Conditionally load interactions based on a threshold.",
          value = 400,
          min = 50,
          max = 1000,
          step = 10
        ),
        
        checkboxInput("graph1", "Similarity Graph 📊", value = FALSE),
        numericInput(
          "Pearson_correlation",
          " Remove edges below absolute Pearson correlation",
          value = 0.5,
          min = 0.1,
          max = 0.99,
          step = 0.1
        ),
        
        
        column(
          width = 4,
          actionButton("run_button1", "Plotting Graphs", style = "position: absolute; right: -145px; bottom: -72px")
        ),
      ),
      # Main panel for displaying outputs ----
      mainPanel(
        tableOutput("text3"),
        
        
        
        column(
          width = 3,
          plotOutput("PPInetwork", width = "400%", height = "800"),
          visNetworkOutput("graph", width = "400%", height = "800")
          
        )
      )
      
    ),
    
    
  )
