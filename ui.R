# Shared upload limits and validation helpers.
source("./Scripts/InputValidation.R", local = TRUE)

#Load the necessary libraries
libs <-c("shiny",
         "shinyalert",
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
         "scran",
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
         "M3Drop",
         "ComplexHeatmap",
         "igraph",
         "visNetwork",
         "SingleR",
         "shinyjs",
         "STRINGdb",
         "fastshap",
         "xgboost",
         "randomForest"
)
lapply(libs, require, character.only = TRUE)
source("Scripts/PublicationUI.R", local = TRUE)



options(repos = BiocManager::repositories())
#options(download.file.method = "libcurl")



websiteLive <- TRUE
#Increase the size of the acceptable import dataset
options(shiny.maxRequestSize = SCGENES_MAX_UPLOAD_BYTES)





# Define UI for data upload app ----

ui <- navbarPage(
  tags$img(
    src = "scGenesElite.png",
    class = "brand-mark",
    alt = "scGeneFinder"
  ),
  id = "main_nav",
  theme = shinythemes::shinytheme("flatly"),
  windowTitle = "scGeneFinder - Single-Cell Gene Analysis Platform",
  collapsible = TRUE,
  header = tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=Space+Grotesk:wght@500;600;700&display=swap"),
    tags$link(rel = "stylesheet", href = "style.css?v=20260906-1"),
    tags$script(src = "radio-buttons.js?v=20260906-1"),
    tags$script(src = "analysis-run-control.js?v=20260906-1")
  ),
  # Home/Overview Tab ----
  tabPanel(
    title = "Home",
    value = "home",
    id = "home",
    div(
      class = "home-page",
      tags$section(
        class = "home-hero",
        div(
          class = "hero-copy",
          div(class = "eyebrow", "Single-cell gene discovery"),
          h1("Find robust marker genes from scRNA-seq data."),
          p(
            class = "hero-lede",
            "Compare more than 15 feature-selection techniques, build ensemble methods, and move from expression data to interpretable biomarkers in one guided workflow."
          ),
          div(
            class = "hero-actions",
            tags$a(
              href = "#",
              class = "btn btn-primary btn-lg sc-nav-link",
              `data-nav-target` = "upload",
              icon("upload"),
              "Upload a dataset"
            ),
            tags$a(
              href = "#workflow",
              class = "btn btn-secondary btn-lg",
              icon("diagram-project"),
              "Explore the workflow"
            )
          ),
          div(
            class = "hero-proof",
            span(icon("layer-group"), " 15+ selection methods"),
            span(icon("chart-line"), " Predictive validation"),
            span(icon("network-wired"), " Biological interpretation")
          )
        ),
        div(
          class = "hero-visual",
          div(
            class = "visual-frame",
            tags$img(
              src = "overview.jpg",
              alt = "scGeneFinder analysis workflow overview"
            )
          ),
          p(class = "visual-caption", "From uploaded counts to ranked genes, pathways and interaction networks.")
        )
      ),
      tags$section(
        class = "capabilities-section",
        div(class = "section-kicker", "What you can do"),
        h2("One workspace for gene selection and interpretation"),
        div(
          class = "capabilities-grid",
          div(
            class = "capability-card",
            div(class = "capability-icon icon-blue", icon("filter")),
            h3("Select robust features"),
            p("Compare variable-gene, statistical, machine-learning and SHAP-based methods with configurable thresholds.")
          ),
          div(
            class = "capability-card",
            div(class = "capability-icon icon-violet", icon("code-branch")),
            h3("Build ensembles"),
            p("Combine complementary techniques and rank consensus biomarkers using a unified analysis pipeline.")
          ),
          div(
            class = "capability-card",
            div(class = "capability-icon icon-teal", icon("microscope")),
            h3("Interpret biology"),
            p("Explore enrichment results, KEGG pathway maps, heatmaps and gene-interaction networks.")
          )
        )
      ),
      tags$section(
        id = "workflow",
        class = "workflow-section",
        div(class = "section-kicker", "Guided analysis"),
        h2("A clear path from data to biological insight"),
        div(
          class = "workflow-steps",
          div(class = "workflow-step", span(class = "step-number", "1"), div(h3("Upload"), p("Import CSV or RDS data and verify its structure."))),
          div(class = "workflow-step", span(class = "step-number", "2"), div(h3("Configure"), p("Choose preprocessing and gene-selection methods."))),
          div(class = "workflow-step", span(class = "step-number", "3"), div(h3("Analyze"), p("Rank genes and validate their predictive performance."))),
          div(class = "workflow-step", span(class = "step-number", "4"), div(h3("Interpret"), p("Explore pathways, enrichment and networks.")))
        )
      ),
      tags$footer(
        class = "app-footer",
        div(
          div(class = "footer-brand", "scGeneFinder"),
          p("A focused workspace for single-cell marker discovery.")
        ),
        div(
          class = "footer-contact",
          span("Questions or support?"),
          a(href = "mailto:p.paplomatas@hotmail.com", icon("envelope"), " p.paplomatas@hotmail.com")
        )
      )
    )
  ),
  tabPanel(
    title = "Data Upload",
    value = "upload",
    div(
      class = "page-header",
      div(
        div(class = "eyebrow", "Step 1 of 4"),
        h1("Upload and validate your dataset"),
        p("Choose the organism and identifier format, then upload a CSV or RDS expression matrix. The preview helps you confirm the structure before analysis.")
      ),
      downloadButton(
        "Example",
        "Download example CSV",
        class = "btn btn-secondary",
        icon = icon("download")
      )
    ),
    tags$details(
      class = "requirements-card",
      tags$summary(icon("circle-info"), " Dataset requirements"),
      div(
        class = "requirements-content",
        tags$ul(
          tags$li("Rows represent cells and expression columns represent genes."),
          tags$li("The final column must contain the class or state annotation."),
          tags$li("Both binary and multiclass annotations are supported."),
          tags$li("CSV and RDS files are accepted; normalized values are recommended.")
        ),
        p(strong("Expected shape:"), " N x E, where N is the number of cells and E - 1 is the number of genes.")
      )
    ),
    sidebarLayout(
      sidebarPanel(
        width = 4,
        div(
          class = "section-card upload-settings",
          div(class = "section-heading", div(h2("Dataset settings"), p("These choices are used by downstream annotation and pathway tools.")), icon("sliders")),
          div(
            class = "form-section",
            h3("Organism"),
            radioButtons(
              "organismus",
              label = NULL,
              choices = c("H. sapiens" = "Human", "M. musculus" = "Mouse"),
              selected = "Mouse",
              inline = TRUE
            )
          ),
          div(
            class = "form-section",
            h3("Gene identifier"),
            radioButtons(
              "GENEid",
              label = NULL,
              choices = c("Gene Symbol" = "SYMBOL", "Ensembl ID" = "EnsemblGenes", "Entrez ID" = "ENTREZID"),
              selected = "SYMBOL",
              inline = TRUE
            )
          ),
          div(
            class = "form-section upload-files",
            h3("Expression file"),
            p(class = "form-help", "Upload either one RDS file or one CSV file."),
            div(
              class = "file-field",
              div(class = "file-type", icon("file-code"), span("RDS")),
              fileInput("rdsFile", label = NULL, buttonLabel = "Choose RDS", placeholder = "No file selected", multiple = FALSE, accept = ".rds")
            ),
            div(class = "upload-separator", span("or")),
            div(
              class = "file-field",
              div(class = "file-type", icon("file-csv"), span("CSV")),
              fileInput(
                "file1",
                label = NULL,
                buttonLabel = "Choose CSV",
                placeholder = "No file selected",
                multiple = FALSE,
                accept = c("text/csv", "text/comma-separated-values,text/plain", ".csv")
              )
            )
          ),
          tags$details(
            class = "advanced-options",
            tags$summary(icon("wrench"), " CSV parsing options"),
            div(
              class = "advanced-options-content",
              radioButtons("header", "Header row", choices = list("Yes" = TRUE, "No" = FALSE), selected = TRUE, inline = TRUE),
              radioButtons("sep", "Separator", choices = c("Comma" = ",", "Semicolon" = ";", "Tab" = "\t"), selected = ",", inline = TRUE),
              radioButtons("quote", "Quote character", choices = c("None" = "", "Double quote" = '"', "Single quote" = "'"), selected = '"', inline = TRUE),
              radioButtons("disp", "Preview size", choices = c("First 10 rows" = "head", "Up to 1,000 rows" = "all"), selected = "head", inline = FALSE)
            )
          )
        )
      ),
      mainPanel(
        width = 8,
        div(
          class = "section-card preview-card",
          div(
            class = "section-heading",
            div(h2("Data preview"), p("Verify genes, annotations and table orientation before continuing.")),
            span(class = "status-chip", icon("shield"), " Validation enabled")
          ),
          conditionalPanel(
            condition = "output.uploadState !== 'ready'",
            div(
              class = "empty-state upload-empty-state",
              div(class = "empty-state-icon", icon("table-cells-large")),
              h3("Your dataset preview will appear here"),
              p("Upload a CSV or RDS file using the sidebar to preview your dataset here.")
            )
          ),
          conditionalPanel(
            condition = "output.uploadState === 'ready'",
            div(
              class = "data-preview-output",
              DT::dataTableOutput("contents"),
              DT::dataTableOutput("Rvalue")
            )
          )
        ),
        div(
          class = "next-step-card",
          div(icon("arrow-right"), div(h3("Next: configure gene selection"), p("After validating the preview, continue to Run Analysis."))),
          tags$a(href = "#", class = "btn btn-primary sc-nav-link", `data-nav-target` = "analysis", "Continue")
        )
      )
    )
  ),

    # Run Analysis Tab
    tabPanel(
      title = "Run Analysis",
      value = "analysis",
      div(
        class = "page-header analysis-page-header",
        div(
          div(class = "eyebrow", "Steps 2 and 3 of 4"),
          h1("Configure and run gene selection"),
          p("Choose preprocessing, select one method or build an ensemble, then generate ranked biomarkers and validation plots.")
        ),
        div(
          class = "page-header-help",
          icon("lightbulb"),
          span("Start with one method from a single category, then compare with an ensemble.")
        )
      ),
      div(
        class = "section-card preprocessing-card",
        div(
          class = "section-heading",
          div(h2("Preprocessing"), p("Prepare the expression matrix before feature selection.")),
          div(
            class = "info-actions",
            dropMenu(
              dropdownButton("Filtering help", status = "info", size = "sm", icon = icon("circle-info")),
              h3("Variance filtering"),
              p("The low-variance filter identifies zero-variance predictors and features with very few distinct values relative to sample frequency."),
              p("The cutoff controls how aggressively low-information genes are removed."),
              placement = "bottom-end",
              arrow = TRUE,
              theme = "material",
              maxWidth = 520
            ),
            dropMenu(
              dropdownButton("Workflow help", status = "info", size = "sm", icon = icon("circle-question")),
              h3("Analysis workflow"),
              p("Variable-gene, statistical, machine-learning and SHAP methods use different selection strategies. Ensembles combine choices made in the method tabs."),
              p("Reducing uninformative features first can improve reliability and computational performance."),
              placement = "bottom-end",
              arrow = TRUE,
              theme = "material",
              maxWidth = 560
            )
          )
        ),
        div(
          class = "preprocessing-grid",
          div(
            class = "form-block",
            h3(icon("filter"), " Variance filter"),
            radioButtons(
              "VarFilter",
              label = NULL,
              choices = c("Remove low variance" = "Strict_Filter", "No filter" = "Unselect"),
              selected = "Strict_Filter",
              inline = TRUE
            ),
            numericInput("uniqueCut", "Low-variance cutoff (%)", 15, 0, 100, width = "100%")
          ),
          div(
            class = "form-block",
            h3(icon("chart-simple"), " Normalization"),
            radioButtons(
              "Norm",
              label = NULL,
              choices = c("LogNormalize" = "Normal", "Keep uploaded values" = "No_Normal"),
              selected = "Normal",
              inline = TRUE
            ),
            p(class = "form-help", "LogNormalize uses a scaling factor derived from the expression matrix.")
          )
        )
      ),
      
      # Main layout with sidebar and main content
      sidebarLayout(
        # Sidebar Panel for Genes Selection Methods
        sidebarPanel(
          width = 4,
          div(class = "method-panel-heading", div(class = "eyebrow", "Method configuration"), h2("Gene-selection methods"), p("Choose a category and configure its parameters.")),
          tabsetPanel(
            type = "tabs",
            
            # Variable Genes Tab
            tabPanel(
              "Variable Genes",
              
              # HVGs Methods
              radioButtons(
                'VariableM',
                label = "HVGs Methods:",
                choices = list(
                  "SCMarker"       = "SCMarker",
                  "scran"       = "DUBStepR",
                  "ScPNMF"         =  "ScPNMF",
                  "VST"          = "SelfE" ,
                  "M3Drop"          = "M3Drop" ,
                  "None" = "NoMethod"
                ),
                selected = "NoMethod"
              ),
              
              conditionalPanel(
                condition = "input.VariableM != 'NoMethod'",
                h3(class = "parameter-heading", "Method parameters")
              ),
              conditionalPanel(
                condition = "input.VariableM == 'SCMarker'",
                fluidRow(
                  column(4, div("SCMarker")),
                  column(4, numericInput("geneK", "geneK", value = 20, min = 0, step = 1)),
                  column(4, numericInput("cellK", "cellK", value = 20, min = 0, step = 1))
                )
              ),
              conditionalPanel(
                condition = "input.VariableM == 'DUBStepR'",
                fluidRow(
                  column(4, div("scran")),
                  column(8, numericInput("np", "Number of genes", value = 300, min = 1, step = 1))
                )
              ),
              conditionalPanel(
                condition = "input.VariableM == 'ScPNMF'",
                fluidRow(
                  column(4, div("ScPNMF")),
                  column(8, numericInput("gM", "Number of genes", value = 300, min = 1, step = 1))
                )
              ),
              conditionalPanel(
                condition = "input.VariableM == 'M3Drop'",
                fluidRow(
                  column(4, div("M3Drop")),
                  column(4, numericInput("M3dropthreshold", "Threshold", value = 0.001, min = 0, max = 1, step = 0.01)),
                  column(4, selectInput("M3Method", "Adjustment", c("Bonferroni" = "bon", "FDR" = "fdr")))
                )
              ),
              conditionalPanel(
                condition = "input.VariableM == 'SelfE'",
                fluidRow(
                  column(4, div("VST")),
                  column(4, numericInput("n", "Number of features", value = 300, min = 1, step = 1)),
                  column(4, selectInput("distMethod", "Distance method", c("Euclidean" = "EucDist", "KL" = "KL", "DPNMF" = "DPNMF")))
                )
              )
            ),
            
            # DEGs Tab
            tabPanel(
              "Statistical",
              
              radioButtons (
                "P_method",
                "DEGs Methods:",
                c(
                  "Wilcoxon rank sum test" = "Seurat_method",
                  "Beta-Poisson generalized linear model" = "BPSC_metchod",
                  "Wald test" = "MAST_method",
                  "Likelihood Ratio Test" = "DESeq2_method",
                  "None" = "Empty"
                ),
                selected = "Empty"
              ),
              
              conditionalPanel(
                condition = "input.P_method != 'Empty'",
                numericInput(
                  "PvalueNum",
                  "P-value threshold",
                  value = 0.01,
                  min = 0.001,
                  max = 0.99,
                  step = 0.01
                ),
                numericInput(
                  "logfc",
                  "LogFC threshold",
                  value = 1,
                  min = 0,
                  max = 5,
                  step = 0.1
                )
              )
            ),
            
            # Machine Learning Feature Selection Tab
            tabPanel(
              "Machine Learning",
              
              selectInput(
                "ML_Method",
                "Tree-based ML Feature Selection Methods:",
                c(
                  "Random Forest Algorithm" = "rf",
                  "eXtreme Gradient Boosting" = "xgbTree",
                  "Bagged CART" = "treebag",
                  "Recursive Partitioning and Regression Trees" = "rpart",
                  "C5.0 Decision Trees and Rule-Based Models" = "C5.0",
                  "None" = "Empty"
                ),
                selected = "Empty"
              ),
              conditionalPanel(
                condition = "input.ML_Method != 'Empty'",
                numericInput(
                  "importanceLimit",
                  "Importance threshold",
                  value = 10,
                  min = 0,
                  max = 100,
                  step = 1
                ),
                div(class = "inline-note", icon("circle-info"), span("Tree-based methods rank genes by feature importance."))
              )
            ),
            
            # SHAP Values Tab (Separate)
            tabPanel(
              "SHAP Values",
              
              #Info Button for SHAP Values
              dropMenu(
                dropdownButton(
                  "SHAP Info",
                  status = 'info',
                  size = "xs",
                  icon = icon('info-circle')
                ),
                h3(strong('SHAP Values Information')),
                br(),
                h5(
                  "SHAP Values (SHapley Additive exPlanations) provide a unified approach to explain the output of any machine learning model. 
                They are based on game theory and offer several advantages over traditional feature importance methods:
                
                1. **Model Interpretability**: SHAP values explain how each feature contributes to the model's prediction
                2. **Feature Interactions**: They account for interactions between features
                3. **Consistency**: SHAP values are consistent across different models
                4. **Robustness**: More stable than other feature importance methods
                
                In scGeneFinder, SHAP values are calculated for Random Forest and XGBoost models to identify the most important genes for classification.
                The method splits data into training (75%) and test (25%) sets, trains the model, and calculates REAL SHAP values using Monte Carlo simulations for feature importance ranking."
                ),
                
                placement = "bottom-start",
                arrow = TRUE,
                theme = "material",
                maxWidth = 1000
              ),
              
              selectInput(
                "SHAP_Method",
                "SHAP Values Methods:",
                c(
                  "SHAP Values (Random Forest)" = "shap_rf",
                  "SHAP Values (XGBoost)" = "shap_xgb",
                  "None" = "Empty"
                ),
                selected = "Empty"
              ),
              conditionalPanel(
                condition = "input.SHAP_Method != 'Empty'",
                numericInput(
                  "SHAP_importanceLimit",
                  "SHAP importance threshold",
                  value = 0.01,
                  min = 0,
                  max = 100,
                  step = 0.001
                ),
                div(class = "inline-note", icon("lightbulb"), span("Use 0.001–0.1 for absolute SHAP values, or an integer for a top-N cutoff."))
              )
            ),
            
            # Ensemble Approach Tab
            tabPanel(
              "Ensemble",
              div(class = "alert alert-info",
                HTML('<strong>Note:</strong> The selection parameters of the methods, which were selected in the previous tabs, remain consistent for the Ensemble approach!')
              ),
              div(class = "alert alert-info",
                HTML('<strong>Tip:</strong> SHAP Values methods are now available for ensemble analysis, providing interpretable feature importance based on game theory principles.')
              ),
              
              radioButtons(
                'ensembleVar',
                label    = "HVGs Methods:",
                choices = list(
                  "SCMarker"       = "SCMarker",
                  "scran"       = "DUBStepR",
                  "ScPNMF"         =  "ScPNMF",
                  "M3Drop"          = "M3Drop" ,
                  "VST"          = "SelfE" ,
                  "None" = "NoMethod"
                ),
                selected = "NoMethod"
              ),
              
              radioButtons(
                'ensemblePvalue',
                label    = "DEGs Methods:",
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
                label    = "Machine Learning Methods:",
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
                'ensembleSHAP',
                label    = "SHAP Values Methods:",
                choices = list(
                  "SHAP Values (Random Forest)" = "shap_rf",
                  "SHAP Values (XGBoost)" = "shap_xgb",
                  "None" = "NoMethod"
                ),
                selected = "NoMethod"
              )
            )
          )
        ),
        
        # Main panel for displaying outputs
        mainPanel(
          width = 8,
          div(
            class = "section-card analysis-options",
            div(class = "section-heading", div(h2("Output settings"), p("Control the size and organization of the generated result set.")), icon("chart-column")),
            numericInput("genes", "Number of ranked genes", value = 20, min = 2, max = 500, step = 1, width = "100%"),
            div(class = "inline-note", icon("triangle-exclamation"), span("Organism and gene ID settings from the upload step are reused here."))
          ),
          div(
            class = "analysis-action-bar",
            div(
              h3("Ready to run?"),
              p(class = "analysis-run-hint is-idle", "The runtime depends on dataset size and selected methods."),
              p(class = "analysis-run-hint is-running", "Analysis in progress. Press Stop analysis to discard this run."),
              p(class = "analysis-run-hint is-stopping", "Stopping. The current step has to finish before the run is discarded.")
            ),
            div(
              class = "analysis-actions",
              selectInput(
                "downloadDpi", NULL,
                choices = c("300 dpi" = "300", "600 dpi" = "600"),
                selected = "600"
              ),
              conditionalPanel(
                condition = "!output.analysisReady",
                tags$span(class = "btn btn-secondary disabled", icon("download"), " Download results")
              ),
              conditionalPanel(
                condition = "output.analysisReady",
                downloadLink(
                  "downloadData",
                  tagList(icon("download"), " Download results"),
                  class = "btn btn-secondary",
                  title = "Download all four panels and composite as PDF/PNG, plotted values, predictions and settings"
                )
              ),
              actionButton(
                "click",
                tagList(icon("play"), " Run analysis"),
                class = "btn btn-primary analysis-run-btn",
                `aria-label` = "Run the gene selection analysis"
              ),
              actionButton(
                "stopAnalysis",
                tagList(icon("stop"), " Stop analysis"),
                class = "btn btn-danger analysis-stop-btn",
                `aria-label` = "Stop the running gene selection analysis"
              )
            )
          ),
          div(class = "analysis-status", tableOutput("text")),
          conditionalPanel(
            condition = "!output.analysisReady",
            div(
              class = "empty-state analysis-empty-state",
              div(class = "empty-state-icon", icon("flask")),
              h3("Results will appear after analysis"),
              p("Select a method, review the output settings and run the analysis to generate ranked genes, expression plots and within-dataset classification.")
            )
          ),
          conditionalPanel(
            condition = "output.analysisReady",
            div(
              class = "analysis-results",
              div(
                class = "results-grid biomarkers-grid",
                div(
                  class = "result-card result-card-wide",
                  div(class = "result-card-heading", div(h2("Ranked genes"), p("Genes ordered by the selected method's score.")), icon("ranking-star")),
                  girafeOutput("TheBarPlot", height = "460px")
                ),
                div(
                  class = "result-card result-list-card",
                  div(class = "result-card-heading", div(h2("Gene list"), p("Selected potential biomarkers.")), icon("list-ol")),
                  DT::dataTableOutput("GenesList")
                )
              ),
              publication_results_ui()
            )
          )
        )
      )
    ),
    
    # Results workspace
    navbarMenu(
      title = "Explore Results",
      icon = icon("chart-pie"),
      tabPanel(
      title = "Enrichment",
      value = "enrichment",
      div(
        class = "page-header",
        div(
          div(class = "eyebrow", "Step 4 of 4 - Biological interpretation"),
          h1("Enrichment analysis"),
          p("Test top-ranked biomarkers against pathway, ontology and disease-drug libraries from Enrichr.")
        ),
        dropMenu(
          dropdownButton("About enrichment", status = "info", size = "sm", icon = icon("circle-info")),
          h3("How enrichment works"),
          p("Choose how many of the highest-scoring genes to test, then select one or more knowledge-base categories."),
          p("Enable all ontology terms for a broader scan, or choose focused libraries for a more targeted interpretation."),
          placement = "bottom-end",
          arrow = TRUE,
          theme = "material",
          maxWidth = 560
        )
      ),
      div(
        class = "section-card enrichment-controls",
        div(class = "section-heading", div(h2("Analysis setup"), p("Select the gene set size and knowledge bases.")), icon("book-medical")),
        div(
          class = "enrichment-grid",
          numericInput("genes1", "Top-ranked genes", value = 50, min = 10, max = 9000, step = 10, width = "100%"),
          div(
            class = "form-block",
            h3("Use all ontology terms"),
            radioButtons("all", label = NULL, choices = list("Yes" = TRUE, "No" = FALSE), selected = FALSE, inline = TRUE)
          ),
          selectInput(
            "BP",
            "Biological pathway",
            c(
              "Select a library" = "-",
              "KEGG 2021 Human" = "KEGG_2021_Human",
              "WikiPathway 2021 Human" = "WikiPathway_2021_Human",
              "BioPlanet 2019" = "BioPlanet_2019",
              "BioCarta 2016" = "BioCarta_2016",
              "MSigDB Hallmark 2020" = "MSigDB_Hallmark_2020",
              "Reactome 2016" = "Reactome_2016"
            )
          ),
          selectInput(
            "BO",
            "Biological ontology",
            c(
              "Select a library" = "-",
              "GO Biological Process 2021" = "GO_Biological_Process_2021",
              "GO Molecular Function 2021" = "GO_Molecular_Function_2021",
              "GO Cellular Component 2021" = "GO_Cellular_Component_2021",
              "MGI Mammalian Phenotype Level 4 2021" = "MGI_Mammalian_Phenotype_Level_4_2021",
              "Human Phenotype Ontology" = "Human_Phenotype_Ontology",
              "Jensen DISEASES" = "Jensen_DISEASES"
            )
          ),
          selectInput(
            "DD",
            "Diseases and drugs",
            c(
              "Select a library" = "-",
              "DisGeNET" = "DisGeNET",
              "DSigDB" = "DSigDB",
              "DrugMatrix" = "DrugMatrix",
              "OMIM Disease" = "OMIM_Disease",
              "HDSigDB Human 2021" = "HDSigDB_Human_2021",
              "COVID-19 Related Gene Sets 2021" = "COVID-19_Related_Gene_Sets_2021"
            )
          )
        ),
        div(
          class = "section-action-row",
          div(class = "inline-note", icon("circle-info"), span("Enrichment uses the ranking produced in Run Analysis.")),
          actionButton("click1", tagList(icon("play"), " Run enrichment"), class = "btn btn-primary", `aria-label` = "Run the enrichment analysis")
        )
      ),
      conditionalPanel(
        condition = "!input.click1",
        div(
          class = "empty-state results-empty-state",
          div(class = "empty-state-icon", icon("chart-bar")),
          h3("No enrichment results yet"),
          p("Run gene selection first, choose at least one knowledge base, and start enrichment to populate this workspace."),
          tags$a(href = "#", class = "btn btn-secondary sc-nav-link", `data-nav-target` = "analysis", "Go to Run Analysis")
        )
      ),
      conditionalPanel(
        condition = "input.click1 > 0",
        div(
          class = "result-card result-card-full enrichment-results",
          div(
            class = "result-card-heading",
            div(h2("Enriched terms"), p("Significant pathway and ontology associations for the selected genes.")),
            div(
              class = "result-card-heading-actions",
              icon("chart-bar"),
              tags$button(
                type = "button",
                class = "fullscreen-btn",
                `data-target` = "BioBarPlot",
                `aria-label` = "View the enriched terms chart full screen",
                title = "View full screen",
                icon("up-right-and-down-left-from-center")
              ),
              selectInput(
                "enrichmentDpi", NULL,
                choices = c("300 dpi" = "300", "600 dpi" = "600"),
                selected = "600"
              ),
              downloadLink(
                "downloadEnrichment",
                tagList(icon("download"), " Download"),
                class = "btn btn-secondary btn-sm"
              )
            )
          ),
          plotOutput("BioBarPlot", width = "100%", height = "680px"),
          tableOutput("Enrichment")
        )
      )
    ),
    
    # KEGG Maps Tab
    tabPanel(
      title = "KEGG Maps",
      value = "kegg",
      div(
        class = "page-header",
        div(
          div(class = "eyebrow", "Step 4 of 4 - Pathway context"),
          h1("KEGG pathway maps"),
          p("Map selected biomarkers onto a KEGG pathway to inspect their biological context and relationships.")
        ),
        dropMenu(
          dropdownButton("How it works", status = "info", size = "sm", icon = icon("circle-info")),
          h3("Visualizing a KEGG pathway"),
          p("Find a pathway in the reference table, copy its five-digit identifier, and enter it in the Pathway ID field."),
          p("The organism selected during upload determines whether human or mouse KEGG maps are used."),
          placement = "bottom-end",
          arrow = TRUE,
          theme = "material",
          maxWidth = 520
        )
      ),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          div(
            class = "section-card kegg-controls",
            div(class = "section-heading", div(h2("Pathway selection"), p("Enter a KEGG identifier from the table.")), icon("map")),
            textInput("inText", "Pathway ID", placeholder = "e.g. 04110"),
            div(class = "inline-note", icon("dna"), span("Uses the organism selected during data upload.")),
            actionButton(
              "click2",
              tagList(icon("eye"), " Visualize pathway"),
              class = "btn btn-primary btn-block",
              `aria-label` = "Visualize KEGG pathway"
            )
          )
        ),
        mainPanel(
          width = 8,
          div(
            class = "result-card kegg-reference-card",
            div(class = "result-card-heading", div(h2("KEGG pathway reference"), p("Search by pathway name or copy a pathway identifier.")), icon("table-list")),
            dataTableOutput("KEGG")
          ),
          conditionalPanel(
            condition = "!input.click2",
            div(
              class = "empty-state kegg-empty-state",
              div(class = "empty-state-icon", icon("map-location-dot")),
              h3("Choose a pathway to create a map"),
              p("Search the reference table, enter its Pathway ID and select Visualize pathway.")
            )
          ),
          conditionalPanel(
            condition = "input.click2 > 0",
            div(
              class = "result-card kegg-map-card",
              div(
                class = "result-card-heading",
                div(h2("Pathway visualization"), p("Selected genes highlighted on the KEGG pathway map.")),
                div(
                  class = "result-card-heading-actions",
                  icon("route"),
                  tags$button(
                    type = "button",
                    class = "fullscreen-btn",
                    `data-target` = "KEGGmap",
                    `aria-label` = "View the KEGG pathway map full screen",
                    title = "View full screen",
                    icon("up-right-and-down-left-from-center")
                  ),
                  downloadLink(
                    "downloadKeggMap",
                    tagList(icon("download"), " Download"),
                    class = "btn btn-secondary btn-sm",
                    title = "Downloads the original KEGG pathway image, at its full native resolution"
                  )
                )
              ),
              plotOutput("KEGGmap", width = "100%", height = "680px")
            )
          )
        )
      )
    ),
    
    # Graph Analysis Tab
    tabPanel(
      title = "Networks",
      value = "graphs",
      div(
        class = "page-header",
        div(
          div(class = "eyebrow", "Step 4 of 4 - Network interpretation"),
          h1("Gene interaction networks"),
          p("Generate protein-protein interaction and expression-similarity networks from the highest-ranked genes.")
        ),
        dropMenu(
          dropdownButton("About networks", status = "info", size = "sm", icon = icon("circle-info")),
          h3("Network options"),
          p("PPI analysis uses STRINGdb combined scores. Higher thresholds retain stronger evidence-backed interactions."),
          p("The similarity network links genes whose expression profiles exceed the selected absolute Pearson correlation."),
          placement = "bottom-end",
          arrow = TRUE,
          theme = "material",
          maxWidth = 560
        )
      ),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          div(
            class = "section-card network-controls",
            div(class = "section-heading", div(h2("Network setup"), p("Choose the scope and network types.")), icon("share-nodes")),
            numericInput("Genes", "Top-ranked genes", value = 50, min = 2, max = 3000, step = 10, width = "100%"),
            div(
              class = "network-option",
              div(h3("Protein-protein interactions"), p("Connect genes using STRINGdb evidence scores.")),
              radioButtons("PPInetwork1", label = NULL, choices = list("Generate" = TRUE, "Skip" = FALSE), selected = FALSE, inline = TRUE),
              numericInput("Score_Threshold_PPI", "Minimum STRING score", value = 400, min = 50, max = 1000, step = 10, width = "100%")
            ),
            div(
              class = "network-option",
              div(h3("Expression similarity"), p("Connect genes with similar expression profiles.")),
              radioButtons("graph1", label = NULL, choices = list("Generate" = TRUE, "Skip" = FALSE), selected = FALSE, inline = TRUE),
              numericInput("Pearson_correlation", "Minimum absolute Pearson correlation", value = 0.5, min = 0.1, max = 0.99, step = 0.1, width = "100%")
            ),
            actionButton(
              "run_button1",
              tagList(icon("share-nodes"), " Generate networks"),
              class = "btn btn-primary btn-block",
              `aria-label` = "Generate PPI and similarity graphs"
            )
          )
        ),
        mainPanel(
          width = 8,
          conditionalPanel(
            condition = "!input.run_button1",
            div(
              class = "empty-state network-empty-state",
              div(class = "empty-state-icon", icon("circle-nodes")),
              h3("Configure a network to begin"),
              p("Select at least one network type, adjust its threshold and generate the visualization.")
            )
          ),
          conditionalPanel(
            condition = "input.run_button1 > 0 && input.PPInetwork1 == 'FALSE' && input.graph1 == 'FALSE'",
            div(class = "alert alert-warning", icon("triangle-exclamation"), strong(" Select at least one network type before generating results."))
          ),
          conditionalPanel(
            condition = "input.run_button1 > 0 && input.PPInetwork1 == 'TRUE'",
            div(
              class = "result-card network-result-card",
              div(class = "result-card-heading", div(h2("Protein-protein interaction network"), p("STRINGdb interactions above the selected score threshold.")), icon("diagram-project")),
              plotOutput("PPInetwork", width = "100%", height = "660px")
            )
          ),
          conditionalPanel(
            condition = "input.run_button1 > 0 && input.graph1 == 'TRUE'",
            div(
              class = "result-card network-result-card",
              div(class = "result-card-heading", div(h2("Expression-similarity network"), p("Interactive gene modules based on expression correlation.")), icon("circle-nodes")),
              visNetworkOutput("graph", width = "100%", height = "660px")
            )
          )
        )
      )
    )
  )
  )
