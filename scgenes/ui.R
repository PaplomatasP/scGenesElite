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
         "shinycustomloader",
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



options(repos = BiocManager::repositories())
#options(download.file.method = "libcurl")



websiteLive <- TRUE
#Increase the size of the acceptable import dataset
options(shiny.maxRequestSize = 256 * 2048 ^ 2)





# Define UI for data upload app ----

ui <- navbarPage(
  tags$img(src='scGenesElite.png', width = '190px', height = '80px', 
           style = "margin-top: -30px; border-radius: 15px; position: relative; z-index: 1000; display: block !important; visibility: visible !important; opacity: 1 !important;"
  ),
  theme = shinythemes::shinytheme("flatly"),
  windowTitle = "scGeneFinder - Single-Cell Gene Analysis Platform",
  collapsible = TRUE,
  header = tags$head(
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&family=Space+Grotesk:wght@500;600;700&display=swap"),
    tags$style(HTML("

      /* GC Petros S 16: Additional aggressive CSS for material theme dropMenu */
      /* Target all possible dropMenu/popover containers with highest specificity */
      div[class*=\"dropdown\"], div[class*=\"popover\"], div[class*=\"material\"],
      [data-theme=\"material\"], .material, .bs-popover, .popover-content {
        background-color: #ffffff !important;
        background: #ffffff !important;
        color: #1e293b !important;
      }
      
      /* Force all text to be dark */
      div[class*=\"dropdown\"] *, div[class*=\"popover\"] *, div[class*=\"material\"] *,
      [data-theme=\"material\"] *, .material *, .bs-popover *, .popover-content * {
        color: #1e293b !important;
      }
      
      /* Ensure headings are bold and dark */
      div[class*=\"dropdown\"] h3, div[class*=\"dropdown\"] h4, div[class*=\"dropdown\"] h5,
      div[class*=\"popover\"] h3, div[class*=\"popover\"] h4, div[class*=\"popover\"] h5,
      [data-theme=\"material\"] h3, [data-theme=\"material\"] h4, [data-theme=\"material\"] h5 {
        color: #1e293b !important;
        font-weight: bold !important;
      }
      
      /* Force divs with inline color styles to be visible */
      div[style*='color:blue'], div[style*='color: #3b82f6'] {
        color: #3b82f6 !important;
        font-weight: 600 !important;
      }
            /* ========================================
         GLOBAL STYLES & FOUNDATION
         ======================================== */
      * {
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
        box-sizing: border-box;
      }
      
      body {
        background: linear-gradient(135deg, #0f172a 0%, #1e293b 50%, #334155 100%);
        background-attachment: fixed;
        min-height: 100vh;
        position: relative;
      }
      
      body::before {
        content: '';
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background: 
          radial-gradient(circle at 20% 50%, rgba(59, 130, 246, 0.1) 0%, transparent 50%),
          radial-gradient(circle at 80% 80%, rgba(139, 92, 246, 0.1) 0%, transparent 50%);
        pointer-events: none;
        z-index: 0;
      }
      
      /* ========================================
         NAVIGATION BAR
         ======================================== */
      .navbar {
        background: linear-gradient(135deg, rgba(15, 23, 42, 0.95) 0%, rgba(30, 41, 59, 0.95) 100%) !important;
        backdrop-filter: blur(10px);
        border-bottom: 1px solid rgba(59, 130, 246, 0.2);
        box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4), 0 2px 8px rgba(59, 130, 246, 0.2);
        padding: 15px 0;
      }
      
      .navbar-default .navbar-nav > li > a {
        color: #e2e8f0 !important;
        font-weight: 600;
        font-size: 15px;
        padding: 12px 24px;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        border-radius: 10px;
        margin: 0 6px;
        position: relative;
        overflow: hidden;
      }
      
      .navbar-default .navbar-nav > li > a::before {
        content: '';
        position: absolute;
        top: 0;
        left: -100%;
        width: 100%;
        height: 100%;
        background: linear-gradient(90deg, transparent, rgba(59, 130, 246, 0.3), transparent);
        transition: left 0.5s;
      }
      
      .navbar-default .navbar-nav > li > a:hover::before {
        left: 100%;
      }
      
      .navbar-default .navbar-nav > li > a:hover {
        background: linear-gradient(135deg, rgba(59, 130, 246, 0.2) 0%, rgba(139, 92, 246, 0.2) 100%) !important;
        transform: translateY(-2px);
        color: #ffffff !important;
      }
      
      .navbar-default .navbar-nav > .active > a {
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%) !important;
        color: #ffffff !important;
        box-shadow: 0 4px 15px rgba(59, 130, 246, 0.4);
      }
      
      /* ========================================
         CONTENT CONTAINERS
         ======================================== */
      .tab-content {
        background: rgba(255, 255, 255, 0.98);
        border-radius: 20px;
        box-shadow: 
          0 20px 60px rgba(0, 0, 0, 0.3),
          0 0 0 1px rgba(255, 255, 255, 0.1);
        padding: 30px 20px 30px 15px;
        margin: 20px 15px 20px 10px;
        min-height: calc(100vh - 200px);
        position: relative;
      }
      
      /* Move all content more to the left */
      .tab-content > * {
        margin-left: -5px;
      }
      
      /* Move all content more to the left - Global adjustments */
      .well, .info-card {
        padding-left: 10px !important;
        margin-left: 0 !important;
      }
      
      .sidebar-panel, .main-panel {
        padding-left: 10px !important;
        margin-left: 0 !important;
      }
      
      .form-group, .form-control, .selectize-input {
        margin-left: 0 !important;
        padding-left: 8px !important;
      }
      
      ul, ol {
        padding-left: 10px !important;
        margin-left: 0 !important;
      }
      
      .radio, .checkbox {
        margin-left: 0 !important;
        padding-left: 0 !important;
      }
      
      .radio label, .checkbox label {
        padding-left: 3px !important;
      }
      
      /* Sidebar layout adjustments */
      .sidebar-layout .sidebar-panel {
        padding-left: 10px !important;
      }
      
      .sidebar-layout .main-panel {
        padding-left: 15px !important;
      }
      
      /* Column adjustments */
      [class*='col-'] {
        padding-left: 8px !important;
        padding-right: 8px !important;
      }
      
      /* Fluid row adjustments */
      .row {
        margin-left: -5px !important;
        margin-right: -5px !important;
      }
      
      /* Home Page Layout - Equal Height Columns */
      #home .row {
        display: flex;
        align-items: stretch;
      }
      
      #home .row > [class*='col-'] {
        display: flex;
        flex-direction: column;
      }
      
      #home .row > [class*='col-'] > div {
        flex: 1;
        display: flex;
        flex-direction: column;
      }
      
      /* ========================================
         TYPOGRAPHY
         ======================================== */
      h1 {
        text-align: center;
        font-size: 42px;
        font-weight: 800;
        font-family: 'Space Grotesk', sans-serif !important;
        margin-bottom: 25px;
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 50%, #ec4899 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        text-shadow: 0 4px 20px rgba(59, 130, 246, 0.3);
        letter-spacing: -0.5px;
      }
      
      h2 {
        font-size: 28px;
        font-weight: 700;
        color: #1e293b;
        margin-bottom: 20px;
        font-family: 'Space Grotesk', sans-serif !important;
      }
      
      h3 {
        font-size: 24px;
        font-weight: 700;
        background: linear-gradient(135deg, #1e293b 0%, #475569 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        margin-bottom: 18px;
        font-family: 'Space Grotesk', sans-serif !important;
      }
      
      h4 {
        font-size: 20px;
        font-weight: 600;
        color: #334155;
        margin-bottom: 15px;
      }
      
      h5 {
        font-size: 16px;
        font-weight: 600;
        color: #475569;
        line-height: 1.6;
        margin-bottom: 12px;
      }
      
      p {
        font-size: 15px;
        line-height: 1.8;
        color: #64748b;
      }
      
      strong {
        font-weight: 600;
        color: #1e293b;
      }
      
      /* ========================================
         PANELS & CARDS
         ======================================== */
      .well {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
        border: 1px solid #e2e8f0;
        border-radius: 16px;
        box-shadow: 
          0 4px 20px rgba(15, 23, 42, 0.08),
          0 0 0 1px rgba(203, 213, 225, 0.5);
        padding: 30px;
        margin-bottom: 25px;
        transition: all 0.3s ease;
      }
      
      .well:hover {
        box-shadow: 
          0 8px 30px rgba(15, 23, 42, 0.12),
          0 0 0 1px rgba(59, 130, 246, 0.3);
        border-color: rgba(59, 130, 246, 0.3);
      }
      
      /* ========================================
         FORM CONTROLS
         ======================================== */
      .form-group label {
        font-weight: 600;
        color: #1e293b;
        margin-bottom: 10px;
        font-size: 14px;
        letter-spacing: 0.3px;
      }
      
      .form-control, .selectize-input {
        border-radius: 12px;
        border: 2px solid #cbd5e1;
        padding: 12px 16px;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        font-size: 14px;
        background: #ffffff;
        box-shadow: 0 2px 8px rgba(15, 23, 42, 0.04);
      }
      
      .form-control:focus, .selectize-input.focus {
        border-color: #3b82f6;
        box-shadow: 
          0 0 0 4px rgba(59, 130, 246, 0.1),
          0 4px 12px rgba(59, 130, 246, 0.15);
        outline: none;
        transform: translateY(-1px);
      }
      
      /* ========================================
         RADIO BUTTONS & CHECKBOXES
         ======================================== */
      .radio, .checkbox {
        margin-top: 12px;
        margin-bottom: 12px;
      }
      
      input[type='checkbox'] {
        display: none !important;
      }
      
      .radio label {
        display: inline-block;
        padding: 10px 24px;
        margin: 0 8px 8px 0;
        background: linear-gradient(135deg, #f8fafc 0%, #f1f5f9 100%);
        border-radius: 12px;
        cursor: pointer;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        border: 2px solid #e2e8f0;
        min-width: 80px;
        text-align: center;
        color: #475569;
        font-weight: 500;
        box-shadow: 0 2px 8px rgba(15, 23, 42, 0.04);
        position: relative;
      }
      
      .radio label:hover {
        background: linear-gradient(135deg, #e2e8f0 0%, #cbd5e1 100%);
        border-color: #3b82f6;
        transform: translateY(-2px);
        box-shadow: 0 4px 12px rgba(59, 130, 246, 0.15);
      }
      
      /* Selected radio button styling - multiple selectors for Shiny compatibility */
      .radio input[type='radio']:checked + label,
      .radio label:has(input[type='radio']:checked),
      .radio input[type='radio']:checked ~ label {
        background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%) !important;
        color: #ffffff !important;
        border-color: #3b82f6 !important;
        box-shadow: 
          0 6px 20px rgba(59, 130, 246, 0.35),
          0 0 0 4px rgba(59, 130, 246, 0.1),
          inset 0 2px 4px rgba(0, 0, 0, 0.1) !important;
        transform: translateY(-2px);
        font-weight: 600 !important;
      }
      
      /* Alternative approach: style the parent div when input is checked */
      .radio:has(input[type='radio']:checked) label {
        background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%) !important;
        color: #ffffff !important;
        border-color: #3b82f6 !important;
        box-shadow: 
          0 6px 20px rgba(59, 130, 246, 0.35),
          0 0 0 4px rgba(59, 130, 246, 0.1),
          inset 0 2px 4px rgba(0, 0, 0, 0.1) !important;
        transform: translateY(-2px);
        font-weight: 600 !important;
      }
      
      .radio input[type='radio'] {
        display: none;
      }
      
      /* Additional visual indicator for checked state */
      .radio input[type='radio']:checked {
        display: none;
      }
      
      /* Fallback: Use JavaScript to add class for better browser compatibility */
      .radio.checked label {
        background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%) !important;
        color: #ffffff !important;
        border-color: #3b82f6 !important;
        box-shadow: 
          0 6px 20px rgba(59, 130, 246, 0.35),
          0 0 0 4px rgba(59, 130, 246, 0.1),
          inset 0 2px 4px rgba(0, 0, 0, 0.1) !important;
        transform: translateY(-2px);
        font-weight: 600 !important;
      }
      
      /* ========================================
         BUTTONS
         ======================================== */
      .btn {
        border-radius: 12px;
        padding: 14px 28px;
        font-weight: 600;
        font-size: 15px;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        border: none;
        box-shadow: 0 4px 15px rgba(0,0,0,0.15);
        letter-spacing: 0.3px;
        position: relative;
        overflow: hidden;
      }
      
      .btn::before {
        content: '';
        position: absolute;
        top: 50%;
        left: 50%;
        width: 0;
        height: 0;
        border-radius: 50%;
        background: rgba(255, 255, 255, 0.3);
        transform: translate(-50%, -50%);
        transition: width 0.6s, height 0.6s;
      }
      
      .btn:hover::before {
        width: 300px;
        height: 300px;
      }
      
      .btn-primary {
        background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%);
        color: #ffffff;
      }
      
      .btn-primary:hover {
        transform: translateY(-3px);
        box-shadow: 0 8px 25px rgba(59, 130, 246, 0.4);
        background: linear-gradient(135deg, #2563eb 0%, #1d4ed8 100%);
      }
      
      .btn-info {
        background: linear-gradient(135deg, #14b8a6 0%, #0d9488 100%);
        color: #ffffff;
      }
      
      .btn-info:hover {
        transform: translateY(-3px);
        box-shadow: 0 8px 25px rgba(20, 184, 166, 0.4);
        background: linear-gradient(135deg, #0d9488 0%, #0f766e 100%);
      }
      
      .btn-success {
        background: linear-gradient(135deg, #10b981 0%, #059669 100%);
        color: #ffffff;
      }
      
      .btn-success:hover {
        transform: translateY(-3px);
        box-shadow: 0 8px 25px rgba(16, 185, 129, 0.4);
        background: linear-gradient(135deg, #059669 0%, #047857 100%);
      }
      
      .btn-danger {
        background: linear-gradient(135deg, #ef4444 0%, #dc2626 100%);
        color: #ffffff;
      }
      
      .btn-danger:hover {
        transform: translateY(-3px);
        box-shadow: 0 8px 25px rgba(239, 68, 68, 0.4);
        background: linear-gradient(135deg, #dc2626 0%, #b91c1c 100%);
      }
      
      /* Special Action Buttons */
      #click, #click1, #click2, #run_button1 {
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%);
        color: #ffffff;
        font-size: 17px;
        font-weight: 700;
        padding: 16px 40px;
        border-radius: 14px;
        box-shadow: 
          0 8px 25px rgba(59, 130, 246, 0.35),
          0 0 0 3px rgba(59, 130, 246, 0.1);
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
      }
      
      #click:hover, #click1:hover, #click2:hover, #run_button1:hover {
        transform: translateY(-4px);
        box-shadow: 
          0 12px 35px rgba(59, 130, 246, 0.5),
          0 0 0 5px rgba(59, 130, 246, 0.15);
        background: linear-gradient(135deg, #2563eb 0%, #7c3aed 100%);
      }
      
      #downloadData, #Example {
        background: linear-gradient(135deg, #10b981 0%, #059669 100%);
        color: #ffffff;
        padding: 14px 28px;
        border-radius: 12px;
        text-decoration: none;
        font-weight: 600;
        box-shadow: 0 6px 20px rgba(16, 185, 129, 0.3);
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        display: inline-block;
      }
      
      #downloadData:hover, #Example:hover {
        transform: translateY(-3px);
        box-shadow: 0 8px 25px rgba(16, 185, 129, 0.45);
        background: linear-gradient(135deg, #059669 0%, #047857 100%);
        text-decoration: none;
      }
      
      /* ========================================
         SLIDERS
         ======================================== */
      .irs-bar {
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%);
        border: none;
        height: 6px;
      }
      
      .irs-line {
        background: #e2e8f0;
        border: none;
        height: 6px;
      }
      
      .irs-handle {
        border: 3px solid #3b82f6;
        background: #ffffff;
        box-shadow: 0 4px 12px rgba(59, 130, 246, 0.4);
        width: 20px;
        height: 20px;
        top: 22px;
      }
      
      .irs-handle:hover {
        background: #3b82f6;
      }
      
      .irs-single, .irs-from, .irs-to {
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%);
        color: #ffffff;
        font-weight: 600;
      }
      
      /* ========================================
         TABLES
         ======================================== */
      table {
        border-radius: 12px;
        overflow: hidden;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        border: 1px solid #e2e8f0;
      }
      
      .table thead {
        background: linear-gradient(135deg, #1e293b 0%, #334155 100%);
        color: #ffffff;
      }
      
      .table thead th {
        border: none;
        padding: 16px;
        font-weight: 700;
        letter-spacing: 0.5px;
        text-transform: uppercase;
        font-size: 13px;
      }
      
      .table tbody tr {
        transition: all 0.2s ease;
        border-bottom: 1px solid #f1f5f9;
      }
      
      .table tbody tr:hover {
        background: linear-gradient(135deg, #f8fafc 0%, #f1f5f9 100%);
        transform: scale(1.01);
      }
      
      .table tbody td {
        padding: 14px 16px;
        color: #475569;
      }
      
      /* ========================================
         INFO CARDS
         ======================================== */
      .info-card {
        background: linear-gradient(135deg, #ffffff 0%, #f8fafc 100%);
        border-radius: 16px;
        padding: 30px 30px 30px 15px;
        margin-bottom: 25px;
        box-shadow: 
          0 8px 30px rgba(15, 23, 42, 0.1),
          0 0 0 1px rgba(203, 213, 225, 0.3);
        border-left: 5px solid transparent;
        border-image: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%) 1;
        transition: all 0.3s ease;
      }
      
      .info-card:hover {
        transform: translateY(-5px);
        box-shadow: 
          0 12px 40px rgba(15, 23, 42, 0.15),
          0 0 0 1px rgba(59, 130, 246, 0.3);
      }
      
      /* ========================================
         FOOTER
         ======================================== */
      .footer {
        background: linear-gradient(135deg, rgba(248, 250, 252, 0.95) 0%, rgba(241, 245, 249, 0.95) 100%);
        backdrop-filter: blur(10px);
        color: #475569;
        text-align: center;
        padding: 35px;
        border-radius: 16px;
        margin-top: 40px;
        box-shadow: 
          0 8px 30px rgba(15, 23, 42, 0.1),
          0 0 0 1px rgba(203, 213, 225, 0.3);
        border: 1px solid rgba(226, 232, 240, 0.5);
      }
      
      .footer h5 {
        color: #1e293b;
        font-weight: 700;
        margin-bottom: 18px;
        font-size: 20px;
      }
      
      .footer p {
        color: #64748b;
        margin-bottom: 12px;
      }
      
      .footer a {
        color: #3b82f6;
        text-decoration: none;
        font-weight: 600;
        transition: all 0.3s ease;
        font-size: 18px;
        position: relative;
      }
      
      .footer a::after {
        content: '';
        position: absolute;
        bottom: -2px;
        left: 0;
        width: 0;
        height: 2px;
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%);
        transition: width 0.3s ease;
      }
      
      .footer a:hover::after {
        width: 100%;
      }
      
      .footer a:hover {
        color: #2563eb;
      }
      
      /* ========================================
         ALERTS & NOTIFICATIONS
         ======================================== */
      .alert {
        border-radius: 12px;
        border: none;
        padding: 18px 24px;
        margin-bottom: 25px;
        box-shadow: 0 4px 15px rgba(0,0,0,0.08);
        font-weight: 500;
      }
      
      .alert-info {
        background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%);
        color: #1e40af;
        border-left: 5px solid #3b82f6;
      }
      
      .alert-warning {
        background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%);
        color: #92400e;
        border-left: 5px solid #f59e0b;
      }
      
      .alert-success {
        background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%);
        color: #065f46;
        border-left: 5px solid #10b981;
      }
      
      /* ========================================
         DROPDOWN MENUS
         ======================================== */
      .dropdown-menu {
        border-radius: 12px;
        box-shadow: 
          0 10px 40px rgba(0,0,0,0.2),
          0 0 0 1px rgba(203, 213, 225, 0.3);
        border: none;
        padding: 20px;
        background: rgba(255, 255, 255, 0.98);
        backdrop-filter: blur(10px);
      /* GC Petros S 15: Fix info modal/dropMenu text color for better visibility */
      .bslib-dropdown, .shiny-material-dropdown {
        background-color: #ffffff !important;
        color: #1e293b !important;
      }
      
      .bslib-dropdown h3, .bslib-dropdown h4, .bslib-dropdown h5,
      .shiny-material-dropdown h3, .shiny-material-dropdown h4, .shiny-material-dropdown h5 {
        color: #1e293b !important;
        font-weight: bold;
      }
      
      .bslib-dropdown p, .bslib-dropdown div, .bslib-dropdown span,
      .shiny-material-dropdown p, .shiny-material-dropdown div, .shiny-material-dropdown span {
        color: #334155 !important;
      }
      
      .bslib-dropdown strong, .shiny-material-dropdown strong {
        color: #1e293b !important;
      }
      
      /* Fix blue text color in info modals - make it visible */
      .bslib-dropdown div[style*='color:blue'], 
      .shiny-material-dropdown div[style*='color:blue'],
      div[style*='color:blue'] {
        color: #3b82f6 !important;
        font-weight: 600 !important;
      }
      

      /* GC Petros S 16: Additional aggressive CSS for material theme dropMenu */
      div[class*=\"dropdown\"], div[class*=\"popover\"], div[class*=\"material\"],
      [data-theme=\"material\"], .material, .bs-popover, .popover-content {
        background-color: #ffffff !important;
        background: #ffffff !important;
        color: #1e293b !important;
      }
      
      div[class*=\"dropdown\"] *, div[class*=\"popover\"] *, div[class*=\"material\"] *,
      [data-theme=\"material\"] *, .material *, .bs-popover *, .popover-content * {
        color: #1e293b !important;
      }
      
      div[class*=\"dropdown\"] h3, div[class*=\"dropdown\"] h4, div[class*=\"dropdown\"] h5,
      div[class*=\"popover\"] h3, div[class*=\"popover\"] h4, div[class*=\"popover\"] h5,
      [data-theme=\"material\"] h3, [data-theme=\"material\"] h4, [data-theme=\"material\"] h5 {
        color: #1e293b !important;
        font-weight: bold !important;
      }
      
      div[style*='color:blue'], div[style*='color: #3b82f6'] {
        color: #3b82f6 !important;
        font-weight: 600 !important;
      }
            /* ========================================
         PROGRESS BARS
         ======================================== */
      .progress {
        height: 10px;
        border-radius: 10px;
        background: #e2e8f0;
        box-shadow: inset 0 2px 5px rgba(0,0,0,0.08);
        overflow: hidden;
      }
      
      .progress-bar {
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%);
        border-radius: 10px;
        box-shadow: 0 0 10px rgba(59, 130, 246, 0.5);
        animation: progressShine 2s ease-in-out infinite;
      }
      
      @keyframes progressShine {
        0% { background-position: -200% center; }
        100% { background-position: 200% center;
      }
      
      /* ========================================
         PLOT CONTAINERS
         ======================================== */
      .shiny-plot-output {
        border-radius: 16px;
        box-shadow: 
          0 8px 30px rgba(0,0,0,0.1),
          0 0 0 1px rgba(203, 213, 225, 0.3);
        background: #ffffff;
        padding: 25px;
        margin: 25px 0;
        transition: all 0.3s ease;
      }
      
      .shiny-plot-output:hover {
        box-shadow: 
          0 12px 40px rgba(0,0,0,0.15),
          0 0 0 1px rgba(59, 130, 246, 0.3);
      }
      
      /* ========================================
         SCROLLBAR STYLING
         ======================================== */
      ::-webkit-scrollbar {
        width: 12px;
        height: 12px;
      }
      
      ::-webkit-scrollbar-track {
        background: #f1f5f9;
        border-radius: 10px;
      }
      
      ::-webkit-scrollbar-thumb {
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%);
        border-radius: 10px;
        border: 2px solid #f1f5f9;
      }
      
      ::-webkit-scrollbar-thumb:hover {
        background: linear-gradient(135deg, #2563eb 0%, #7c3aed 100%);
      }
      
      /* ========================================
         BADGES & LABELS
         ======================================== */
      .badge {
        background: linear-gradient(135deg, #3b82f6 0%, #8b5cf6 100%);
        color: #ffffff;
        padding: 6px 14px;
        border-radius: 8px;
        font-weight: 600;
        font-size: 13px;
        box-shadow: 0 2px 8px rgba(59, 130, 246, 0.3);
      }
      
      /* ========================================
         FILE INPUT
         ======================================== */
      .shiny-input-container {
        margin-bottom: 25px;
      }
      
      .btn-file {
        background: linear-gradient(135deg, #14b8a6 0%, #0d9488 100%);
        color: #ffffff;
        border-radius: 12px;
        padding: 12px 24px;
        font-weight: 600;
        transition: all 0.3s ease;
      }
      
      .btn-file:hover {
        background: linear-gradient(135deg, #0d9488 0%, #0f766e 100%);
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(20, 184, 166, 0.4);
      }
      
      /* ========================================
         LIST STYLING
         ======================================== */
      ul {
        list-style: none !important;
        padding-left: 0 !important;
      }
      
      ul li {
        list-style: none !important;
        padding-left: 0 !important;
        margin-bottom: 12px;
        position: relative;
        padding-left: 8px;
      }
      
      ul li::before {
        content: '▸';
        position: absolute;
        left: -10px;
        color: #3b82f6;
        font-weight: bold;
        font-size: 18px;
      }
      
      /* Extra left alignment for info-card lists */
      .info-card ul li {
        padding-left: 5px;
      }
      
      .info-card ul li::before {
        left: -15px;
      }
      
      /* ========================================
         ANIMATIONS
         ======================================== */
      .fade-in {
        animation: fadeInUp 0.6s cubic-bezier(0.4, 0, 0.2, 1);
      }
      
      @keyframes fadeInUp {
        from {
          opacity: 0;
          transform: translateY(30px);
        }
        to {
          opacity: 1;
          transform: translateY(0);
        }
      }
      
      .card-hover {
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
      }
      
      .card-hover:hover {
        transform: translateY(-8px) scale(1.02);
        box-shadow: 0 20px 50px rgba(0,0,0,0.15);
      }
      
      /* ========================================
         ERROR MESSAGES
         ======================================== */
      .shiny-output-error-validation {
        color: #dc2626;
        font-weight: 600;
        padding: 15px 20px;
        background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%);
        border-radius: 12px;
        border-left: 5px solid #dc2626;
        box-shadow: 0 4px 15px rgba(220, 38, 38, 0.15);
      }
      
      /* ========================================
         RESPONSIVE DESIGN
         ======================================== */
      @media (max-width: 768px) {
        .tab-content {
          margin: 15px;
          padding: 20px;
        }
        
        h1 {
          font-size: 28px;
        }
        
        .navbar-default .navbar-nav > li > a {
          padding: 10px 16px;
          font-size: 14px;
        }
        
        .btn {
          padding: 12px 20px;
          font-size: 14px;
        }
        
        /* Home page responsive */
        #home .row {
          flex-direction: column;
        }
        
        #home .row > [class*='col-'] {
          width: 100% !important;
          min-height: auto !important;
      }
      
      /* ========================================
         SPECIAL EFFECTS
         ======================================== */
      @keyframes shimmer {
        0% {
          background-position: -1000px 0;
        }
        100% {
          background-position: 1000px 0;
        }
      }
      
      .shimmer {
        background: linear-gradient(
          90deg,
          rgba(255, 255, 255, 0) 0%,
          rgba(255, 255, 255, 0.3) 50%,
          rgba(255, 255, 255, 0) 100%
        );
        background-size: 1000px 100%;
        animation: shimmer 2s infinite;
      }
    ")),
    tags$script(HTML("
      // Function to update radio button checked state styling
      function updateRadioButtons() {
        // Find all radio buttons
        document.querySelectorAll('input[type=\"radio\"]').forEach(function(radio) {
          var radioDiv = radio.closest('.radio');
          if (radioDiv) {
            if (radio.checked) {
              radioDiv.classList.add('checked');
            } else {
              radioDiv.classList.remove('checked');
            }
          }
        });
      }
      
      // Update on page load
      document.addEventListener('DOMContentLoaded', function() {
        updateRadioButtons();
      });
      
      // Update when Shiny is ready
      if (window.Shiny) {
        $(document).on('shiny:connected', function() {
          updateRadioButtons();
        });
      }
      
      // Update on any radio button change
      $(document).on('change', 'input[type=\"radio\"]', function() {
        updateRadioButtons();
      });
      
      // Use MutationObserver to catch dynamically added radio buttons
      var observer = new MutationObserver(function(mutations) {
        updateRadioButtons();
      });
      
      // Start observing when document is ready
      $(document).ready(function() {
        observer.observe(document.body, {
          childList: true,
          subtree: true
        });
        updateRadioButtons();
      });
    "))
  ),
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
  # Home/Overview Tab ----
  tabPanel(
    title = "Home",
    value = "home",
    id = "home",
    fluidRow(
      column(7,
        div(style = "display: flex; flex-direction: column; height: 100%; min-height: 600px; padding: 20px;",
          strong(
            h1(
              "A web platform that facilitates the identification of leading genes from scRNA-seq data."
            )
          ),
          tags$hr(),
          br(),
          p(
            "Single-cell RNA-sequencing has transformed biomedical research, yet it faces computational analysis challenges.
            Navigating the vast data dimensions poses several issues, with gene selection methods being paramount. 
            Our platform, scGeneFinder, adeptly pinpoints dominant genes within scRNA-seq datasets, integrating over 15 tailored feature extraction techniques. 
            A standout feature of this app is its ensemble approach, which empowers users to craft their unique method by amalgamating one technique from each category. This tool offers an in-depth analysis of top genes, gauging their predictive accuracy and their association with various biological and drug-related ontologies. Additionally, visual aids like KEGG pathways and PPI networks offer a holistic perspective on the identified genes. 
            Given its wide-ranging functionalities, scGeneFinder stands as a comprehensive guide for uncovering and understanding transcriptional markers for intricate diseases through scRNA-seq research."
          ),
          div(style = "margin-top: auto;",
            div(
              class = "footer",
              h5("Contact Us"),
              p(
                a(href = "mailto:p.paplomatas@hotmail.com", 
                  "p.paplomatas@hotmail.com", 
                  style = "color: #3b82f6; font-size: 18px; font-weight: 500;"),
                br(),
                "Feel free to contact us for any questions or support"
              )
            )
          )
        )
      ),
      column(5,
        div(style = "display: flex; align-items: center; justify-content: center; height: 100%; min-height: 600px; padding: 20px;",
          tags$img(
            src = "overview.jpg",
            style = "width: 100%; height: auto; max-height: 100%; border-radius: 12px; box-shadow: 0 4px 15px rgba(0,0,0,0.1); object-fit: contain;"
          )
        )
      )
    )
  ),
  tabPanel(
    title = "Data Upload",
    value = "upload",
    tags$head(tags$style(
      HTML("
          .shiny-output-error-validation {
            color: #dc2626;
            font-weight: 600;
            padding: 15px 20px;
            background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%);
            border-radius: 12px;
            border-left: 5px solid #dc2626;
            box-shadow: 0 4px 15px rgba(220, 38, 38, 0.15);
          }
        ")
    )),
    fluidRow(
      column(12,
             div(class = "info-card",
                 style = "background: linear-gradient(135deg, rgba(59, 130, 246, 0.05) 0%, rgba(139, 92, 246, 0.05) 100%);",
                 h3(style = "color: #1e293b;", "📤 Data Upload Information"),
                 p(strong("Data Format Requirements:")),
                 tags$ul(
                   tags$li("Single-cell RNA-sequencing data in the form of read counts with annotations"),
                   tags$li("Matrix format: N×E, where N = cell samples and E-1 = gene expressions"),
                   tags$li("X(i,j) represents the expression value of gene j for cell i"),
                   tags$li("The last column must contain annotations (e.g., health-disease)"),
                   tags$li("Supports binary classification (control-state) or multiple classes"),
                   tags$li("Works best with normalized data")
                 ),
                 p(strong("Note:"), "The app can process both CSV and RDS file formats. Make sure to select the correct organism and gene ID format.")
             )
      )
    ),
    
    # Sidebar layout with input and output definitions ----
    sidebarLayout(
      sidebarPanel(
        width = 4,
        div(class = "well",
            h4(style = "color: #3b82f6; margin-bottom: 20px;", "💾 Download Example Data"),
            downloadButton("Example", "Example DataSet", 
                           class = "btn btn-primary btn-block", 
                           style = "width:100%; margin-bottom: 25px;"),
            
            tags$hr(style = "border-top: 2px solid #e2e8f0; margin: 25px 0;"),
            
            h4(style = "color: #8b5cf6; margin-bottom: 20px;", "⚙️ Configuration"),
            
            h5("Organism Selection"),
            radioButtons(
              "organismus",
              label = NULL,
              choices = c("H. sapiens" = "Human",
                          "M. musculus" = "Mouse"),
              selected = "Mouse",
              inline = TRUE
            ),
            
            br(),
            
            h5("Gene ID Format"),
            radioButtons(
              "GENEid",
              label = NULL,
              choices = c(
                "Gene Symbol" = "SYMBOL",
                "Ensembl ID" = "EnsemblGenes",
                "Entrez ID" = "ENTREZID"
              ),
              selected = "SYMBOL",
              inline = TRUE
            ),
            
            tags$hr(style = "border-top: 2px solid #e2e8f0; margin: 25px 0;"),
            
            h4(style = "color: #10b981; margin-bottom: 20px;", "📁 File Upload"),
            
            h5("RDS File"),
            fileInput("rdsFile", label = NULL, 
                      buttonLabel = "Browse...",
                      placeholder = "No file selected",
                      multiple = FALSE,
                      accept = ".rds"),
            
            tags$hr(style = "border-top: 1px solid #e2e8f0; margin: 20px 0;"),
            
            h5("CSV File"),
            fileInput(
              "file1",
              label = NULL,
              buttonLabel = "Browse...",
              placeholder = "No file selected",
              multiple = FALSE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")
            ),
            
            tags$hr(style = "border-top: 2px solid #e2e8f0; margin: 25px 0;"),
            
            h4(style = "color: #3b82f6; margin-bottom: 20px;", "🔧 CSV Options"),
            
            h5("Header Row:"),
            radioButtons("header", label = NULL,
                         choices = list("Yes" = TRUE, "No" = FALSE),
                         selected = TRUE,
                         inline = TRUE),
            
            br(),
            
            h5("Separator:"),
            radioButtons(
              "sep",
              label = NULL,
              choices = c(
                "Comma (,)" = ",",
                "Semicolon (;)" = ";",
                "Tab" = "\t"
              ),
              selected = ",",
              inline = FALSE
            ),
            
            h5("Quote Character:"),
            radioButtons(
              "quote",
              label = NULL,
              choices = c(
                "None" = "",
                'Double Quote (")' = '"',
                "Single Quote (')" = "'"
              ),
              selected = '"',
              inline = FALSE
            ),
            
            tags$hr(style = "border-top: 2px solid #e2e8f0; margin: 25px 0;"),
            
            h5("Display Options:"),
            radioButtons(
              "disp",
              label = NULL,
              choices = c("Head (first 10 rows)" = "head",
                          "All rows" = "all"),
              selected = "head",
              inline = FALSE
            )
        )
      ),
      
      # Main panel for displaying outputs ----
      mainPanel(
        width = 8,
        div(class = "well",
            h4(style = "color: #1e293b; margin-bottom: 20px;", "Data Preview"),
            br(),
            div(class = "shiny-plot-output",
                tableOutput("contents"),
                br(),
                tableOutput("Rvalue")
            )
        )
      )
    )
  ),

    # Run Analysis Tab
    tabPanel(
      title = "Run Analysis",
      value = "analysis",
      # Top section: Data Filtering and Normalization (horizontal layout)
      fluidRow(
        column(12,
          div(class = "well", style = "padding: 20px; margin-bottom: 20px;",
            fluidRow(
              column(3,
                h5(style = "margin-top: 0; color: #1e293b;", "Data Filtering:"),
                radioButtons(
                  "VarFilter",
                  label = NULL,
                  choices = c(
                    "Remove Low Variance" = "Strict_Filter",
                    "No Filter" = "Unselect"
                  ),
                  selected = "Strict_Filter",
                  inline = TRUE
                ),
                numericInput("uniqueCut", "Low Variance cutoff", 15, 0, 100, width = "100%")
              ),
              column(3,
                h5(style = "margin-top: 0; color: #1e293b;", "Normalization:"),
                radioButtons(
                  "Norm",
                  label = NULL,
                  c("Normalization" = "Normal",
                    "No-Normalization" = "No_Normal"),
                  selected = "Normal",
                  inline = TRUE
                )
              ),
              column(3,
                br(),
                dropMenu(
                  dropdownButton(
                    "Info",
                    status = 'info',
                    size = "xs",
                    icon = icon('info-circle')
                  ),
                  h3(strong('Variance Filter Information')),
                  br(),
                  h5(
                    div(style = "color: #3b82f6; font-weight: 600;",
                        strong("Remove Low Variance:"), ),
                    "
Using the nearZeroVar function from R package:Caret identifies predictors with one unique value (zero variance predictors) or predictors with both of the following characteristics: they have a small number of unique values compared to the number of samples and a high frequency of the most frequent value. A threshold option is provided as a cutoff for the percentage of distinct values in relation to the total number of samples in the dataset.
" ,
                    br(),
                    div(style = "color: #3b82f6; font-weight: 600;",
                        strong("Normalization :"), ),
                    "
we utilize the Seurat package for data normalization, a crucial step in single-cell RNA sequencing analysis. Specifically, we determine a scaling factor based on the mean of the column sums of the dataset. This scaling factor is then used in the LogNormalize function from Seurat to perform log normalization on the data.
"
                  ),
                  placement = "bottom-start",
                  arrow = TRUE,
                  theme = "material",
                  maxWidth = 1000
                )
              ),
              column(3,
                br(),
                dropMenu(
                  dropdownButton(
                    "Info",
                    status = 'info',
                    size = "xs",
                    icon = icon('info-circle')
                  ),
                  h3(strong('Information')),
                  br(),
                  h5(
                    "In this suite of analysis tools, a wide range of strategies are provided for the analysis of the data. These include methods for identifying variable genes based on single cell RNA sequencing, statistical approaches based on p-value and LogFC threshold, and Feature Selection through tree-based ML models. A unique feature of this suite is the ability for the user to create their own ensemble method by combining one or more of these methods. Additionally, various threshold options are available for each method (for more information, consult the tutorial provided). It should be noted that each method employs a distinct approach and may require varying amounts of computational resources, depending on the size of the data. Therefore, it is recommended that the user employs preprocessing filters to reduce the data size and improve the reliability of the results. The normalization button allows the user to normalize their data prior to analysis. On the Ensemble tab, the options selected from the previous tabs for the various methods are applied. The genes identified through the analysis are depicted in a visual format utilizing both a barplot and a heatmap. The heatmap specifically uses cell type predict or State labeling as a method to organize and present the data. The state retains the labels provided in the dataset by the user, or the cell type prediction utilizes the singleR package for cell type annotation. Furthermore, a classification K-nearest neighbors (K-nn) model is executed to assess the ability of the model to accurately classify based solely on the isolated genes. The results of the K-NN analysis are presented in a comprehensive manner through the utilization of a confusion matrix and a data table."
                  ),
                  placement = "bottom-start",
                  arrow = TRUE,
                  theme = "material",
                  maxWidth = 1000
                )
              )
            )
          )
        )
      ),
      
      # Main layout with sidebar and main content
      sidebarLayout(
        # Sidebar Panel for Genes Selection Methods
        sidebarPanel(
          width = 3,
          strong("Genes Selection Methods", style = "color:black; font-size: 16px;"),
          br(),
          br(),
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
              
              fluidRow(
                column(4, div(style = "color:black", "SCMarker:")),
                column(4, numericInput("geneK", "geneK", value = 20, min = 0, step = 1)),
                column(4, numericInput("cellK", "cellK", value = 20, min = 0, step = 1))
              ),
              fluidRow(
                column(4, div(style = "color:black", "scran:")),
                column(8, numericInput("np", "num.Genes", value = 300, min = 1, step = 1))
              ),
              fluidRow(
                column(4, div(style = "color:black", "ScPNMF:")),
                column(8, numericInput("gM", "M Genes #", value = 300, min = 1, step = 1))
              ),
              fluidRow(
                column(4, div(style = "color:black", "M3Drop:")),
                column(4, numericInput("M3dropthreshold", "mt_threshold", value = 0.001, min = 0, max = 1, step = 0.01)),
                column(4, selectInput("M3Method", "mt_method", c("bon" = "bon", "fdr" = "fdr")))
              ),
              fluidRow(
                column(4, div(style = "color:black", "VST:")),
                column(4, numericInput("n", "# Features", value = 300, min = 1, step = 1)),
                column(4, selectInput("distMethod", "distMethod", c("EucDist" = "EucDist", "KL" = "KL", "DPNMF" = "DPNMF")))
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
              
              numericInput(
                "PvalueNum",
                "P-Value Threshold",
                value = 0.01,
                min = 0.001,
                max = 0.99,
                step = 0.01
              ),
              numericInput(
                "logfc",
                "LogFC Threshold",
                value = 1,
                min = 0,
                max = 5,
                step = 0.1
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
              numericInput(
                "importanceLimit",
                "Importance Threshold",
                value = 10,
                min = 0,
                max = 100,
                step = 1
              ),
              HTML('<div style="margin-top: 10px; padding: 10px; background-color: #f8f9fa; border-left: 4px solid #007bff; border-radius: 4px;">
                <small><strong>Note:</strong> Tree-based ML methods use feature importance scores to select the most important genes for classification.</small>
              </div>')
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
              numericInput(
                "SHAP_importanceLimit",
                "SHAP Importance Threshold",
                value = 0.01,
                min = 0,
                max = 100,
                step = 0.001
              ),
              HTML('<div style="margin-top: 10px; padding: 10px; background-color: #fff3cd; border-left: 4px solid #ffc107; border-radius: 4px;">
                <small><strong>Tip:</strong> For SHAP values, use small thresholds (0.001-0.1) for absolute SHAP, or integers (1-100) for "top N features". SHAP values are in their original scale (mean absolute SHAP).</small>
              </div>')
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
          width = 9,
          
          # Action Buttons Section
          fluidRow(
            column(12, align = "center",
              div(style = "margin: 30px 0;",
                actionButton("click", 
                            "Run Analysis",
                            class = "btn btn-primary btn-lg",
                            style = "margin-right: 20px; padding: 15px 40px; font-size: 18px;"),
                downloadLink(
                  "downloadData",
                  "Download Results",
                  class = "btn btn-success btn-lg",
                  style = "padding: 15px 40px; font-size: 18px; text-decoration: none; display: inline-block;"
                )
              )
            )
          ),
          
          # Analysis Controls
          fluidRow(
            column(12,
              div(class = "well", style = "padding: 20px; margin-bottom: 20px;",
                fluidRow(
                  column(4,
                    sliderInput(
                      "genes",
                      "Number of genes:",
                      min = 2,
                      max = 500,
                      value = 100
                    )
                  ),
                  column(4,
                    div(class = "alert alert-warning", style = "margin-top: 25px; padding: 10px;",
                      HTML('<strong>Important:</strong> Don\'t forget to choose the right Organism and Gene ID!')
                    )
                  ),
                  column(4,
                    h5("Generate Heatmap:", style = "margin-top: 10px; margin-bottom: 5px;"),
                    radioButtons("HeatMap1", label = NULL,
                                choices = list("Yes" = TRUE, "No" = FALSE),
                                selected = FALSE,
                                inline = TRUE)
                  ),
                  column(4,
                    h5("Heatmap Clustering:", style = "margin-top: 10px; margin-bottom: 5px;"),
                    radioButtons(
                      "clustering",
                      label = NULL,
                      choices = c("Between Groups" = "cluster_between_groups",
                                  "None" = "None_Clustering"),
                      selected = "None_Clustering",
                      inline = TRUE
                    )
                  ),
                  column(4,
                    h5("Heatmap Splitting:", style = "margin-top: 10px; margin-bottom: 5px;"),
                    radioButtons(
                      "Split",
                      label = NULL,
                      choices = c("Cell Type Predict" = "CellType",
                                  "State" = "labels"),
                      selected = "CellType",
                      inline = TRUE
                    )
                  )
                )
              )
            )
          ),
          
          # Error message
          fluidRow(
            column(12,
              tableOutput("error")
            )
          ),
          
          # Status text
          fluidRow(
            column(12,
              tableOutput("text")
            )
          ),
          
          # Report Section - Bar Plot and Genes List
          fluidRow(
            column(8,
              div(class = "shiny-plot-output",
                h4(style = "color: #1e293b; margin-bottom: 15px;", 
                   "Most Important Genes which Act as Potential Biomarkers"),
                plotOutput("TheBarPlot", height = "500px")
              )
            ),
            column(4,
              div(class = "well", style = "height: 500px; overflow-y: auto; padding: 15px;",
                h4(style = "color: #1e293b; margin-bottom: 15px;", 
                   "# OF GENES WHICH OPERATE AS BIOMARKERS"),
                tableOutput("GenesList")
              )
            )
          ),
          
          br(),
          
          # K-NN Classifier
          fluidRow(
            column(12,
              div(class = "shiny-plot-output",
                h4(style = "color: #1e293b; margin-bottom: 15px;", 
                   "k-NN Classification Results"),
                plotOutput("KnnClassifier", height = "500px")
              )
            )
          ),
          
          br(),
          
          # Heatmap Section
          fluidRow(
            column(8,
              div(class = "shiny-plot-output",
                h4(style = "color: #1e293b; margin-bottom: 15px;", 
                   "Heatmap Visualization"),
                plotOutput("HeatMap", height = "600px")
              )
            ),
            column(4,
              div(class = "well", style = "height: 600px; overflow-y: auto; padding: 15px;",
                h4(style = "color: #1e293b; margin-bottom: 15px;", 
                   "Heatmap Gene List"),
                tableOutput("HeatmapList")
              )
            )
          ),
          
          br()
        )
      )
    ),
    
    # Enrichment Analysis Tab
    tabPanel(
      title = "Enrichment Analysis",
      value = "enrichment",
      
      #Info Button
      dropMenu(
        dropdownButton(
          "Info",
          status = 'info',
          size = "xs",
          icon = icon('info-circle'),
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
            h5("All Available Ontology Terms:"),
            radioButtons("all", label = NULL,
                         choices = list("Yes" = TRUE, "No" = FALSE),
                         selected = FALSE,
                         inline = TRUE)
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
                 actionButton("click1", "Run Enrichment Analysis",
                             class = "btn btn-primary btn-block")),
           column(
             7,
          # 
          #   downloadButton("GP", "Export Plots"),
          #   downloadButton("GT", "Export Tables"),
          ),
          br(),
          fluidRow(
            column(12,
              div(class = "shiny-plot-output",
                plotOutput("BioBarPlot", width = "100%", height = "1300px")
              ),
              br(),
              tableOutput('Enrichment')
            )
          )
        )
      )
    ),
    
    # KEGG Maps Tab
    tabPanel(
      title = "KEGG Maps",
      value = "kegg",
      dropMenu(
        dropdownButton(
          "Info",
          status = 'info',
          size = "xs",
          icon = icon('info-circle')
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
            12,
            align = "center",
            actionButton("click2", "Visualize Pathway",
                        class = "btn btn-primary btn-lg",
                        style = "margin-top: 20px;")
          ),
        ),
        
        # Main panel for displaying outputs ----
        mainPanel(
          dataTableOutput("KEGG"),
          textOutput ("text1"),
          
          sidebarPanel(
            width = 12,
            
            
            plotOutput("KEGGmap", width = "100%",height = 1000) ,#, width = "100%", height = 500
            
            
          )
        )
        
        
      )
    ),
    
    # Graph Analysis Tab
    tabPanel(
      title = "Graph Analysis",
      value = "graphs",
      
      dropMenu(
        dropdownButton(
          "Info",
          status = 'info',
          size = "xs",
          icon = icon('info-circle')
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
              value = 50
            ) ) ),
         
        h5("Protein–protein interaction (PPI) network analysis:"),
        radioButtons("PPInetwork1", label = NULL,
                    choices = list("Yes" = TRUE, "No" = FALSE),
                    selected = FALSE,
                    inline = TRUE),
        
        numericInput(
          "Score_Threshold_PPI",
          "Conditionally load interactions based on a threshold.",
          value = 400,
          min = 50,
          max = 1000,
          step = 10
        ),
        
        h5("Similarity Graph:"),
        radioButtons("graph1", label = NULL,
                    choices = list("Yes" = TRUE, "No" = FALSE),
                    selected = FALSE,
                    inline = TRUE),
        numericInput(
          "Pearson_correlation",
          " Remove edges below absolute Pearson correlation",
          value = 0.5,
          min = 0.1,
          max = 0.99,
          step = 0.1
        ),
        
        
        column(
          width = 12,
          align = "center",
          actionButton("run_button1", "Generate Graphs",
                      class = "btn btn-primary btn-lg",
                      style = "margin-top: 20px;")
        ),
      ),
      # Main panel for displaying outputs ----
      mainPanel(
        width = 8,
        tableOutput("text3"),
        
        div(class = "shiny-plot-output",
          plotOutput("PPInetwork", width = "100%", height = "850px"),
          br(),
          visNetworkOutput("graph", width = "100%", height = "850px")
        )
      )
    )
  )

