# File uploads are delivered to the server, so browser conditional panels must
# use a server output rather than input.file1/input.rdsFile.
register_upload_preview <- function(input, output) {
  output$uploadState <- shiny::renderText({
    if (!is.null(input$file1) || !is.null(input$rdsFile)) "ready" else "empty"
  })
  shiny::outputOptions(output, "uploadState", suspendWhenHidden = FALSE)
  #Read the import .rds dataset and visualize it.
  output$Rvalue <- DT::renderDataTable({
    shiny::req(input$rdsFile)

    RDS_file <- input$rdsFile
    if (is.null(RDS_file)) {
      return()
    } else {
      rd <- read_uploaded_rds(RDS_file)
      DT::datatable(preview_expression_data(rd),
        options = list(pageLength = 10, scrollX = TRUE, dom = 'frtip'),
        class = 'cell-border stripe')
    }

  })
  
  #Read the import .csv dataset and visualize it.
  output$contents <- DT::renderDataTable({
    shiny::req(input$file1)
    header_val <- if(is.null(input$header)) TRUE else as.logical(input$header)
    sep_val <- if(is.null(input$sep)) "," else input$sep
    quote_val <- if(is.null(input$quote)) "\"" else input$quote

    CSV_file <- read_uploaded_csv(
      input$file1,
      header = header_val,
      sep = sep_val,
      quote = quote_val
    )
    data_to_show <- if (is.null(input$disp) || input$disp == "head") {
      preview_expression_data(CSV_file)
    } else {
      preview_expression_data(CSV_file, max_rows = 1000L, max_columns = 100L)
    }
    DT::datatable(data_to_show,
      options = list(pageLength = 10, scrollX = TRUE, dom = 'frtip'),
      class = 'cell-border stripe')
  })
  
 

}
