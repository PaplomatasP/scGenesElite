library(shiny)
source("Scripts/AnalysisRunner.R")

testServer(function(input, output, session) {
  state <- new.env()
  state$calls <- 0L
  state$fail <- TRUE
  result <- create_analysis_runner(input, function() {
    state$calls <- state$calls + 1L
    if (state$fail) stop("lazy-load database 'rtracklayer.rdb' is corrupt")
    list(ig = data.frame(score = 1, row.names = "GeneA"),
         newdata = data.frame(GeneA = 1:4))
  })
  observeEvent(input$click, result())
  output$first <- renderText(if (is.null(result())) "failed" else "ready")
  output$second <- renderText(if (is.null(result())) "failed" else "ready")
}, {
  session$flushReact()
  stopifnot(state$calls == 0L)
  session$setInputs(click = 1)
  stopifnot(identical(output$first, "failed"), identical(output$second, "failed"))
  stopifnot(state$calls == 1L, !session$isClosed())
  # Re-reading a failed run must not retry it from another output.
  stopifnot(is.null(result()), state$calls == 1L)
  state$fail <- FALSE
  session$setInputs(click = 2)
  stopifnot(identical(output$first, "ready"), identical(output$second, "ready"))
  stopifnot(state$calls == 2L, !session$isClosed())
  session$setInputs(genes = 20)
  stopifnot(state$calls == 2L)
})
cat("Analysis errors are contained; each click runs once and a later click recovers.\n")
