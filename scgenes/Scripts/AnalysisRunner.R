# One cached computation per Run Analysis click, isolated to the current session.
# Catch ordinary errors here: an uncaught observer error closes a Shiny session.
create_analysis_runner <- function(input, run_method) {
  shiny::eventReactive(input$click, {
    shiny::removeNotification("analysis-error")
    tryCatch({
      result <- run_method()
      if (!is.list(result) || is.null(result$ig) || is.null(result$newdata)) {
        stop("No gene-selection result was returned. Check the selected method and input data.")
      }
      result
    }, error = function(error) {
      if (inherits(error, "shiny.silent.error")) {
        message_text <- conditionMessage(error)
        if (!nzchar(message_text)) return(NULL)
      } else {
        message_text <- conditionMessage(error)
        message("Gene selection failed: ", message_text)
      }
      if (grepl("lazy-load database|R_decompress|namespace.*failed", message_text, ignore.case = TRUE)) {
        message_text <- paste(
          "A required R package could not be loaded. Restart the R session and run the app again.",
          "If the error persists, reinstall the package named in the R console."
        )
      }
      shiny::showNotification(message_text, type = "error", duration = NULL,
                              id = "analysis-error")
      NULL
    })
  }, ignoreInit = TRUE)
}

# Run/Stop control for the analysis button pair.
#
# R is single threaded: while `analysis_result()` computes, the session cannot
# read anything from the browser. The button swap is therefore done client side
# (www/analysis-run-control.js) so a run can never be started twice, and a click
# on Stop is only processed once the run has released the process - it then
# discards that run and puts the UI back to its pre-run state.
register_run_controls <- function(input, output, session, analysis_result) {
  state <- shiny::reactiveValues(running = FALSE, stopRequested = FALSE, ready = FALSE)

  notify <- function() {
    session$sendCustomMessage("scgenes-analysis-state", list(
      running = isTRUE(state$running),
      stopping = isTRUE(state$running) && isTRUE(state$stopRequested)
    ))
  }

  # Drives every "results are available" panel. Unlike the click counter this
  # can go back to FALSE, which is what makes a clean re-run possible.
  # A plain function (not renderText) so the browser receives a JSON boolean;
  # the `...` keeps it callable both by a live session and by testServer.
  output$analysisReady <- function(...) {
    isTRUE(state$ready)
  }
  shiny::outputOptions(output, "analysisReady", suspendWhenHidden = FALSE)

  # The priority keeps this ahead of the result outputs, so the state is
  # settled before anything renders.
  #
  # No ignoreInit here on purpose: Shiny reports the button as 0 when the
  # browser connects, and that pass has to reach the event reactive to prime
  # its own ignoreInit. Without it the first real click would be swallowed.
  shiny::observeEvent(input$click, {
    priming <- !isTRUE(input$click > 0)
    if (!priming) {
      shiny::removeNotification("analysis-stopped")
      state$running <- TRUE
      state$stopRequested <- FALSE
      state$ready <- FALSE
    }
    tryCatch({
      result <- analysis_result()
      if (!priming) {
        state$ready <- is.list(result) && !is.null(result$ig) && !is.null(result$newdata)
      }
    }, finally = {
      if (!priming) {
        state$running <- FALSE
        notify()
      }
    })
  }, priority = 100)

  shiny::observeEvent(input$stopAnalysis, {
    state$running <- FALSE
    state$stopRequested <- TRUE
    state$ready <- FALSE
    notify()
    shiny::showNotification(
      paste("Analysis stopped. The results of that run were discarded -",
            "press Run analysis to start again."),
      type = "warning", duration = 8, id = "analysis-stopped"
    )
  }, ignoreInit = TRUE)

  invisible(state)
}
