library(shiny)
source("Scripts/AnalysisRunner.R")

testServer(function(input, output, session) {
  state <- new.env()
  state$calls <- 0L
  state$fail <- FALSE
  result <- create_analysis_runner(input, function() {
    state$calls <- state$calls + 1L
    if (state$fail) stop("gene selection blew up")
    list(ig = data.frame(score = 1, row.names = "GeneA"),
         newdata = data.frame(GeneA = 1:4))
  })
  runState <- register_run_controls(input, output, session, result)
}, {
  # Both action buttons report 0 when the browser connects.
  session$setInputs(click = 0, stopAnalysis = 0)
  stopifnot(
    "connecting must not start an analysis" = state$calls == 0L,
    "no results before the first run" = isFALSE(output$analysisReady),
    "the run button starts released" = isFALSE(runState$running)
  )

  session$setInputs(click = 1)
  stopifnot(
    "one click runs the analysis once" = state$calls == 1L,
    "results become visible after a run" = isTRUE(output$analysisReady),
    "the button is released once the run ends" = isFALSE(runState$running)
  )

  # Stop is only read after the run released the process; it discards the run.
  session$setInputs(stopAnalysis = 1)
  stopifnot(
    "a stop hides the results again" = isFALSE(output$analysisReady),
    "a stop releases the button" = isFALSE(runState$running),
    "a stop is recorded" = isTRUE(runState$stopRequested),
    "stopping does not trigger another run" = state$calls == 1L
  )

  # A stopped analysis can be started again.
  session$setInputs(click = 2)
  stopifnot(
    "a stopped analysis can be re-run" = state$calls == 2L,
    "the re-run shows its results" = isTRUE(output$analysisReady),
    "the stop flag is cleared on re-run" = isFALSE(runState$stopRequested),
    "the re-run releases the button" = isFALSE(runState$running),
    "the session survives a re-run" = !session$isClosed()
  )

  # A failing run must release the button instead of leaving it stuck on Stop.
  state$fail <- TRUE
  session$setInputs(click = 3)
  stopifnot(
    "a failing run is still attempted once" = state$calls == 3L,
    "a failed run shows no results" = isFALSE(output$analysisReady),
    "a failed run still releases the button" = isFALSE(runState$running),
    "a failed run does not close the session" = !session$isClosed()
  )

  # Unrelated inputs never touch the run state.
  session$setInputs(genes = 20)
  stopifnot("unrelated inputs do not re-run the analysis" = state$calls == 3L)
})
cat("Run/Stop: one run per click, stop clears the results, re-run and failures reset the button.\n")
