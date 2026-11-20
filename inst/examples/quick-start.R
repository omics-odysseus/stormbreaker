# Run the packaged binary and survival demos without sharing their variables.
quick_start <- lapply(c(binary = "07-route-workflow.R", survival = "04-survival.R"),
  function(script) {
    example <- new.env(parent = globalenv())
    sys.source(system.file("examples", "scenarios", script,
      package = "stormbreaker", mustWork = TRUE), envir = example)
    example$example_result
  })
