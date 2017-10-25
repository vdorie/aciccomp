## dirs should be a named list with at least 'data'
collectRunCasesFromDirectory <- function(dirs) {
  runCaseNames <- character()
  runCaseIterations <- integer()
  
  topLevelDirs <- list.dirs(dirs$data, recursive = FALSE, full.names = FALSE)
  for (dir in topLevelDirs) {
    files <- strsplit(list.files(file.path(dirs$data, dir), "^[0-9]+\\.csv", recursive = TRUE), .Platform$file.sep)
    
    ## force '/' as separator
    names <- sapply(files, function(file) if (length(file) > 1L) paste0(dir, "/", paste0(head(file, -1L), collapse = "/")) else dir)
    iterations <- as.integer(sub("^([0-9]+)\\.csv$", "\\1", sapply(files, function(file) tail(file, 1L))))
    runCaseNames      <- c(runCaseNames, names)
    runCaseIterations <- c(runCaseIterations, iterations)
  }
  
  runCases <- data.frame(name = runCaseNames, iter = runCaseIterations)
  ## enforce an official order that is hopefully platform independent
  runCaseNames <- sort(unique(runCases$name))
  runCaseIndices <- lapply(runCaseNames, function(x) which(runCases$name == x))
  
  data.frame(name = rep(runCaseNames, sapply(runCaseIndices, length)),
             iter = unlist(lapply(runCaseIndices, function(x) sort(runCases$iter[x]))))
}
