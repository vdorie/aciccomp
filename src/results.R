createResults <- function(runStatus, valueNames)
{
  results <- vector("list", length(valueNames))
  names(results) <- valueNames
 
  ## numMethods * numSettings * numIters
  ## in case numIters isn't shared across all settings, it takes the max
  
  numIters <- max(sapply(runStatus[[1L]], nrow))
  for (i in seq_along(results)) {
    results[[i]] <- array(NA_real_, c(length(runStatus), length(runStatus[[1L]]), numIters),
                          dimnames = list(names(runStatus), names(runStatus[[1L]]), NULL))
  }
  results
}

addMethodsToResults <- function(results, methods) {
  results.old <- results
  numMethods.old <- dim(results[[1L]])[1L]
  
  dims.new <- dim(results[[1L]])
  dims.new[1L] <- dims.new[1L] + length(methods)
  methodNames <- c(dimnames(results[[1L]])[[1L]], methods)
  settingNames <- dimnames(results[[1L]])[[2L]]
  
  for (i in seq_along(results)) {
    results[[i]] <- array(NA_real_, dims.new, dimnames = list(methodNames, settingNames, NULL))
    results[[i]][seq_len(numMethods.old),,] <- results.old[[i]]
  }
  results
}

removeMethodsFromResults <- function(results, methods) {
  keepIndices <- !(dimnames(results[[1L]])[[1L]] %in% methods)
  
  for (i in seq_along(results))
    results[[i]] <- results[[i]][keepIndices,,drop = FALSE]
  
  results
}

updateResults <- function(runStatus, results, dirs, functions, runMethods = NULL)
{
  if (length(runStatus) == 0L) return(results)
  
  for (i in seq_along(runStatus)) {
    methodName <- names(runStatus)[i]
    if (!is.null(runMethods) && length(runMethods) > 0L && !(methodName %in% runMethods)) next
    
    cat("updating results for method '", methodName, "'\n", sep = "")
    
    for (j in seq_along(runStatus[[i]])) {
      completedIndices <- runStatus[[i]][[j]]$status == "complete"
      unevaluatedIndices <- is.na(results[[1L]][i,j,]) & completedIndices
      
      if (!any(unevaluatedIndices)) next
      
      runCaseName <- names(runStatus[[i]])[j]
      
      dgp <- read.csv(file.path(dirs$data, runCaseName, "dgp.csv"))
      
      for (iter in which(unevaluatedIndices)) {
        dataFile        <- file.path(dirs$data, runCaseName, paste0(iter, ".csv"))
        resultsFile     <- file.path(dirs$results, methodName, runCaseName, paste0(iter, ".csv"))
        resultsFile.ind <- file.path(dirs$results, methodName, runCaseName, paste0(iter, "_ind.csv"))
        
        data <- read.csv(dataFile)
        results.ij <- read.csv(resultsFile)
        results.ij.ind <- if (file.exists(resultsFile.ind)) read.csv(resultsFile.ind) else NULL
        
        for (k in seq_along(functions)) {
          results[[names(functions)[k]]][i,j,iter] <- functions[[k]](data, dgp, results.ij, results.ij.ind)
        }
      }
    }
  }
  results
}
