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

updateResults <- function(runStatus, results, dirs, functions, runMethods = NULL)
{
  if (length(runStatus) == 0L) return(results)
  
  for (i in seq_along(runStatus)) {
    methodName <- names(runStatus)[i]
    if (!is.null(runMethods) && !(methodName %in% runMethods)) next
    
    cat("updating results for method '", methodName, "'\n")
    
    for (j in seq_along(runStatus[[i]])) {
      completedIndices <- runStatus[[i]][[j]]$status == "complete"
      unevaluatedIndices <- is.na(results[[1L]][i,j,]) & completedIndices
      
      if (!any(unevaluatedIndices)) next
      
      runCaseName <- names(runStatus[[i]])[j]
      
      for (iter in which(unevaluatedIndices)) {
        dataFile        <- file.path(dirs$data, runCaseName, paste0(iter, ".csv"))
        resultsFile     <- file.path(dirs$results, methodName, runCaseName, paste0(iter, ".csv"))
        resultsFile.ind <- file.path(dirs$results, methodName, runCaseName, paste0(iter, "_ind.csv"))
        
        data <- read.csv(dataFile)
        results.ij <- read.csv(resultsFile)
        results.ij.ind <- if (file.exists(resultsFile.ind)) read.csv(resultsFile.ind) else NULL
        
        for (k in seq_along(functions)) {
          results[[names(functions)[k]]][i,j,iter] <- functions[[k]](data, results.ij, results.ij.ind)
        }
      }
    }
  }
  results
}

if (FALSE) {
  source("site_setup.R")
  load("runStatus.Rdata")
  
  results <- createResults(runStatus, c("bias", "coverage"))
  save(results, file = "results.Rdata")
}

if (FALSE) {
  source("site_setup.R")
  load("runStatus.Rdata")
  load("results.Rdata")
  
  ## not actually bias, obviously
  biasFn <- function(data, results, results.ind) {
    (mean(data$y[data$z == 1]) - mean(data$y[data$z == 0])) - results$est
  }
  ## and obviously not coverage
  coverageFn <- function(data, results, results.ind) {
    te <- mean(data$y[data$z == 1]) - mean(data$y[data$z == 0])
    results$ci_lower <= te & te <= results$ci_upper
  }
  
  results <- updateResults(runStatus, results, dirs, list(bias = biasFn, coverage = coverageFn))
}
