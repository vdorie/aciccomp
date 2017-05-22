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

updateResults <- function(runStatus, results, methods, dirs, functions, runMethods = NULL)
{
  if (length(runStatus) == 0L) return(results)
  
  for (i in seq_along(runStatus)) {
    methodName <- names(runStatus)[i]
    if (!is.null(runMethods) && length(runMethods) > 0L && !(methodName %in% runMethods)) next
    method <- methods[which.max(methods$name %in% methodName),]
    
    cat("updating results for method '", methodName, "'\n", sep = "")
    
    for (j in seq_along(runStatus[[i]])) {
      completedIndices <- runStatus[[i]][[j]]$status == "complete"
      unevaluatedIndices <- rep.int(FALSE, length(completedIndices))
      for (k in seq_along(results)) {
        if (names(results)[k] == "indiv" && method$individual_effects == 0L) next
        unevaluatedIndices <- unevaluatedIndices | is.na(results[[k]][i,j,])
      }
      unevaluatedIndices <- unevaluatedIndices & completedIndices
      
      if (!any(unevaluatedIndices)) next
      
      runCaseName <- gsub("/", .Platform$file.sep, names(runStatus[[i]])[j])
      
      dgp <- read.csv(file.path(dirs$data, runCaseName, "dgp.csv"))
      
      for (iter in which(unevaluatedIndices)) {
        dataFile        <- file.path(dirs$data, runCaseName, paste0(iter, ".csv"))
        resultsFile     <- file.path(dirs$results, methodName, runCaseName, paste0(iter, ".csv"))
        resultsFile.ind <- file.path(dirs$results, methodName, runCaseName, paste0(iter, "_ind.csv"))
        
        respHasHeaders <- grepl("['\"]z['\"]\\s*,\\s*['\"]y['\"]", readLines(dataFile, n = 1L), perl = TRUE)
        resp <- if (respHasHeaders) read.csv(dataFile) else read.csv(dataFile, header = FALSE, col.names = c("z", "y"))
        results.ij <- read.csv(resultsFile, header = method$headers_out == 1L)
        results.ij.ind <- if (file.exists(resultsFile.ind)) read.csv(resultsFile.ind, header = method$headers_out == 1L) else NULL
        
        if (method$headers_out != 1L) {
          colnames(results.ij) <- c("est", "ci_lower", "ci_upper")
          if (!is.null(results.ij.ind)) colnames(results.ij.ind) <- colnames(results.ij)
        }
        
        for (k in seq_along(functions)) {
          results[[names(functions)[k]]][i,j,iter] <- functions[[k]](resp, dgp, results.ij, results.ij.ind)
        }
      }
    }
  }
  results
}
