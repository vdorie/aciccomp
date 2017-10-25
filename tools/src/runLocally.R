## dirs should be a named list with members data, methods, results, and log
runLocally <- function(runStatus, methods, dirs, runMethods = NULL)
{
  if (length(runStatus) == 0L) return(invisible(NULL))
  
  runCaseNames <- names(runStatus[[1L]])
  
  x.comp <- read.csv(file.path(dirs$data, "x.csv"))
  
  df <- data.frame(z = integer(nrow(x.comp)), y = numeric(nrow(x.comp)))
  for (i in seq_along(x.comp)) df[[names(x.comp)[i]]] <- x.comp[[i]]
  rm(x.comp)
  
  inFile <- tempfile(fileext = ".csv")
  
  for (i in seq_len(nrow(methods))) {
    method <- methods[i,]
    if (!is.null(runMethods) && length(runMethods) > 0L && !(method$name %in% runMethods)) next
    
    cat("fitting method '", method$name, "'\n", sep = "")
    
    command <- file.path(dirs$methods, method$name)
    
    for (j in seq_along(runCaseNames)) {
      iters <- runStatus[[i]][[j]]$iter
      runCaseDir <- gsub("\\/", .Platform$file.sep, runCaseNames[j])
      
      logDir <- file.path(dirs$log, method$name)
      if (!dir.exists(logDir)) dir.create(logDir, recursive = TRUE)
      
      logFile <- file.path(logDir, paste0(gsub("\\/", "_", runCaseNames[j]), ".o"))
      errorFile <- file.path(logDir, paste0(gsub("\\/", "_", runCaseNames[j]), ".e"))
      logCon <- file(logFile, "wt")
      errorCon <- file(errorFile, "wt")
      
      resultsDir <- file.path(dirs$results, method$name, runCaseDir)
      if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)
    
      currentResultIndices <- which(runStatus[[i]][[j]]$status %in% c("complete", "hung"))
      if (length(currentResultIndices) == length(iters)) next
    
      cat("fitting method '", method$name, "' in setting '", runCaseNames[j], "'\n", sep = "", file = logCon)
      
      dataDir <- file.path(dirs$data, runCaseDir)
      
      for (k in seq_along(iters)) {
        if (k %in% currentResultIndices) next
        
        dataFile <- file.path(dataDir, paste0(iters[k], ".csv"))
        
        respHasHeaders <- grepl("['\"]z['\"]\\s*,\\s*['\"]y['\"]", readLines(dataFile, n = 1L), perl = TRUE)
        resp <- if (respHasHeaders) read.csv(dataFile) else read.csv(dataFile, header = FALSE, col.names = c("z", "y"))
        df$z <- resp$z
        df$y <- resp$y
        #df$y <- ifelse(resp$z == 0L, resp$y_0, resp$y_1)
        write.table(df, inFile, sep = ",", dec = ".",
                    row.names = FALSE, col.names = method$headers_in == 1L)
        
        outFile <- file.path(resultsDir, paste0(iters[k], ".csv"))
        
        args <- c(inFile, outFile)
    
        if (method$individual_effects != 0L) {
          outFile.ind <- file.path(resultsDir, paste0(iters[k], "_ind.csv"))
          args <- c(args, outFile.ind)
        }
        
        cat(runCaseNames[j], ",", iters[k], ",", sep = "", file = logCon)
        startTime <- proc.time()
        consoleLog <- system2(command, args, stdout = TRUE, stderr = TRUE)
        timeDiff <- proc.time() - startTime
        
        if (!is.null(attr(consoleLog, "status"))) {
          cat("error in: '", runCaseNames[j], "' iter ", iters[k], "\n", sep = "", file = errorCon)
          cat(consoleLog, sep = "\n", file = errorCon)
          cat("\n\n", file = errorCon)
          
          cat("NA\n", file = logCon)
          result <- data.frame(est = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
          write.table(result, file = outFile, sep = ",", dec = ".", row.names = FALSE, method$headers_out == 1L)
          
          if (method$individual_effects != 0L) {
            result <- data.frame(
              est =  rep_len(NA_real_, nrow(df)),
              ci_lower = rep_len(NA_real_, nrow(df)),
              ci_upper = rep_len(NA_real_, nrow(df)))
            
            write.table(result, file = outFile.ind, sep = ",", dec = ".", row.names = FALSE, method$headers_out == 1L)
          }
        } else {
          cat(timeDiff[["sys.child"]], "\n", sep = "", file = logCon)
        }
      }
      close(logCon)
      close(errorCon)
    }
  }
  unlink(inFile)
}
