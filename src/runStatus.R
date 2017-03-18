createRunStatus <- function(runCases, methods)
{
  runStatus <- vector("list", nrow(methods))
  names(runStatus) <- methods$name
  
  runCaseNames <- sort(unique(runCases$name))
  statuses <- c("NA", "missing", "broken", "failed", "hung", "complete")
  
  for (i in seq_len(nrow(methods))) {
    runStatus[[i]] <- vector("list", length(runCaseNames))
    names(runStatus[[i]]) <- runCaseNames
    
    for (j in seq_along(runCaseNames)) {
      iters <- with(runCases, iter[name == runCaseNames[j]])
      
      runStatus[[i]][[j]] <- data.frame(iter = iters,
                                        runTime = rep_len(NA_real_, length(iters)),
                                        status  = factor(rep_len(1L, length(iters)), levels = seq_along(statuses), labels = statuses))
    }
  }
  runStatus
}

addMethodsToRunStatus <- function(runStatus, methods)
{
  emptyStatus <- runStatus[[1L]]
  for (i in seq_along(emptyStatus)) {
    emptyStatus[[i]][,"runTime"] <- NA_real_
    emptyStatus[[i]][,"status"] <- "NA"
  }
  
  for (i in seq_along(methods)) {
    runStatus[[length(runStatus) + 1L]] <- emptyStatus
    names(runStatus)[length(runStatus)] <- methods[i]
  }
  
  runStatus
}

removeMethodsFromRunStatus <- function(runStatus, methods)
{
  for (method in methods) runStatus[[method]] <- NULL
  runStatus
}

## dirs should be a named list with at least results and log
updateRunStatus <- function(runStatus, dirs, runMethods = NULL)
{
  if (length(runStatus) == 0L) return(invisible(NULL))
  
  methodNames <- names(runStatus)
  runCaseNames <- names(runStatus[[1L]])
  
  
  currentJobs <- data.frame(id = integer(0L), status = character(0L), name = character(0L), stringsAsFactors = FALSE)
  
  qstatPath <- suppressWarnings(system2("command", c("-v", "qstat"), stdout = TRUE))
  if (length(qstatPath) > 0L) {
    jobsFile <- tempfile()
    system2("qstat", stdout = jobsFile)
    if (file.info(jobsFile)$size > 0L) {
      currentJobs <- system2("grep", c("-Ee", "'^\\s+[0-9]+\\s+.*'"),
                             stdout = TRUE,
                             stdin  = jobsFile)
      if (length(currentJobs) > 0L) {
        currentJobs <- sapply(currentJobs, function(job) strsplit(trimws(job), "\\s+"))
        currentJobs <- data.frame(id     = unname(sapply(currentJobs, function(job) job[1L])),
                                  status = unname(sapply(currentJobs, function(job) job[5L])), stringsAsFactors = FALSE)
        currentJobs$name <- if (nrow(currentJobs) > 0L) sapply(seq_len(nrow(currentJobs)), function(i) {
          system2("sed", c("-nEe", "'s/^job_name:\\s+(\\S+)$/\\1/p'"), stdout = TRUE,
                  input = system2("qstat", c("-j", currentJobs$id[i]), stdout = TRUE))
        }) else character(0L)
      }
    }
    unlink(jobsFile)
  }
  
  for (i in seq_along(methodNames)) {
    methodName <- methodNames[i]
    if (!is.null(runMethods) && length(runStatus) > 0L && !(methodName %in% runMethods)) next
  
    cat("updating run status for method '", methodName, "'\n", sep = "")
    
    for (j in seq_along(runCaseNames)) {
      resultsDir <- file.path(dirs$results, methodName, gsub("\\/", .Platform$file.sep, runCaseNames[j]))
      iters <- runStatus[[i]][[j]]$iter
      
      jobName <- paste0(methodName, "_", gsub("/", "_", runCaseNames[j]))
      if (jobName %in% currentJobs$name) next
      
      if (dir.exists(resultsDir)) {
	for (k in seq_along(iters)) {
	  if (runStatus[[i]][[j]]$status[k] %in% c("complete", "hung")) next
	  
	  resultsFile <- file.path(resultsDir, paste0(runStatus[[i]][[j]]$iter[k], ".csv"))
	  if (file.exists(resultsFile)) {
	    result <- read.csv(resultsFile)
	    if (is.null(result$est) || length(result$est) == 0L) {
	      runStatus[[i]][[j]]$status[k] <- "broken"
	    } else if (is.na(result$est)) {
	      runStatus[[i]][[j]]$status[k] <- "failed"
	    } else {
	      runStatus[[i]][[j]]$status[k] <- "complete"
	    }
	  } else {
	    runStatus[[i]][[j]]$status[k] <- "missing"
	  }
	}
      }
      
      logFileName <- paste0(gsub("/", "_", runCaseNames[j]), ".o")
      logFile <- file.path(dirs$log, methodName, logFileName)
      if (file.exists(logFile)) {
	suppressWarnings(rawRunTimes <- system2("grep", c("-E", paste0("'^", runCaseNames[j], ",[0-9]+,([0-9.]+|(NA))$'")),
			                        stdin = logFile, stdout = TRUE))
	if (length(rawRunTimes) == 0L || !is.null(attr(rawRunTimes, "status"))) {
	  cat("  warning: could not parse ", paste0(methodName, .Platform$file.sep, logFileName), "\n", sep = "")
	} else {
	  stringConnection <- textConnection(rawRunTimes)
	  rawRunTimes.j <- read.csv(stringConnection, header = FALSE, col.names = c("setting", "iter", "runTime"))
	  close(stringConnection)
	  
	  for (k in seq_len(nrow(rawRunTimes.j))) {
	    if (is.na(rawRunTimes.j$runTime[k])) next
	    
	    index <- which(runStatus[[i]][[j]]$iter == rawRunTimes.j$iter[k])
	    if (runStatus[[i]][[j]]$status[index] == "complete")
	      runStatus[[i]][[j]]$runTime[index] <- rawRunTimes.j$runTime[k]
	  }
	}
	
	if (runStatus[[i]][[j]]$status[length(iters)] == "missing") {
	  lastLine <- system2("tail", c("-n", "1", logFile), stdout = TRUE)
	  if (grepl(paste0("^", runCaseNames[j], ",[0-9]+,$"), lastLine)) {
	    lastIter <- as.integer(strsplit(lastLine, ",")[[1]][2])
	    lastIndex <- which(runStatus[[i]][[j]]$iter == lastIter)
	    if (all(runStatus[[i]][[j]]$status[seq.int(lastIndex, numItersPerSetting)] == "missing")) {
	      cat("setting ", names(runStatus)[i], ",", runCaseNames[j], ",", lastIter, ": hung\n", sep = "")
	      runStatus[[i]][[j]]$status[lastIndex] <- "hung"
	    }
	  }
	}
      }
    }
  }
  runStatus
}
