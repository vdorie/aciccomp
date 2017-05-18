queueJobs <- function(runStatus, methods, dirs, runMethods = NULL, dryRun = FALSE)
{
  if (length(runStatus) == 0L) return(invisible(NULL))
  
  runCaseNames <- names(runStatus[[1L]])
  
  if (!dir.exists(dirs$job)) dir.create(dirs$job, recursive = TRUE)
  if (!dir.exists(dirs$log)) dir.create(dirs$log, recursive = TRUE)
  
  for (i in seq_len(nrow(methods))) {
    method <- methods[i,]
    if (!is.null(runMethods) && length(runMethods) > 0L && !(method$name %in% runMethods)) next
    
    if (file.exists(file.path(dirs$results, paste0(method$name, ".tar.gz")))) next
    
    if (!dir.exists(file.path(dirs$job, method$name))) dir.create(file.path(dirs$job, method$name))
    if (!dir.exists(file.path(dirs$log, method$name))) dir.create(file.path(dirs$log, method$name))
    
    for (j in seq_along(runCaseNames)) {
      if (all(runStatus[[i]][[j]]$status %in% c("complete", "hung"))) next
      
      runCaseName <- gsub("/", "_", runCaseNames[j])
      jobName <- paste0(method$name, "_", runCaseName)
      
      jobFile <- file.path(dirs$job, method$name, paste0(runCaseName, ".pbs"))
      
      args <- c("-r", "-e",
                paste0("'s|_METHOD_|", method$name, "|;s|_RUNCASE_NAME_|", runCaseName, "|;",
                       "s|_JOBNAME_|", jobName, "|'"))
      system2("sed", args,
              stdin  = "template.pbs",
              stdout = jobFile)
      cat("queueing ", jobName, ":\n", sep = "")
      with(runStatus[[i]][[j]], cat("  ", sum(status == "NA"), " NA, ", sum(status == "missing"), " missing, ", sum(status == "failed"), " failed\n", sep = ""))
      if (!dryRun)
        system2("qsub", c("-C", "'#PBS'", "-q", "short", jobFile))
    }
  }
}
