queueJobs <- function(runStatus, methods, dirs, runMethods = NULL, dryRun = FALSE)
{
  if (length(runStatus) == 0L) return(invisible(NULL))
  
  runCaseNames <- names(runStatus[[1L]])
  
  if (!dir.exists(dirs$job)) dir.create(dirs$job, recursive = TRUE)
  
  for (i in seq_len(nrow(methods))) {
    method <- methods[i,]
    if (!is.null(runMethods) && !(method$name %in% runMethods)) next
    
    if (file.exists(file.path(dirs$results, paste0(method$name, ".tar.gz")))) next
    
    moduleString <- "module load r/intel/3.2.2"
    if (method$language == "python") {
      moduleString <- paste0(moduleString, "; ",
                             "module switch gcc gcc/4.9.2; ",
                             "module load scikit-learn/intel/0.17b1")
    } else if (method$language == "stata") {
      moduleString <- paste0(moduleString, "; ",
                             "module load stata/14.0")
    }
    
    for (j in seq_along(runCaseNames)) {
      if (all(runStatus[[i]][[j]]$status %in% c("complete", "hung"))) next
      
      jobName <- paste0(method$name, "_", runCaseNames[j])
      
      jobFile <- file.path(dirs$job, paste0(jobName, ".pbs"))
      
      args <- paste0("'s|_METHOD_|", method$name, "|;s|_RUNCASE_|", j, "|;",
                     "s|_DATA_DIR_|", dirs$data, "|;s|_METHOD_DIR_|", dirs$methods,
                     "s|_JOBNAME_|", jobName, "|;s|_MODULE_|", moduleString, "|;'")
      system2("sed", args,
              stdin  = "template.pbs",
              stdout = jobFile)
      if (dryRun)
        cat("queueing ", jobName, ":\n", sep = "")
      print(format(subset(runStatus[[i]][[j]], status != "complete" & status != "hung")))
      if (!dryRun)
        system2("qsub", jobFile)
      cat("\n")
    }
  }
}
