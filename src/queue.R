getCurrentJobs <- function() {
  currentJobs <- data.frame(id = integer(0L), status = character(0L), name = character(0L), stringsAsFactors = FALSE)
  
  qstatPath <- suppressWarnings(system2("command", c("-v", "qstat"), stdout = TRUE))
  
  if (length(qstatPath) == 0L) {
    warning("unable to find qstat")
    return(currentJobs)
  }
  
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
  
  currentJobs
}
