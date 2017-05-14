## gets important values from the environment, expects to be run in a directory
## where runCases.csv, methods.csv, runStatus.Rdata, and site_setup.R all exist
method <- Sys.getenv("METHOD")
runCaseName <- strsplit(Sys.getenv("RUNCASE_NAME"), "_")[[1L]]
runCaseName <- paste0(paste0(head(runCaseName, -1L), collapse = "_"), "/", tail(runCaseName, 1L))

source("site_setup.R")

runCases <- read.csv("runCases.csv")
methods <- read.csv("methods.csv", stringsAsFactors = FALSE)
load("runStatus.Rdata")

method <- methods[methods$name == method,]
runStatus <- runStatus[[which(methods$name == method$name)]][[which(names(runStatus[[1L]]) == runCaseName)]]

x.comp <- read.csv(file.path(dirs$data, "x.csv"))

df <- data.frame(z = integer(nrow(x.comp)), y = numeric(nrow(x.comp)))
for (i in seq_along(x.comp)) df[[names(x.comp)[i]]] <- x.comp[[i]]
rm(x.comp)

inFile <- tempfile(fileext = ".csv")

command <- file.path(dirs$methods, method$name)

iters <- runStatus$iter
runCaseDir <- gsub("\\/", .Platform$file.sep, runCaseName)

resultsDir <- file.path(dirs$results, method$name, runCaseDir)
if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)

currentResultIndices <- which(runStatus$status %in% c("complete", "hung"))
if (length(currentResultIndices) == length(iters)) { unlink(inFile); q("no") }

cat("fitting method '", method$name, "' in setting '", runCaseName, "'\n", sep = "")

dataDir <- file.path(dirs$data, runCaseDir)

for (i in seq_along(iters)) {
  if (i %in% currentResultIndices) next
  
  dataFile <- file.path(dataDir, paste0(iters[i], ".csv"))
  
  resp <- read.csv(dataFile, header = FALSE, col.names = c("z", "y"))
  df$z <- resp$z
  df$y <- resp$y
  
  write.table(df, file = inFile, sep = ",", dec = ".", row.names = FALSE, col.names = method$headers_in == 1L)
  
  outFile <- file.path(resultsDir, paste0(iters[i], ".csv"))
  
  args <- c(inFile, outFile)
  
  if (method$individual_effects != 0) {
    outFile.ind <- file.path(resultsDir, paste0(iters[i], "_ind.csv"))
    args <- c(args, outFile.ind)
  }
  
  cat(gsub("/", "_", runCaseName), ",", iters[i], ",", sep = "")
  startTime <- proc.time()
  consoleLog <- system2(command, args, stdout = TRUE, stderr = TRUE)
  timeDiff <- proc.time() - startTime
  
  if (!is.null(attr(consoleLog, "status"))) {
    cat("error in: '", runCaseName, "' iter ", iters[i], "\n", sep = "", file = stderr())
    cat(consoleLog, sep = "\n", file = stderr())
    cat("\n\n", file = stderr())
    
    cat("NA\n", sep = "")
    result <- data.frame(est = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_)
    write.table(result, file = outFile, sep = ",", dec = ".", row.names = FALSE, col.names = method$headers_out == 1L)
    
    if (method$individual_effects != 0L) {
      result <- data.frame(
        est =  rep_len(NA_real_, nrow(df)),
        ci_lower = rep_len(NA_real_, nrow(df)),
        ci_upper = rep_len(NA_real_, nrow(df)))
      
      write.table(result, file = outFile.ind, sep = ",", dec = ".", row.names = FALSE, col.names = method$headers_out == 1L)
    }
  } else {
    cat(timeDiff[["sys.child"]] + timeDiff[["user.child"]], "\n", sep = "")
  }
}

unlink(inFile)
