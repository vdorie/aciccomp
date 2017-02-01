## gets important values from the environment, expects to be run in a directory
## where runCases.csv, methods.csv, and runStatus.Rdata all exist
method <- Sys.getenv("METHOD")
runCaseName <- Sys.getenv("RUNCASE")

dataDir <- Sys.getenv("DATA_DIR")
methodDir <- Sys.getenv("METHOD_DIR")

runCases <- read.csv("runCases.csv")
methods <- read.csv("methods.csv", stringsAsFactors = TRUE)
load("runStatus.Rdata")

runStatus <- runStatus[[which(methods$name == method$name)]][[which(names(runStatus[[1L]]) == runCaseName)]]

x.comp <- read.csv(file.path(dataDir, "x.csv"))

df <- data.frame(z = integer(nrow(x.comp)), y = numeric(nrow(x.comp)))
for (i in seq_along(x.comp)) df[[names(x.comp)[i]]] <- x.comp[[i]]
rm(x.comp)

inFile <- tempfile(fileext = ".csv")

method <- methods[methods$name == method,]
command <- file.path(methodDir, method$name)

iters <- runStatus$iter
runCaseDir <- gsub("\\/", .Platform$file.sep, runCaseName)


resultsDir <- file.path(resultsDir, method$name, runCaseDir)
if (!dir.exists(resultsDir)) dir.create(resultsDir, recursive = TRUE)

currentResultIndices <- which(runStatus$status %in% c("complete", "hung"))
if (length(currentResultIndices) == length(iters)) { unlink(inFile); q("no") }

cat("fitting method '", method$name, "' in setting '", runCaseName, "'\n", sep = "")

dataDir <- file.path(dataDir, runCaseDir)

for (i in seq_along(iters)) {
  if (i %in% currentResultIndices) next
  
  dataFile <- file.path(dataDir, paste0(iters[i], ".csv"))
  
  resp <- read.csv(dataFile)
  df$z <- resp$z
  df$y <- resp$y
  #df$y <- ifelse(resp$z == 0L, resp$y_0, resp$y_1)
  write.csv(df, file = inFile, row.names = FALSE)
  
  outFile <- file.path(resultsDir, paste0(iters[i], ".csv"))
  
  args <- c(inFile, outFile)
  
  if (method$individual_effects != 0) {
    outFile.ind <- file.path(resultsDir, paste0(iters[i], "_ind.csv"))
    args <- c(args, outFile.ind)
  }
  
  cat(setting, ",", iters[i], ",", sep = "")
  startTime <- proc.time()
  consoleLog <- system2(command, args, stdout = TRUE, stderr = TRUE)
  timeDiff <- proc.time() - startTime
  
  if (!is.null(attr(consoleLog, "status"))) {
    cat("error in: '", runCaseName, "' iter ", iters[i], "\n", sep = "", file = stderr())
    cat(consoleLog, sep = "\n", file = stderr())
    cat("\n\n", file = stderr())
    
    cat("NA\n", sep = "")
    write.csv(data.frame(est = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_), file = outFile, row.names = FALSE)
    
    if (method$individual_effects != 0L) {
      result <- data.frame(
        est =  rep_len(NA_real_, nrow(df)),
        ci_lower = rep_len(NA_real_, nrow(df)),
        ci_upper = rep_len(NA_real_, nrow(df)))
      
      write.csv(result, file = outFile.ind, row.names = FALSE)
    }
  } else {
    cat(timeDiff[["sys.child"]], "\n", sep = "")
  }
}

unlink(inFile)
