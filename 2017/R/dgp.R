dgp_2017 <- function(parameters, random.seed)
{
  if (is.numeric(parameters)) {
    ## parameter is a run-case number
    
    if (length(parameters) != 1L || parameters < 1L || parameters > nrow(aciccomp2017::parameters_2017))
      stop("numeric 'parameters' argument must specify row of 'parameters_2017'")
    runCaseNumber <- parameters
    parameters <- aciccomp2017::parameters_2017[runCaseNumber,]
    
    if (length(random.seed) == 1L) {
      ## random.seed is an iteration number
      
      if (random.seed < 1L || random.seed > length(randomSeeds[[runCaseNumber]]))
        stop("random.seed for parameters_2017 must be in [1, ", length(randomSeeds[[runCaseNumber]]), "]")
      
      runCaseIter  <- random.seed
      random.seed  <- randomSeeds[[runCaseNumber]][[runCaseIter]]
    }
  }
  ## in case parameters is a row of a data.frame, convert it into an ordinary named list
  parameters <- as.list(parameters)
  factorPars <- sapply(parameters, is.factor)
  parameters[factorPars] <- sapply(parameters[factorPars], as.character)
  #evalx(parameters[sapply(parameters, is.factor)], x <- sapply(x, as.character))
  
  if (is.integer(random.seed) && length(random.seed) > 1L) {
    ## random.seed is .Random.seed the object
    .GlobalEnv$.Random.seed <- random.seed
  } else {
    if (is.numeric(random.seed)) random.seed <- list(seed = random.seed, kind = "Mersenne-Twister", normal.kind = "Inversion")
    
    if (!is.list(random.seed))
      stop("random.seed must be a .Random.seed integer vector, an integer scalar, or a list containing the arguments to set.seed")
    if (any(names(random.seed) == "seed")) {
      seedPosition <- which.max(names(random.seed) == "seed")
    } else if (any(names(random.seed == ""))) {
      seedPosition <- which.max(names(random.seed) == "")
    } else {
      stop("random.seed does not conform to formals of set.seed")
    }
    
    do.call("set.seed", random.seed)
  }
  
  error       <- parameters$error
  magnitude   <- parameters$magnitude
  noise       <- parameters$noise
  confounding <- parameters$confounding
  
  ## we cache just the columns of the data set that we actually use, since the transformation isn't
  ## straightforward from the data given to contestants
  
  ## colsInUse <- c(1, 3, 10, 14, 15, 21, 24, 43)
  x <- transformedData_2017
  n <- nrow(x)
  
  effectSize <- if (magnitude <= 0) 1 / 3 else 2
  if (confounding <= 0) {
    beta0 <- 0
    beta1 <- 0.5
  } else {
    beta0 <- -1
    beta1 <- 3
  }
  snr <- if (noise <= 0) 0.25 else 1.25
  
  reffectv <- if (error %in% "group_corr") 0.1 else 0.0
  
  p <- function(x) {
    result <- x[,"x_1"] + x[,"x_43"] + 0.3 * (2 - as.numeric(x[,"x_10"]))
    1 / (1 + exp(beta0 + beta1 * result))
  }
  
  mu <- function(x) -sin(qnorm(p(x))) + x[,"x_43"]
  
  alpha <- function(x) {
    result <- as.numeric(x[,"x_3"] == 'leq_0') * as.numeric(x[,"x_24"] == 'B') + (2 - as.numeric(x[,"x_14"])) - (2 - as.numeric(x[,"x_15"]))
    effectSize * result
  }
  
  linpart <- mu(x) + p(x) * alpha(x)
  sigma_y <- snr * sd(linpart)
   
  if (error %in% c("nonadditive")) {
    b <- sqrt(sigma_y^2 + var(linpart))
    a <- mean(linpart)
    muvec1 <- (mu(x) + alpha(x) - a) / (1.25 * b)
    muvec0 <- (mu(x) - a) / (1.25 * b)
    
    alphavec <- 13 * pnorm(muvec1 / sqrt(sigma_y^2 / (1.25 * b)^2 + 1)) - 13 * pnorm(muvec0 / sqrt(sigma_y^2/(1.25 * b)^2 + 1))
    dgp <- data.frame(alpha = alphavec, mu = rep(0, n))
  } else {
    dgp <- data.frame(alpha = alpha(x), mu = mu(x))
  }
  
  het <- if (error %in% "heteroskedastic") seq(0.4, 1.4, length.out = 16) else rep(1, 16)
  
  ## part below random, above not
  z <- rbinom(n, 1, p(x))
  y.temp <- mu(x) + het[as.numeric(x[,"x_21"])] * (1 - reffectv) * sigma_y * rnorm(n)
  y.obs <- y.temp +       z * alpha(x)
  y.cf  <- y.temp + (1 - z) * alpha(x)
  
  ref <- reffectv * sigma_y * rnorm(nlevels(x[,"x_21"]))
  
  y.obs <- y.obs + ref[as.numeric(x[,"x_21"])]
  y.cf  <- y.cf  + ref[as.numeric(x[,"x_21"])]
  
  if (error %in% "nonadditive") {
    eps <- rnorm(n, 0, 0.01)
    y.obs <- 13 * pnorm(y.obs, a, 1.25 * b) - 6 + eps
    y.cf  <- 13 * pnorm(y.cf,  a, 1.25 * b) - 6 + eps
  } 
  
  #data.frame(z = z, y = y.obs, y.cf = y.cf, alpha = dgp$alpha)
  data.frame(z = z, y = y.obs, alpha = dgp$alpha)
}


if (FALSE) {
## this re-generates data and subsetting procedure
X_2017 <- read.csv("~/Downloads/ACIC_debrief/contest_data/X.csv", header = FALSE)

newData <- matrix(sapply(X_2017, as.integer), nrow = nrow(X_2017))
for (col in c(17L, 22L, 38L)) newData[,col] <- newData[,col] - 1L
oldData <- matrix(sapply(aciccomp2016::input_2016, as.integer), nrow = nrow(aciccomp2016::input_2016))

matches <- vector("list", nrow(newData))
for (i in seq_len(nrow(newData))) {
  newRow <- newData[i,]
  matches[[i]] <- which(apply(oldData, 1L, function(oldRow) all(oldRow == newRow)))
}
oldToNewMap <- sapply(matches, function(x) x)

input_2017 <- X_2017
colnames(input_2017) <- paste0("x_", seq_len(ncol(input_2017)))

save(input_2017, file = "~/Repositories/aciccomp/2017/data/input_2017.RData")

transformedData_2017 <- aciccomp2016:::transformData(aciccomp2016::input_2016)[oldToNewMap,c(1, 3, 10, 14, 15, 21, 24, 43)]

}

if (FALSE) {
## checks that all of the files match explictly
setting <- 1L
rootDir <- "~/Downloads/ACIC_debrief/contest_data"
for (error in c("group_corr", "heteroskedastic", "iid", "nonadditive")) {
  for (magnitude in c(0, 1)) {
    for (noise in c(0, 1)) {
      for (confounding in c(0, 1)) {
        dgp.file <- read.csv(file.path(rootDir, error, paste0(magnitude, noise, confounding), "dgp.csv"))
        for (iter in seq_len(250L)) {
          dataFile <- file.path(rootDir, error, paste0(magnitude, noise, confounding), paste0(iter, ".csv"))
          dataHasHeaders <- grepl("['\"]z['\"]\\s*,\\s*['\"]y['\"]", readLines(dataFile, n = 1L), perl = TRUE)
          data.file <- if (dataHasHeaders) read.csv(dataFile) else read.csv(dataFile, header = FALSE, col.names = c("z", "y"))
          
          
          data <- aciccomp2017::dgp_2017(setting, iter)
          if (mean(abs(data$y - data.file$y)) > 1e-10 || any(data$z != data.file$z) || mean(abs(data$alpha - dgp.file$alpha)) > 1e-10)
            stop("mismatch for parameters ", error, " ", magnitude, noise, confounding, ": ", iter, "\n")
        }
        setting <- setting + 1L
      }
    }
  }
}

## test that individual effects match dgp$alpha
setting <- 1L
rootDir <- "~/Downloads/ACIC_debrief/contest_data"
for (error in c("group_corr", "heteroskedastic", "iid", "nonadditive")) {
  for (magnitude in c(0, 1)) {
    for (noise in c(0, 1)) {
      for (confounding in c(0, 1)) {
        dgp.file <- read.csv(file.path(rootDir, error, paste0(magnitude, noise, confounding), "dgp.csv"))
        for (iter in seq_len(250L)) {
          dataFile <- file.path(rootDir, error, paste0(magnitude, noise, confounding), paste0(iter, ".csv"))
          dataHasHeaders <- grepl("['\"]z['\"]\\s*,\\s*['\"]y['\"]", readLines(dataFile, n = 1L), perl = TRUE)
          data.file <- if (dataHasHeaders) read.csv(dataFile) else read.csv(dataFile, header = FALSE, col.names = c("z", "y"))
          
          
          data <- aciccomp2017::dgp_2017(setting, iter)
          if (mean(abs(data$y - data.file$y)) > 1e-10 || any(data$z != data.file$z) || mean(abs(data$alpha - dgp.file$alpha)) > 1e-10)
            stop("mismatch for parameters ", error, " ", magnitude, noise, confounding, ": ", iter, "\n")
        }
        setting <- setting + 1L
      }
    }
  }
}
}

