dgp_2016 <- function(x, parameters, random.seed, constants = constants_2016(), extraInfo = FALSE)
{
  C <- constants
  C$RUN_BUGGED <- ""
  
  if (is.matrix(x) || is.numeric(x)) x <- as.data.frame(x)
  if (!is.data.frame(x)) stop("'x' argument must be coercible to data.frame")
  if (NROW(x) == 0L || NCOL(x) == 0L) stop("'x' argument must have positive dimensions")
  
  if (is.numeric(parameters)) {
    ## parameter is a run-case number
    #data("parameters_2016", package = "aciccomp2016")
    
    if (length(parameters) != 1L || parameters < 1L || parameters > nrow(aciccomp2016::parameters_2016))
      stop("numeric 'parameters' argument must specify row of 'parameters_2016'")
    runCaseNumber <- parameters
    parameters <- aciccomp2016::parameters_2016[runCaseNumber,]
    
    if (length(random.seed) == 1L) {
      ## random.seed is an iteration number
      
      if (random.seed < 1L || random.seed > nrow(randomSeeds[[runCaseNumber]]))
        stop("random.seed for parameters_2016 must be in [1, ", nrow(randomSeeds[[runCaseNumber]]), "]")
      
      runCaseIter  <- random.seed
      C$RUN_BUGGED <- strsplit(as.character(randomSeeds[[runCaseNumber]][runCaseIter, "run.bugged"]), ":")[[1L]]
      random.seed  <- randomSeeds[[runCaseNumber]][runCaseIter, "random.seed"]
    }
  }
  ## in case parameters is a row of a data.frame, convert it into an ordinary named list
  parameters <- as.list(parameters)
  evalx(parameters[sapply(parameters, is.factor)], x <- sapply(x, as.character))
  
  if (is.integer(random.seed) && length(random.seed) > 1L) {
    ## random.seed is .Random.seed the object
    .GlobalEnv$.Random.seed <- random.seed
  } else {
    if (is.numeric(random.seed)) random.seed <- list(seed = random.seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rounding")
    
    if (!is.list(random.seed))
      stop("random.seed must be a .Random.seed integer vector, an integer scalar, or a list containing the arguments to set.seed")
    
    if (any(names(random.seed) == "seed")) {
      seedPosition <- which.max(names(random.seed) == "seed")
    } else if (any(names(random.seed == ""))) {
      seedPosition <- which.max(names(random.seed) == "")
    } else {
      stop("random.seed does not conform to formals of set.seed")
    }
    
    ## seed manipulation in original
    random.seed[[seedPosition]] <- random.seed[[seedPosition]] * 5L + if (random.seed[[seedPosition]] <= 500L) 565L else 7565L
    
    suppressWarnings(do.call("set.seed", random.seed))
  }
  
  if (any(c("model.trt", "root.trt", "overlap.trt", "model.rsp", "alignment", "te.hetero") %not_in% names(parameters)))
    stop("parameters must be named list with members 'model.trt', 'root.trt', 'overlap.trt', 'model.rsp', 'alignment', and 'te.hetero'")
  
  model.trt   <- parameters$model.trt
  root.trt    <- parameters$root.trt
  overlap.trt <- parameters$overlap.trt
  model.rsp   <- parameters$model.rsp
  alignment   <- parameters$alignment
  te.hetero   <- parameters$te.hetero
  
  x.orig <- x
  x <- transformData(x.orig)
  
  isContinuousColumn <- sapply(x, function(x.i) !is.discrete(x.i))
  isSparseDiscreteColumn <- sapply(seq_along(x), function(i) !isContinuousColumn[i] && { tab <- table(x[[i]]) / nrow(x); length(tab) == 2L && max(tab) >= 0.95 })
  
  continuousColumnNames <- colnames(x)[isContinuousColumn]
  
  termProb <- rep_len(C$DEFAULT_COVARIATE_WEIGHT, ncol(x))
  termProb[isSparseDiscreteColumn] <- C$SPARSE_COVARIATE_WEIGHT
  termProb[ isContinuousColumn] <- C$CONTINUOUS_COVARIATE_WEIGHT
  termProb <- termProb / sum(termProb)
  
  x.bias <- x[,colnames(x) %in% continuousColumnNames]
  
  C$dist.lin@termInclusionProb[[2L]] <- termProb
  C$dist.int@termInclusionProb[[2L]] <- termProb
  C$dist.poly@termInclusionProb[[2L]] <- termProb
  C$dist.pure.poly@termInclusionProb[[2L]] <- termProb
  C$dist.step@termInclusionProb[[2L]] <- termProb
  C$dist.exp@termInclusionProb[[2L]] <- termProb
  
  dist.trt <- if (is.character(model.trt)) {
    switch(model.trt,
      linear          = C$dist.lin,
      interaction     = C$dist.int,
      polynomial      = C$dist.poly,
      pure.polynomial = C$dist.pure.poly,
      step            = C$dist.step,
      test            = C$dist.test)
  } else model.trt
  
  C$BF_DF <- C$TRT_BF_DF
  
  gaf.trt <- generateGeneralizedAdditiveFunction(C, x, "identity", c(1, 0), ## (1, 0) gives one 1 additive function
                                                 c(C$TRT_INPUT_SCALE, C$TRT_OUTPUT_SCALE),
                                                 dist.trt)
  baselineScale <- C$TRT_BASELINE_SHIFT(root.trt)
  gaf.trt <- setBaselineForGeneralizedAdditiveFunction(C, gaf.trt, x, qlogis(rnorm(1L, root.trt, baselineScale)))
  
  p.score.raw <- evaluate(gaf.trt, x)
  p.score.raw.mean <- mean(p.score.raw)
  p.score.raw.sd <- sd(p.score.raw)
  ## use a beta-prime
  p.score.scale.theta <- rgamma(1L, C$TRT_LINEAR_SCALE_SHAPE_1, C$TRT_LINEAR_SCALE_RATE)
  p.score.scale <- rgamma(1L, C$TRT_LINEAR_SCALE_SHAPE_2 * p.score.scale.theta, p.score.scale.theta)
  
  if (overlap.trt != "full") {
    if (overlap.trt != "two-term") {
      af.bias <- generateBiasingFunction3(C, x.bias, c(1, 10 * C$TRT_BIAS_SCALE * p.score.raw.sd), C$dist.bias1)
      gaf.trt <- addAdditiveFunctionToGeneralizedAdditiveFunction(gaf.trt, af.bias)
    } else {
      af.bias1 <- generateBiasingFunction3(C, x.bias, c(1, C$TRT_BIAS_SCALE * p.score.raw.sd), C$dist.bias2)
      af.bias2 <- generateBiasingFunction3(C, x.bias, c(1, C$TRT_BIAS_SCALE * p.score.raw.sd), C$dist.bias2)
      gaf.trt <- addAdditiveFunctionToGeneralizedAdditiveFunction(gaf.trt, af.bias1, weight = 1 / sqrt(2))
      gaf.trt <- addAdditiveFunctionToGeneralizedAdditiveFunction(gaf.trt, af.bias2, weight = 1 / sqrt(2))
    }
  }
  gaf.trt <- scaleGeneralizedAdditiveFunction(gaf.trt, p.score.raw.mean, p.score.scale / p.score.raw.sd)
  
  p.score.raw <- evaluate(gaf.trt, x)

  p.score <- plogis(p.score.raw)
  z <- factor(rbinom(nrow(x), 1L, p.score), labels = c("ctl", "trt"))
  
  dist.rsp <- switch(model.rsp,
    linear          = C$dist.lin,
    interaction     = C$dist.int,
    polynomial      = C$dist.poly,
    pure.polynomial = C$dist.pure.poly,
    step            = C$dist.step,
    exponential     = C$dist.poly,
    test            = C$dist.test)
  
  C$BF_DF <- C$RSP_BF_DF
  
  if (alignment > 0) {
    gaf.rsp <- generateGeneralizedAdditiveFunctionFromPrototypeAndForceOverlapAlignment(
      C, x, "identity", c(1, 0),
      c(C$RSP_INPUT_SCALE, C$RSP_OUTPUT_SHAPE_2),
      dist.rsp, gaf.trt, alignment)
  } else {
    gaf.rsp <- generateGeneralizedAdditiveFunction(C, x, "identity", c(1, 0),
                                                   c(C$RSP_INPUT_SCALE, C$RSP_OUTPUT_SHAPE_2),
                                                   dist.rsp)
  }
  
  gaf.rsp <- setBaselineForGeneralizedAdditiveFunction(C, gaf.rsp, x, rt(1L, C$BF_DF))

  if (model.rsp == "exponential") {
    gaf.exp <- if (alignment > 0)
        generateGeneralizedAdditiveFunctionFromPrototype(C, x, "identity", c(1, 0), c(C$RSP_INPUT_SCALE, 1),
                                                         C$dist.exp, gaf.trt, alignment / 3)
      else
        generateGeneralizedAdditiveFunction(C, x, "identity", c(1, 0), c(C$RSP_INPUT_SCALE, 1), C$dist.exp)
    
    if (length(gaf.exp@functions) != 1L || length(gaf.exp@functions[[1L]]@terms) != 0L) {
      gaf.exp@functions[[1L]] <- setBaselineForAdditiveFunction(C, gaf.exp@functions[[1L]], x, 0)
      
      expScale <- rgamma(1L, C$RSP_EXP_SCALE_SHAPE, C$RSP_EXP_SCALE_RATE)
      gaf.exp@functions[[1L]]@outputScale <- expScale / sd(evaluate(gaf.exp, x))
      gaf.exp@functions[[1L]]@transformation@type <- transformationsMap[["exp"]]
      
      weights <- rep_len(sd(evaluate(gaf.rsp, x)) * rgamma(1L, C$RSP_EXP_WEIGHT_SHAPE, C$RSP_EXP_WEIGHT_RATE) /
                         sd(evaluate(gaf.exp, x)), length(gaf.exp@weights))
      if (runif(1L) < 0.5) weights <- -weights
      gaf.exp@weights <- as.list(weights)
      
      gaf.rsp <- combineGafs(gaf.rsp, gaf.exp)
    }
  }
  
  if (te.hetero == "med") {
    gaf.rsp <- addTreatmentToGeneralizedAdditiveFunction(C, gaf.rsp, x, z, C$dist.hetero.med)
  } else if (te.hetero == "high") {
    gaf.rsp <- addTreatmentToGeneralizedAdditiveFunction(C, gaf.rsp, x, z, C$dist.hetero.high)
  }
  
  x$.z <- factor(rep_len(0L, nrow(x)), levels = c(0L, 1L), labels = c("ctl", "trt"))
  mu.0 <- evaluate(gaf.rsp, x)
  x$.z <- factor(rep_len(1L, nrow(x)), levels = c(0L, 1L), labels = c("ctl", "trt"))
  mu.1 <- evaluate(gaf.rsp, x)
  
  mu.mean <- c(mean(mu.0), mean(mu.1))
  mu.sd <- c(sd(mu.0), sd(mu.1))
  mu.theta <- rgamma(1L, C$RSP_OUTPUT_SHAPE_1, C$RSP_OUTPUT_RATE)
  mu.scale <- rgamma(1L, C$RSP_OUTPUT_SHAPE_2 * mu.theta, mu.theta)
  
  mu.0 <- mu.scale * (mu.0 - mu.mean[1L]) / mu.sd[1L] + mu.mean[1L]
  mu.1 <- mu.scale * (mu.1 - mu.mean[2L]) / mu.sd[2L] + mu.mean[2L]
  
  estSatt <- mean(mu.1[z == "trt"] - mu.0[z == "trt"])
  trueAtt <- C$RSP_TE_MEAN + C$RSP_TE_SCALE * rt(1L, C$RSP_TE_DF)
  
  mu.1 <- mu.1 + (mu.scale * trueAtt - estSatt)
  
  gaf.rsp <- addTreatmentMainEffectToGeneralizedAdditiveFunction3(C, gaf.rsp, z,
    -mu.mean[1L] + (mu.sd[1L] / mu.scale) * mu.mean[1L],
    -mu.mean[2L] + (mu.sd[2L] / mu.scale) * (mu.mean[2L] + mu.scale * trueAtt - estSatt))
  
  weightFunc <- function(x) scale * ifelse(x$.z == "ctl", mu.scale / mu.sd[1L], mu.scale / mu.sd[2L])
  for (i in seq_along(gaf.rsp@weights)) {
    oldWeight <- gaf.rsp@weights[[i]]
    gaf.rsp@weights[[i]] <- weightFunc
    environment(gaf.rsp@weights[[i]]) <- new.env(parent = baseenv())
    environment(gaf.rsp@weights[[i]])$mu.scale <- mu.scale
    environment(gaf.rsp@weights[[i]])$mu.sd <- mu.sd
    environment(gaf.rsp@weights[[i]])$scale <- oldWeight
  }
  
  y.0 <- rnorm(nrow(x), mu.0, C$RSP_SIGMA_Y)
  y.1 <- rnorm(nrow(x), mu.1, C$RSP_SIGMA_Y)
  
  y <- ifelse(z == "trt", y.1, y.0)
  
  result <- list(z = z, y = y, y.0 = y.0, y.1 = y.1, mu.0 = mu.0, mu.1 = mu.1, e = p.score)
  if (extraInfo) {
    result$f.z <- gaf.trt
    result$f.y <- gaf.rsp
    result$x <- x
    result$x$.z <- NULL
    result$valid <- isValidSim(parameters, x.orig, result)
  } else {
    result$z <- as.integer(result$z) - 1L
  }
  
  result
}

isValidSim <- function(parameters, data, result) {
  treatedRows <- result$z == "trt"
  
  ## overlap as percentage of controls less than min treated p.score
  overlap <- mean(result$e[!treatedRows] < min(result$e[treatedRows]))
  
  ## spearman correlation between linear propensity score and response
  corr.trt.rsp <- cor(qlogis(result$e), result$y, method = "spearman")
  
  pste <- sd(result$mu.1 - result$mu.0)
  
  ## calculate true att
  satt <- mean(result$y.1[treatedRows] - result$y.0[treatedRows])
  #satt.z <- satt / sd(result$y)
  #pate <- mean(result$y.1 - result$y.0)
  
  ## calculate att with untransformed design matrix
  data$z <- result$z
  data$y <- result$y
  fit <- lm(y ~ ., data)
  
  x.pred <- data[treatedRows,]
  x.pred$z <- factor(rep_len(0L, nrow(x.pred)), levels = c(0L, 1L), labels = c("ctl", "trt"))
  x.pred$y <- NULL
  y.hat.0 <- predict(fit, x.pred)
  x.pred$z <- factor(rep_len(1L, nrow(x.pred)), levels = c(0L, 1L), labels = c("ctl", "trt"))
  y.hat.1 <- predict(fit, x.pred)
  
  satt.lm <- mean(y.hat.1 - y.hat.0)
  satt.z.lm <- satt.lm / sd(result$y)
  bias.z.lm <- (satt - satt.lm) / sd(result$y)
  
  rsq.lm <- summary(fit)$r.squared
  
  
  isValidResult <- with(parameters,
    (overlap.trt == "full" || overlap > 0.2) &&
    ((model.trt == "linear" && model.rsp == "linear") || te.hetero == "none" || alignment == 0 || abs(bias.z.lm) > 0.1) &&
    ((te.hetero == "high" && pste > 2) || (te.hetero == "med" && pste > 1) || te.hetero == "none") &&
    (model.rsp == "linear" || rsq.lm < 0.8) &&
    (alignment <= 0.25 || abs(corr.trt.rsp) > 0.15)
  )
  
  isValidResult
}

