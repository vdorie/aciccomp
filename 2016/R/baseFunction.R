newContinuousBaseFunction <- function(C, term, name, pars = NULL) {
  rsign <- function() if (rbinom(1L, 1L, 0.5) > 0L) 1 else -1
  
  if (is.null(pars)) {
    pars <- switch(name,
      constant = rt(1L, C$BF_DF) * C$BF_CONSTANT_SCALE,
      linear   = rt(1L, C$BF_DF) * C$BF_LINEAR_SCALE,
      quadratic = {
        a <- 2 * rbeta(1L, C$BF_QUADRATIC_SHAPE_1, C$BF_QUADRATIC_SHAPE_1) - 1
        b <- C$BF_CONTINUOUS_SCALE * C$BF_QUADRATIC_SCALE * rsign() /
             rgamma(1L, C$BF_QUADRATIC_SHAPE_2, C$BF_QUADRATIC_RATE)
        c(a, b)
      },
      cubic     = {
        a <- 2 * rbeta(2L, C$BF_CUBIC_SHAPE, C$BF_CUBIC_SHAPE) - 1
        f <- function(x) x * (x - a[1L]) * (x - a[2L])
        roots <- sum(a) + c(-1, 1) * sqrt(a[1L]^2 - prod(a) + a[2L]^2)
        extrema <- c(f(roots), f(-1), f(1))
        b <- C$BF_CONTINUOUS_SCALE * C$BF_CUBIC_SCALE * rsign() / 
          rgamma(1L, max(abs(extrema)) * C$BF_CUBIC_RATE, C$BF_CUBIC_RATE)
        c(a, b)
      },
      step.constant = c(2 * rbeta(1L, C$BF_STEP_SHAPE, C$BF_STEP_SHAPE) - 1, C$BF_CONTINUOUS_SCALE * C$BF_STEP_CONSTANT_SCALE * rt(1L, C$BF_DF)),
      step.linear   = c(2 * rbeta(1L, C$BF_STEP_SHAPE, C$BF_STEP_SHAPE) - 1, C$BF_CONTINUOUS_SCALE * C$BF_STEP_LINEAR_SCALE * rt(1L, C$BF_DF)),
      sigmoid       = c(rgamma(1L, C$BF_SIGMOID_SHAPE_1, C$BF_SIGMOID_RATE_1), rsign() * rgamma(1L, C$BF_SIGMOID_SHAPE_2, C$BF_SIGMOID_RATE_2)),
      step.quantile = {
        a <- rbeta(1L, C$BF_QUANTILE_SHAPE_1, C$BF_QUANTILE_SHAPE_2) / 2
        b <- rbinom(1L, 1L, 0.5)
        q <- unname(quantile_old(term, if (b > 0) 1 - a else a)) ## b is which side of 0.5
        c(b, q)
      }
    )
  }
  
  f <- switch(name,
    constant = function(x) rep_len(pars[1L], length(x)),
    linear   = function(x) pars[1L] * x,
    quadratic = function(x) pars[2L] * x * (x - pars[1L]),
    cubic = function(x) pars[3L] * x * (x - pars[1L]) * (x - pars[2L]),
    step.constant = if (pars[1L] < 0) function(x) ifelse(x < pars[1L], pars[2L], 0) else function(x) ifelse(x < pars[1L], 0, pars[2L]),
    step.linear   = if (pars[1L] < 0) function(x) ifelse(x < pars[1L], (x - pars[1L]) * pars[2L], 0) else function(x) ifelse(x < pars[1L], 0, (x - pars[1L]) * pars[2L]),
    sigmoid       = if (pars[2L] < 0) function(x) pnorm(-pars[1L] * (x - pars[2L])) else function(x) pnorm(pars[1L] * (x - pars[2L])),
    step.quantile = if (pars[1L] > 0) function(x) ifelse(x >= pars[2L], 1, 0) else function(x) ifelse(x <= pars[2L], 1, 0)
  )
  environment(f) <- new.env(parent = baseenv())
  environment(f)$pars <- pars
  
  new("BaseFunction", type = baseFunctionsMap[[name]], pars = pars, f = f)
}

newDiscreteBaseFunction <- function(C, term, pars = NULL) {
  uniqueMembers <- unique(term)
  numLevels <- length(uniqueMembers)
  if (numLevels <= 1) browser()
  
  if (is.null(pars)) {
    pars <- 2.3 * rt(numLevels - 1L, C$BF_DF) / (3 * sqrt(numLevels - 1))
    zeroIndex <- sample(numLevels, 1L)
    if (zeroIndex > length(pars)) {
      pars <- c(pars, 0)
    } else if (zeroIndex == 1L) {
      pars <- c(0, pars)
    } else {
      pars <- c(pars[seq.int(1L, zeroIndex - 1L)], 0, pars[seq.int(zeroIndex, length(pars))])
    }
  }
  
  parMap <- names(pars)
  if (is.null(parMap)) {
    if (!is.null(levels(uniqueMembers))) {
      parMap <- levels(uniqueMembers)
    } else {
      parMap <- sort(uniqueMembers)
    }
  }
  attr(pars, "map") <- parMap

  
  f <- function(x) pars[match(x, attr(pars, "map"))]
  environment(f) <- new.env(parent = baseenv())
  environment(f)$pars <- pars
  
  new("BaseFunction", type = baseFunctionsMap[["step.discrete"]], pars = pars, f = f)
}

shuffleBaseline <- function(func)
{
  pars <- func@pars
  parMap <- attr(pars, "map")
  
  swapIndex <- sample(length(parMap), 1L)
  temp <- parMap[swapIndex]
  parMap[swapIndex] <- parMap[1L]
  parMap[1L] <- temp

  attr(pars, "map") <- parMap
  environment(func@f)$pars <- pars
  func@pars <- pars
  func
}

rescaleBaseFunction <- function(func, power)
{
  funcType <- baseFunctions[func@type]
  if (funcType %in% c("constant", "linear")) {
   func@pars <- sign(func@pars) * abs(func@pars)^power
  } else if (funcType %in% c("quadratic", "step.constant", "step.linear")) {
    func@pars[2L] <- sign(func@pars[2L]) * abs(func@pars[2L])^power
  } else if (funcType == "cubic") {
    func@pars[3L] <- sign(func@pars[3L]) * abs(func@pars[3L])^power
  }
  environment(func@f)$pars <- func@pars
  func
}

tweakFunction <- function(C, func)
{
  rsign <- function(x) if (runif(1L) < C$BF_TWEAK_SIGN_PROB) 1 else -1
  funcType <- baseFunctions[func@type]
  if (funcType %in% c("constant", "linear")) {
    func@pars <- rsign() * (func@pars + rnorm(1L, 0, C$BF_TWEAK_NORMAL_SCALE))
  } else if (funcType %in% c("quadratic")) {
    func@pars[2L] <- rsign() * (func@pars[2L] + rgamma(1L, C$BF_TWEAK_GAMMA_SHAPE, C$BF_TWEAK_GAMMA_RATE) / 2 - 1)
  } else if (funcType == "cubic") {
    func@pars[3L] <- rsign() * (func@pars[3L] + rgamma(1L, C$BF_TWEAK_GAMMA_SHAPE, C$BF_TWEAK_GAMMA_RATE) / 2 - 1)
  } else if (funcType %in% c("step.constant", "step.linear")) {
    func@pars[2L] <- rsign() * (func@pars[2L] + rnorm(1L, 0, C$BF_TWEAK_NORMAL_SCALE))
  }
  
  if ("tweakFunction" %not_in% C$RUN_BUGGED)
    environment(func@f) <- new.env(parent = baseenv())
  environment(func@f)$pars <- func@pars
  func
}

setConstantForFunction <- function(func, const)
{
  func@pars <- const
  environment(func@f)$pars <- func@pars
  func
}

absBaseFunction <- function(func)
{
  funcType <- baseFunctions[func@type]
  if (funcType %in% c("constant", "linear")) {
    func@pars <- abs(func@pars)
  } else if (funcType %in% c("quadratic", "step.constant", "step.linear")) {
    func@pars[2L] <- abs(func@pars[2L])
  } else if (funcType == "cubic") {
    func@pars[3L] <- abs(func@pars[3L])
  }
  environment(func@f)$pars <- func@pars
  func
}

