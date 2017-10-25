transformData <- function(x)
{
  ## function definitions
  boxcoxTransform <- function(y)
  {
    boxcox <- function(x, lambda) 
      if (lambda == 0) log(x) else (x^lambda - 1) / lambda
    
    boxcox.likelihood <- function(y, lambda, mu, sigma, log = TRUE) {
      res <- if (lambda == 0)
        sum(dnorm(log(y), mu, sigma, TRUE))
      else
        res <- sum(dnorm((y^lambda - 1) / lambda, mu, sigma, TRUE))
      
      if (log) res else exp(res)
    }
    
    boxcox.profile <- function(y, lambda) {
      n <- length(y)
      y.sum.log <- sum(log(y))
      
      sapply(lambda, function(lambda.i) {
        y.lambda <- if (lambda.i == 0)
          log(y)
        else
          (y^lambda.i - 1) / lambda.i
        
        ssr <- sum((y.lambda - mean(y.lambda))^2)
        n * (log(2 * pi) + log(ssr / n) + 1) - 2 * (lambda.i - 1) * y.sum.log
      })
    }
    
    wrapper <- function(lambda) boxcox.profile(y, lambda)
    optimResults <- optim(0, wrapper, method = "L-BFGS-B", lower = -8, upper = 8)
    boxcox(y, optimResults$par)
  }
  
  x_orig <- x; rm(x)
  
  standardize <- function(x) (x - mean(x)) / sd(x)
  
  factorColumns <- evalx(colnames(x_orig), x[sapply(x, function(y) is.factor(x_orig[[y]]))])
  
  ## rely on quantile picking out order statistics, and these columns bunching up repeatedly at ends
  percentExtreme <- sapply(setdiff(colnames(x_orig), factorColumns), function(x) {
    evalx(x_orig[[x]], mean(x <= quantile_old(x, 0.1) | x >= quantile_old(x, 0.9)))
  })
  continuousColumns <- evalx(percentExtreme, names(x)[abs(x) <= 0.35])
  
  discreteColumns <- setdiff(colnames(x_orig), c(factorColumns, continuousColumns))
  
  cutoffs <- sapply(discreteColumns, function(y) evalx(x_orig[[y]], {
    q <- quantile_old(x, c(0.25, 0.5, 0.75))
    res <- unname(q[if (q[2L] - q[1L] <= q[3L] - q[2L]) 1L else 3L])
    ## cutoff is >= res, but we do > res - 1 for consistency
    if (res > 0L) res - 1L else res
  }))
  
  x_trans <- x_orig
  
  for (colName in continuousColumns)
    x_trans[[colName]] <- standardize(boxcoxTransform(evalx(x_trans[[colName]], if (any(x == 0L)) x + 1L else x)))
  
  for (colName in discreteColumns) {
    if (cutoffs[colName] == 0L) {
      x_trans[[colName]] <- factor(x_trans[[colName]] != 0, labels = c("leq_0", "gt_0"))
    } else {
      x_trans[[colName]] <- factor(x_trans[[colName]] > cutoffs[colName], labels = paste0(c("leq_", "gt_"), cutoffs[colName]))
    }
  }
  
  x_trans
}

