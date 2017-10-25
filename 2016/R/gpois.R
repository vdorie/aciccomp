rgpois <- function(n, m, v)
{
  if (any(m <= 0)) stop("mean must be positive")
  if (any(v < 0)) stop("variance must be non-negative")
  
  n <- as.integer(n)
  m <- rep_len(m, n)
  v <- rep_len(v, n)
  
  result <- rep(NA_integer_, n)
  zero.v <- v == 0
  if (any(zero.v)) result[zero.v] <- m[zero.v]
  
  n.r <- n - sum(zero.v)
  if (n.r == 0L) return(result)
  
  m.r <- m[!zero.v]
  v.r <- v[!zero.v]
  
  
  lambda <- 1 - sqrt(m.r / v.r)
  theta  <- m.r * (1 - lambda)
  
  result[!zero.v] <- sapply(seq_len(n.r), function(i) {
    if (lambda[i] < 0) {
      q <- floor(-theta[i] / lambda[i])
      
      y.vals <- seq.int(1L, q)
      probs <- c(exp(-theta[i]), theta[i] * (theta[i] + lambda[i] * y.vals)^(y.vals - 1) / factorial(y.vals) * exp(-theta[i] - lambda[i] * y.vals))
      y <- sample(seq.int(0L, q), 1L, prob = probs / sum(probs))
    } else {
      u <- runif(1L)
      p <- exp(-theta[i])
      y <- 0L
      while (p <= u) {
        y <- y + 1L
        p <- p + theta * (theta + lambda * y)^(y - 1) / factorial(y) * exp(-theta - lambda * y)
      }
    }
    y
  })
  
  result
}

# n is population size, p is:
#   [0, 1)   - proportion for independent samples
#   [1, n]   - expected value for independent samples
#   (-1, 0)  - negative proportion for sample of fixed size
#   [-n, -1] - negative size of fixed-size sample
#   c(mean, var) - mean and variance in generalized poisson
#   list(p, probs); p can be one of the above
sampleIndices <- function(n, p) {
  if (!is.numeric(n)) stop("n must be numeric")
  if (is.list(p)) {
    prob <- p[[2L]]
    p <- p[[1L]]
  } else {
    prob <- rep_len(1 / n, n)
  }
  if (!is.numeric(p)) stop("p must be numeric")
  if (length(p) != 1L && length(p) != 2L) stop("p must be of length 1 or 2")
  if (!is.numeric(prob)) stop("prob must be numeric")
  if (length(prob) != n) stop("length of probabilities must equal n")
  
  if (length(p) == 1L) {
    if (p >= 0) {
      prob <- min(if (p >= 1) p / n else p, 1)
      return(which(rbinom(n, 1L, prob) == 1))
    }
    
    if (p > -1) {
      nt <- -p * n
      nt <- floor(nt) + if (runif(1L) < nt - floor(nt)) 1L else 0L
    } else {
      nt <- round(-p)
    }
    nt <- min(nt, n)
    return(sample(n, nt, prob = prob))
  }
  
  if (length(p) == 2L) {
    nt <- min(rgpois(1L, p[1L], p[2L]), n)
    return(sample(n, nt, prob = prob))
  }
  
  which(rbinom(n, 1L, prob) == 1)
}

## as above, but only samples the number of terms not their actual indices
sampleSize <- function(n, p) {
  if (!is.numeric(n)) stop("n must be numeric")
  if (is.list(p)) {
    prob <- p[[2L]]
    p <- p[[1L]]
  } else {
    prob <- rep_len(1 / n, n)
  }
  if (!is.numeric(p)) stop("p must be numeric")
  if (length(p) != 1L && length(p) != 2L) stop("p must be of length 1 or 2")
  if (!is.numeric(prob)) stop("prob must be numeric")
  if (length(prob) != n) stop("length of probabilities must equal n")
  
  if (length(p) == 1L) {
    if (p >= 0) {
      prob <- if (p >= 1) p / n else p
      return(sum(rbinom(n, 1L, prob)))
    }
    
    if (p > -1) {
      nt <- -p * n
      nt <- floor(nt) + if (runif(1L) < nt - floor(nt)) 1L else 0L
    } else {
      nt <- round(-p)
    }
    nt <- min(nt, n)
    return(nt)
  }
  
  if (length(p) == 2L) {
    nt <- min(rgpois(1L, p[1L], p[2L]), n)
    return(nt)
  }
  
  sum(rbinom(n, 1L, prob))
}
