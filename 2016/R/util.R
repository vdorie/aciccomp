"%not_in%" <- function(x, table) match(x, table, nomatch = 0L) <= 0L

is.discrete <- function(col) is.factor(col) || length(unique(col)) == 2L

## from lme4
namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1L]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

evalx.recurse <- function(x, e) {
  if (length(e) == 0L || typeof(e) == "symbol") return(e)
  
  for (i in seq_along(e)) {
    if (!is.language(e[[i]])) next
    
    e[[i]] <- if (e[[i]] == "x") x else evalx.recurse(x, e[[i]])
  }
  
  e
}

## evaluates the expression 'e' by after first replacing all instances of 'x' with the expression x
x <- NULL
evalx <- function(x, e) {
  mc <- match.call()
  callingEnv <- parent.frame()
  
  e <- evalx.recurse(mc$x, mc$e)
  eval(e, callingEnv)
}

