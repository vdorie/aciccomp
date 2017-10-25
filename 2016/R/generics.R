setGeneric("evaluate", function(object, x) standardGeneric("evaluate"))

setMethod("evaluate", "UnivariateAdditiveTerm", function(object, x) {
  if (object@name == "constant") return(rep_len(object@functions[[1L]]@f(0), NROW(x)))
  
  col <- x[[object@name]]
  result <- rep(0, length(col))
  for (i in seq_along(object@functions)) {
    result <- result + object@functions[[i]]@f(col)
  }
  result
})

setMethod("evaluate", "InteractionAdditiveTerm", function(object, x) {
  evaluate(object@term1, x) * evaluate(object@term2, x)
})

setMethod("evaluate", "AdditiveFunction", function(object, x) {
  if (!is.data.frame(x) && is.null(names(x))) stop("argument to evaluate must be a data frame")
  
  for (i in seq_along(x))
    if (!is.discrete(x[[i]])) x[[i]] <- object@inputScale * x[[i]]
  
  
  if (object@operator == operatorsMap[["+"]]) {
    result <- rep(0, NROW(x))
    for (i in seq_along(object@terms)) {
      res.i <- evaluate(object@terms[[i]], x)
      if (anyNA(res.i)) break
      result <- result + object@weights[i] * res.i #evaluate(object@terms[[i]], x)
    }
    if (length(object@terms) == 0L || i != length(object@terms)) browser()
  } else {
    result <- rep(1, NROW(x))
    if (length(object@terms) == 0L)
      result <- 0
    else for (i in seq_along(object@terms)) {
      result <- result * object@weights[i] * evaluate(object@terms[[i]], x)
    }
  }
  applyTransformation(object@transformation, object@outputScale * result)
})
setMethod("evaluate", "InteractionFunction", function(object, x) {
  result <- evaluate(object@func1, x) * evaluate(object@func2, x)
  applyTransformation(object@transformation, result)
})

setMethod("evaluate", "GeneralizedAdditiveFunction", function(object, x) {
  if (!is.data.frame(x) && is.null(names(x))) stop("argument to evaluate must be a data frame")
  
  result <- rep(0, NROW(x))
  for (i in seq_along(object@functions)) {
    f.x <- evaluate(object@functions[[i]], x)
    w <- if (is.function(object@weights[[i]])) object@weights[[i]](x) else object@weights[[i]]
    result <- result + w * f.x
  }
  result
})

format.BaseFunction <- function(x, termName, symbolic = TRUE, digits = getOption("digits"), ...) {
  par1 <- signif(x@pars[1L], digits)
  switch(baseFunctions[x@type],
    constant  = if (symbolic) "C" else par1,
    linear    = if (symbolic) termName else paste0(par1, " * ", termName),
    quadratic = if (symbolic) paste0(termName, "^2") else paste0(par1, " * ", termName, "^2"),
    cubic     = if (symbolic) paste0(termName, "^3") else paste0(par1, " * ", termName, "^3"),
    step.constant = {
      if (symbolic)
        paste0("I(", termName, " <= b)")
      else
        paste0(signif(x@pars[2L], digits), " * I(", termName, if (x@pars[1L] < 0)  " <= " else " > ", par1, ")")
    },
    step.linear = {
      if (symbolic)
        paste0(termName, " * I(", termName, " <= b)")
      else
        paste0(signif(x@pars[2L], digits), " * ", termName, " * I(", termName, if (x@pars[1L] < 0) " <= " else " > ", par1, ")")
    },
    step.discrete = {
      if (symbolic) {
        paste0("I(", termName, " == ", if (length(x@pars) == 2L) attr(x@pars, "map")[x@pars != 0] else "b", ")")
      } else {
        pars <- as.numeric(x@pars[x@pars != 0])
        parNames <- attr(x@pars, "map")[x@pars != 0]
        paste0(signif(pars, digits), " * I(", termName, " == ", parNames, ")", collapse = " + ")
      }
    },
    sigmoid = {
      if (symbolic)
        paste0("pnorm(a * (", termName, " - b))")
      else
        paste0("pnorm(", if (x@pars[2L] < 0) "-" else "", par1, " * (", termName, " - ", signif(x@pars[2L], digits), "))")
    },
    step.quantile = {
      if (symbolic)
        paste0("I(", termName, " <= q)")
      else
        paste0("I(", termName, if (x@pars[1L] > 0) " >= " else " <= ", signif(x@pars[2L], digits), ")")
    }
  )
}

format.UnivariateAdditiveTerm <- function(x, ...) {
  if (length(x@functions) == 0L) return("0")
  
  paste0(sapply(x@functions, function(f.i) format(f.i, termName = x@name, ...)), collapse = " + ")
}
format.InteractionAdditiveTerm <- function(x, ...) {
  t1String <- format(x@term1, ...)
  t2String <- format(x@term2, ...)
  
  
  if (is(x@term1, "InteractionAdditiveTerm") || (is(x@term1, "UnivariateAdditiveTerm") && length(x@term1@functions) > 1L))
    t1String <- paste0("(", t1String, ")")
  if (is(x@term2, "InteractionAdditiveTerm") || (is(x@term2, "UnivariateAdditiveTerm") && length(x@term2@functions) > 1L))
    t2String <- paste0("(", t2String, ")")
  
  paste0(t1String, " * ", t2String)
}


format.AdditiveFunction <- function(x, symbolic = TRUE, digits = getOption("digits"), ...) {
  if (length(x@terms) == 0L) return("0")
  
  operator <- if (x@operator == operatorsMap[["+"]]) " + " else " * "
  
  pre <- getPrefixForTransformation(x@transformation, symbolic, digits)
  suf <- getSuffixForTransformation(x@transformation, symbolic, digits)

  paste0(pre, paste0(sapply(x@terms, function(term.i) format(term.i, symbolic, digits, ...)), collapse = operator), suf)
}

format.InteractionFunction <- function(x, symbolic = TRUE, digits = getOption("digits"), ...)
{
  f1String <- format(x@func1, symbolic, digits, ...)
  f2String <- format(x@func2, symbolic, digits, ...)
  
  pre <- getPrefixForTransformation(x@transformation, symbolic, digits)
  suf <- getSuffixForTransformation(x@transformation, symbolic, digits)
  
  if (is(x@func1, "AdditiveFunction") && x@func1@transformation@type == transformationsMap[["identity"]])
    f1String <- paste0("(", f1String, ")")
  if (is(x@func2, "AdditiveFunction") && x@func2@transformation@type == transformationsMap[["identity"]])
    f2String <- paste0("(", f2String, ")")
  
  paste0(pre, f1String, " * ", f2String, suf)
}

format.GeneralizedAdditiveFunction <- function(x, symbolic = TRUE, digits = getOption("digits"), ...)
{
  if (length(x@functions) == 0L) return("0")
  
  pre <- if (symbolic) "" else paste0(signif(x@weights[1L], digits), " * ")
  result <- paste0(pre, format(x@functions[[1L]], symbolic, digits, ...))
  
  if (length(x@functions) > 1L) for (i in seq.int(2L, length(x@functions))) {
    pre <- if (symbolic) "" else paste0(signif(x@weights[i], digits), " * ")
    
    result <- paste0(result, " + ", pre, format(x@functions[[i]], symbolic, digits, ...))
  }
  
  result
}

print.UnivariateAdditiveTerm <- function(x, ...) cat(format(x, ...), "\n", sep = "")
print.InteractionAdditiveTerm <- function(x, ...) cat(format(x, ...), "\n", sep = "")
print.AdditiveFunction <- function(x, ...) cat(format(x, ...), "\n", sep = "")
print.InteractionFunction <- function(x, ...) cat(format(x, ...), "\n", sep = "")
print.GeneralizedAdditiveFunction <- function(x, ...) cat(format(x, ...), "\n", sep = "")

setMethod("show", "UnivariateAdditiveTerm", function(object) print(object))
setMethod("show", "InteractionAdditiveTerm", function(object) print(object))
setMethod("show", "AdditiveFunction", function(object) print(object))
setMethod("show", "InteractionFunction", function(object) print(object))
setMethod("show", "GeneralizedAdditiveFunction", function(object) print(object))

