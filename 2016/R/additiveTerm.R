generateUnivariateTerm <- function(C, name, col, functionDist)
{
  if (is.discrete(col))
    return(new("UnivariateAdditiveTerm", name = name, functions = list(newDiscreteBaseFunction(C, col))))
  
  functionSelector <- rbinom(1L, 1L, functionDist[[1L]])
  functionIndex <- functionSelector
  
  if (length(functionDist) > 1L) for (i in seq.int(2L, length(functionDist))) {
    functionProbs <- as.numeric(functionDist[[i]])
    functionSelector <- c(functionSelector, rbinom(1L, 1L, functionProbs[functionIndex + 1L]))
    
    functionIndex <- 2L * functionIndex + functionSelector[i]
  }
  
  functionNames <- names(functionDist)[ifelse(functionSelector == 1, TRUE, FALSE)]
  functions <- lapply(functionNames, function(fnName) newContinuousBaseFunction(C, col, fnName))

  
  new("UnivariateAdditiveTerm", name = name, functions = functions)
}

generateInteractionTerm <- function(C, terms, cols, interactionRetentionProbs)
{
  t1 <- subsampleTermFunctions(C, terms[[1L]], cols[[1L]], interactionRetentionProbs)
  if (length(terms) == 2L) {
    t2 <- subsampleTermFunctions(C, terms[[2L]], cols[[2L]], interactionRetentionProbs)
  } else {
    t2 <- generateInteractionTerm(C, tail(terms, -1L), cols[seq.int(2L, length(cols))], interactionRetentionProbs)
  }
  
  if (is.null(t1) || is.null(t2)) return(NULL)
  
  new("InteractionAdditiveTerm", term1 = t1, term2 = t2)
}

copyAndTweakTerm <- function(C, term, df)
{
  if (is(term, "InteractionAdditiveTerm")) {
    term@term1 <- copyAndTweakTerm(C, term@term1, df)
    term@term2 <- copyAndTweakTerm(C, term@term2, df)
    return(term)
  }
  col <- if (is.null(df) || term@name == "constant" || term@name == "" || !is.data.frame(df)) df else df[[term@name]]
  
  for (i in seq_along(term@functions))
    term@functions[[i]] <- tweakFunction(C, term@functions[[i]])
  
  term
}

subsampleTermFunctions <- function(C, term, df, interactionRetentionProbs)
{
  if (is(term, "InteractionAdditiveTerm")) {
    term@term1 <- subsampleTermFunctions(C, term@term1, df, interactionRetentionProbs)
    term@term2 <- subsampleTermFunctions(C, term@term2, df, interactionRetentionProbs)
    return(term)
  }
  col <- if (is.null(df) || term@name == "constant" || term@name == "" || !is.data.frame(df)) df else df[[term@name]]
  
  functions <- list()
  for (i in seq_along(term@functions)) {
    func.i <- term@functions[[i]]
    keep.prob <- interactionRetentionProbs[[baseFunctions[func.i@type]]]
    
    if (keep.prob != 1 && rbinom(1L, 1L, keep.prob) == 0L) next
    
    newFunction <-
      if (func.i@type == baseFunctionsMap[["step.discrete"]]) {
        func <- newDiscreteBaseFunction(C, col)
        # func <- shuffleBaseline(func)
      } else {
        newContinuousBaseFunction(C, col, baseFunctions[[func.i@type]])
      }
    
    functions <- append(functions, newFunction)
  }
  
  if (length(functions) > 0L)
    new("UnivariateAdditiveTerm", name = term@name, functions = functions)
  else
    NULL
}

interactTermWithTreatment <- function(C, term, df, z, interactionRetentionProbs)
{
  t1 <- subsampleTermFunctions(C, term, df, interactionRetentionProbs)
  
  if (is.null(t1)) return(NULL)
  
  t2 <- generateUnivariateTerm(C, ".z", z, NULL)
  # t2@functions[[1L]] <- shuffleBaseline(t2@functions[[1L]])
  
  new("InteractionAdditiveTerm", term1 = t1, term2 = t2)
}

isTreatmentInteractionTerm <- function(term)
  is(term, "InteractionAdditiveTerm") && is(term@term2, "UnivariateAdditiveTerm") &&
  length(term@term2@functions) == 1L && term@term2@name == ".z"


getTermSize <- function(term)
{
  if (is.null(term)) return(0L)
  if (is(term, "UnivariateAdditiveTerm")) return(1L)
  getTermSize(term@term1) + getTermSize(term@term2)
}

getTermNames <- function(term)
{
  if (is(term, "InteractionAdditiveTerm")) {
    n1 <- getTermNames(term@term1)
    n2 <- getTermNames(term@term2)
    return(c(n1, n2))
  } else {
    if (term@name != "") return(term@name)
  }
  character(0L)
}

i.countEffectiveNumTerms <- function(term)
{
  if (is(term, "UnivariateAdditiveTerm")) {
    if (term@functions[[1L]]@type == baseFunctionsMap[["step.discrete"]]) return(-1L)
    else return(length(term@functions))
  }
  
  c(i.countEffectiveNumTerms(term@term1), i.countEffectiveNumTerms(term@term2))
}

countEffectiveNumTerms <- function(term)
{
  if (is(term, "UnivariateAdditiveTerm"))
    return(if (term@functions[[1L]]@type == baseFunctionsMap[["step.discrete"]]) 0.5 else length(term@functions))
  
  numTerms <- i.countEffectiveNumTerms(term)
  
  discreteTerms <- numTerms[numTerms < 0]
  continuousTerms <- numTerms[numTerms > 0]
  
  result <- 1
  if (length(continuousTerms) > 0L) result <- result / (10^(length(continuousTerms) - 1))
  if (length(discreteTerms) > 0L) result <- result / 2^length(discreteTerms)
  #numerator <- if (length(continuousTerms) > 0) sqrt(prod(continuousTerms)) / (10^(length(continuousTerms) - 1)) else 1
  #numerator <- if (length(continuousTerms) > 0) 1 / (10^(length(continuousTerms) - 1)) else 1
  #denominator <- if (length(discreteTerms) > 0) sqrt(length(discreteTerms)) else 1
  #denominator <- if (length(discreteTerms) > 0) 1.05 * length(discreteTerms) else 1
  
  #numerator / denominator
  result
}

getTermWeight <- function(term)
{
  if (is(term, "UnivariateAdditiveTerm"))
    return(if (term@functions[[1]]@type == baseFunctionsMap[["step.discrete"]]) 0.5 else 1)
  
  numTerms <- i.countEffectiveNumTerms(term)
  
  discreteTerms <- numTerms[numTerms < 0]
  continuousTerms <- numTerms[numTerms > 0]
  
  result <- 1
  # if (length(continuousTerms) > 0) result <- result / (10^(length(continuousTerms) - 1))
  if (length(discreteTerms) > 0L) result <- result / 2^length(discreteTerms)
  
  result
}

allTermsAreDiscrete <- function(term)
{
  if (is(term, "UnivariateAdditiveTerm"))
    return(term@functions[[1L]]@type == baseFunctionsMap[["step.discrete"]])
  
  allTermsAreDiscrete(term@term1) && allTermsAreDiscrete(term@term2)
}

unsetAllDiscreteCoefficients <- function(term)
{
  if (is(term, "InteractionTerm")) {
    term@term1 <- unsetAllDiscreteCoefficients(term@term1)
    term@term2 <- unsetAllDiscreteCoefficients(term@term2)
  } else {
    term@functions[[1L]]@pars[term@functions[[1L]]@pars != 0] <- 1
  }
  term
}

unsetAllButLeftmostDiscreteCCoefficient <- function(term)
{
  if (is(term, "InteractionTerm")) {
    term@term1 <- unsetAllButLeftmostDiscreteCCoefficient(term@term1)
    term@term2 <- unsetAllDiscreteCoefficients(term@term2)
  }
  term
}

rescaleInteractionTerm <- function(term, numTotalTerms)
{
  # if (allTermsAreDiscrete(term))
  
  if (is(term, "InteractionAdditiveTerm")) {
    term@term1 <- rescaleInteractionTerm(term@term1, numTotalTerms)
    term@term2 <- rescaleInteractionTerm(term@term2, numTotalTerms)
  } else {
    for (i in seq_along(term@functions)) {
      term@functions[[i]] <- rescaleBaseFunction(term@functions[[i]], (1 / numTotalTerms))
    }
  }
  term
}

absInteractionTerm <- function(term)
{
  if (is(term, "InteractionAdditiveTerm")) {
    term@term1 <- absInteractionTerm(term@term1)
    term@term2 <- absInteractionTerm(term@term2)
  } else {
    for (i in seq_along(term@functions)) {
      term@functions[[i]] <- absBaseFunction(term@functions[[i]])
    }
  }
  term
}

evaluateAdditiveTermToDataframe <- function(term, x)
{
  if (is(term, "InteractionAdditiveTerm")) {
    df1 <- evaluateAdditiveTermToDataframe(term@term1, x)
    df2 <- evaluateAdditiveTermToDataframe(term@term2, x)
    
    result <- vector("list", length(df1) * length(df2))
    names <- character(length(df1) * length(df2))
    iter <- 1L
    for (i in seq_along(df1)) {
      for (j in seq_along(df2)) {
        result[[iter]] <- df1[[i]] * df2[[j]]
        names[iter] <- paste0(names(df1)[i], " * ", names(df2)[j])
        iter <- iter + 1L
      }
    }
    names(result) <- names
  } else {
    if (term@name == "constant")
      return({ result <- list(rep_len(term@functions[[1L]]@f(0), NROW(x))); names(result) <- format(term); result })
  
    col <- x[[term@name]]
    if (is.discrete(col) && length(term@functions[[1L]]@pars) > 2L) {
      pars <- term@functions[[1L]]@pars
      isZero <- pars == 0
      parMap <- attr(pars, "map")
      colIndices <- match(col, parMap)
      parMatr <- diag(pars)[,!isZero]
      result <- parMatr[colIndices,]
      colnames(result) <- paste0("I(", term@name, " == ", parMap[!isZero], ")")
      result <- as.data.frame(result)
    } else {
      result <- vector("list", length(term@functions))
      for (i in seq_along(term@functions)) {
        result[[i]] <- term@functions[[i]]@f(col)
      }
      names(result) <- sapply(term@functions, function(func.i) format(func.i, term@name))
    }
  }
  result
}

