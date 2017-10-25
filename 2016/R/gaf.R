if (FALSE) addConstantTermToAdditiveFunction <- function(C, af, mean, scale)
{
  constantTerm <- new("UnivariateAdditiveTerm", name = "constant",
                      functions = list(newContinuousBaseFunction("constant", mean + scale * rt(1L, C$BF_DF))))
  af@terms <- append(constantTerm, af@terms)
  af@weights <- c(1, af@weights)
  af
}

setBaselineForAdditiveFunction <- function(C, af, x, m) {
  f.x <- evaluate(af, x)
  constantTerm <- new("UnivariateAdditiveTerm", name = "constant", functions =
                      list(newContinuousBaseFunction(C, NULL, "constant", m - mean(f.x))))
  af@terms <- append(constantTerm, af@terms)
  af@weights <- c(1, af@weights)
  af
}

setBaselineForGeneralizedAdditiveFunction <- function(C, gaf, x, m)
{
  f.x <- evaluate(gaf, x)
  constantTerm <- new("UnivariateAdditiveTerm", name = "constant", functions =
                      list(newContinuousBaseFunction(C, NULL, "constant", m - mean(f.x))))
  constantFunc <- new("AdditiveFunction", transformation = newTransformation("identity"),
                      inputScale = 1, outputScale = 1,
                      terms = list(constantTerm), weights = 1, operator = operatorsMap[["+"]])
  gaf@functions <- append(constantFunc, gaf@functions)
  gaf@weights   <- append(1, gaf@weights)
  gaf
}

addConstantTermToGeneralizedAdditiveFunction <- function(gaf, const)
{
  constantTerm <- new("UnivariateAdditiveTerm", name = "constant", functions =
                      list(newContinuousBaseFunction(NULL, "constant", const)))
  constantFunc <- new("AdditiveFunction", transformation = newTransformation("identity"),
                      inputScale = 1, outputScale = 1,
                      terms = list(constantTerm), weights = 1, operator = operatorsMap[["+"]])
  gaf@functions <- append(constantFunc, gaf@functions)
  gaf@weights   <- append(1, gaf@weights)
  gaf
}

if (FALSE) newScaleFunction <- function(term, name, values)
{
  scaleTerm <- new("UnivariateAdditiveTerm", name = name, functions =
                   list(newDiscreteBaseFunction(term, values)))
  scaleFunc <- new("AdditiveFunction", transformation = newTransformation("identity"),
                   inputScale = 1, outputScale = 1,
                   terms = list(scaleTerm), weights = 1, operator = operatorsMap[["+"]])
  scaleFunc@isScale <- TRUE
  scaleFunc
}


isConstantAdditiveFunction <- function(af)
  is(af, "AdditiveFunction") && length(af@terms) == 1L &&
  is(af@terms[[1L]], "UnivariateAdditiveTerm") && af@terms[[1L]]@name == "constant"

combineGafs <- function(gaf1, gaf2)
{
  gaf1@functions <- append(gaf1@functions, gaf2@functions)
  gaf1@weights <- append(gaf1@weights, gaf2@weights)
  #numNonConstant <- length(gaf1@functions)
  #for (i in seq_along(gaf1@functions)) {
  #  if (isConstantAdditiveFunction(gaf1@functions[[i]]))
  #    numNonConstant <- numNonConstant - 1L
  #}
  #for (i in seq_along(gaf1@functions)) {
  #  if (!isConstantAdditiveFunction(gaf1@functions[[i]]))
  #    gaf1@weights[i] <- gaf1@weights[i] / sqrt(numNonConstant)
  #}
  gaf1
}

generateAdditiveFunction <- function(C, df, scales, dist)
{
  p <- ncol(df)
  
  termInclusionProb <- if (is.function(dist@termInclusionProb)) dist@termInclusionProb(p) else dist@termInclusionProb
  
  includedTerms <- sampleIndices(p, termInclusionProb)
  
  numTerms <- length(includedTerms)
  termNames <- colnames(df)[includedTerms]
  
  terms <- vector("list", numTerms)
  for (i in seq_len(numTerms))
    terms[[i]] <- generateUnivariateTerm(C, termNames[i], df[[termNames[i]]], dist@baseFunctionDist)
  
  interactionInclusionProbs <- dist@interactionInclusionProbs
  interactionRetentionProbs <- dist@interactionRetentionProbs
  
  interactionTerms <- list()
  if (!is.null(interactionInclusionProbs) && length(interactionInclusionProbs) > 0) for (i in seq_along(interactionInclusionProbs))
  {
    if (numTerms < i + 1L) break
    marginalInteractionProb <- if (is.list(interactionInclusionProbs)) interactionInclusionProbs[[i]] else interactionInclusionProbs[i]
    
    numPossibleInteractions <- choose(numTerms, i + 1L)
    interactionIndices <- sampleIndices(numPossibleInteractions, marginalInteractionProb)
    
    interactions <- combn(numTerms, i + 1L, simplify = TRUE)[,interactionIndices, drop = FALSE]
    
    for (j in seq_len(ncol(interactions))) {
      terms.j <- terms[interactions[,j]]
      cols.j <- df[termNames[interactions[,j]]]
      interactionTerm <- generateInteractionTerm(C, terms.j, cols.j, interactionRetentionProbs)
      if (!is.null(interactionTerm)) {
        interactionTerm <- rescaleInteractionTerm(interactionTerm, i + 1L)
        interactionTerms <- append(interactionTerms, interactionTerm)
      }
    }
  }
  terms <- append(terms, interactionTerms)
  
  #weights <- rep_len(0, length(terms))
  #for (i in seq_along(terms))
  #  weights[i] <- getTermWeight(terms[[i]])
  #weights <- sqrt(weights / sum(weights))
  #weights <- replicate(length(terms), { scale <- rgamma(1L, 4 * 2, 2); rgamma(1L, 1 * scale, scale) })
  #weights <- (weights / sum(weights))^(1/3)
  #effectiveNumTerms <- 0
  #for (i in seq_along(terms))
  #  effectiveNumTerms <- effectiveNumTerms + countEffectiveNumTerms(terms[[i]])
  
  #weights <- rep_len(1 / sqrt(effectiveNumTerms), length(terms))
  # weights <- rt(length(terms), 1) / 4.5
  
  weights <- getWeightsForTerms(terms)
    
  new("AdditiveFunction", inputScale = scales[1], outputScale = scales[2],
               terms = terms, weights = weights, operator = operatorsMap[["+"]])
}

generateAdditiveFunctionFromPrototype <- function(C, df, scales, dist, af, alignmentProb)
{
  p <- ncol(df)
  
  termInclusionProb <- if (is.function(dist@termInclusionProb)) dist@termInclusionProb(p) else dist@termInclusionProb
  
  termNames.af <- lapply(af@terms, getTermNames)
  termLengths.af <- sapply(termNames.af, length)
  
  indices.af <- seq_len(length(termNames.af))[termLengths.af == 1L]
  termIndicesToKeep <- indices.af[runif(length(indices.af)) < alignmentProb]
  numTermsKept <- length(termIndicesToKeep)
  
  terms <- vector("list", numTermsKept)
  
  for (i in seq_len(numTermsKept))
    terms[[i]] <- copyAndTweakTerm(C, af@terms[[termIndicesToKeep[i]]], df)
  
  termsKept <- unlist(termNames.af[termIndicesToKeep])
  isKeptTerm <- colnames(df) %in% termsKept
  
  
  totalNumTerms <- sampleSize(p, termInclusionProb)
  if (totalNumTerms > numTermsKept) {
    if (is.list(termInclusionProb)) {
      termProb <- termInclusionProb[[2L]][!isKeptTerm]
      termInclusionProb <- list(-(totalNumTerms - numTermsKept), termProb)
    } else {
      termInclusionProb <- -(totalNumTerms - numTermsKept)
    }
    includedTerms <- sampleIndices(p - numTermsKept, termInclusionProb)
    termNames <- colnames(df)[!isKeptTerm][includedTerms]
    for (i in seq_along(includedTerms))
      terms <- append(terms, generateUnivariateTerm(C, termNames[i], df[[termNames[i]]], dist@baseFunctionDist))
  }
  
  numTerms <- length(terms)    
  termNames <- sapply(terms, getTermNames)
  
  interactionInclusionProbs <- dist@interactionInclusionProbs
  interactionRetentionProbs <- dist@interactionRetentionProbs
  
  interactionTerms <- list()
  if (!is.null(interactionInclusionProbs) && length(interactionInclusionProbs) > 0) for (i in seq_along(interactionInclusionProbs))
  {
    if (numTerms < i + 1L) break
    marginalInteractionProb <- if (is.list(interactionInclusionProbs)) interactionInclusionProbs[[i]] else interactionInclusionProbs[i]
    
    interactionIndices.af <- seq_len(length(termNames.af))[termLengths.af == i + 1L]
    interactionNames.af <- lapply(interactionIndices.af, function(index) sort(termNames.af[[index]]))
    
    interactionIndicesToKeep <- interactionIndices.af[runif(length(interactionIndices.af)) < alignmentProb]
    numInteractionsKept <- length(interactionIndicesToKeep)
    
    for (j in seq_along(interactionIndicesToKeep))
      terms <- append(terms, copyAndTweakTerm(C, af@terms[[interactionIndicesToKeep[j]]], df))
    
    numPossibleInteractions <- choose(numTerms, i + 1L)
    totalNumInteractions <- sampleSize(numPossibleInteractions, marginalInteractionProb)
    
    if (totalNumInteractions > numInteractionsKept) {
      allInteractions <- combn(numTerms, i + 1L, simplify = TRUE)
      allInteractionsString <-
        sapply(seq_len(ncol(allInteractions)), function(col) paste0(allInteractions[,col], collapse = ","))
      currentInteractionsString <-
        sapply(interactionNames.af, function(interaction) paste0(sort(match(interaction, termNames)), collapse = ","))
      remainingInteractions <- setdiff(allInteractionsString, currentInteractionsString)
      
      if (length(remainingInteractions) == 0L) next
      
      interactions <- remainingInteractions[sampleIndices(length(remainingInteractions), -(totalNumInteractions - numInteractionsKept))]
      for (j in seq_along(interactions)) {
        interactionIndices <- as.integer(strsplit(interactions[j], ",")[[1L]])
        interactionNames <- termNames[interactionIndices]
        
        interactionTerms <- terms[interactionIndices]
        interactionTerm <- generateInteractionTerm(C, interactionTerms, df[interactionNames], interactionRetentionProbs)
        if (!is.null(interactionTerm)) {
          interactionTerm <- rescaleInteractionTerm(interactionTerm, i + 1L)
          terms <- append(terms, interactionTerm)
        }
      }
    }
  }
  
  #weights <- replicate(length(terms), { scale <- rgamma(1L, 4 * 2, 2); rgamma(1L, 1 * scale, scale) })
  #weights <- (weights / sum(weights))^(1/3)
  #weights <- rep_len(0, length(terms))
  #for (i in seq_along(terms))
  #  weights[i] <- getTermWeight(terms[[i]])
  #weights <- sqrt(weights / sum(weights))
  
  #effectiveNumTerms <- 0
  #for (i in seq_along(terms))
  #  effectiveNumTerms <- effectiveNumTerms + countEffectiveNumTerms(terms[[i]])
  #
  #weights <- rep_len(1 / sqrt(effectiveNumTerms), length(terms))
  #weights <- rt(length(terms), 1) / 4.5
  weights <- getWeightsForTerms(terms)
  
  new("AdditiveFunction", inputScale = scales[1], outputScale = scales[2],
               terms = terms, weights = weights, operator = operatorsMap[["+"]])
}

generateGeneralizedAdditiveFunction <- function(C, df, gafTransformations, afProbs, afScales, afDists)
{
  if (!is.matrix(afProbs))  afProbs <- matrix(afProbs, 1)
  if (!is.matrix(afScales)) afScales <- matrix(afScales, 1)
  if (!is.list(afDists))    afDists <- list(afDists)
  
  afs <- list()
  
  for (i in seq_along(gafTransformations)) {
    numFunctions <- rgpois(1L, afProbs[i,1L], afProbs[i,2L])
    afs.i <- list()
    
    if (numFunctions > 0L) for (j in seq_len(numFunctions)) {
      af <- generateAdditiveFunction(C, df, afScales[i,], afDists[[i]])
      if (is.null(af)) next
      
      af@transformation <- newTransformation(gafTransformations[[i]])
      
      afs.i[[length(afs.i) + 1L]] <- af
    }
    afs <- append(afs, afs.i)
  }
  
  new("GeneralizedAdditiveFunction", functions = afs, weights = as.list(rep(1, length(afs))))
}

generateGeneralizedAdditiveFunctionFromPrototype <-
  function(C, df, gafTransformations, afProbs, afScales, afDists, gaf, alignmentProb)
{
  if (!is.matrix(afProbs))  afProbs <- matrix(afProbs, 1)
  if (!is.matrix(afScales)) afScales <- matrix(afScales, 1)
  if (!is.list(afDists))    afDists <- list(afDists)
  
  ## find principal linear additive function to get names of main effects
  af.main <- NULL
  for (i in seq_along(gaf@functions)) {
    if (is(gaf@functions[[i]], "AdditiveFunction") && length(gaf@functions[[i]]@terms) > 1) {
      af.main <- gaf@functions[[i]]
      break
    }
  }
  me.names <- sapply(af.main@terms, function(term.i) if (is(term.i, "UnivariateAdditiveTerm")) term.i@name else "")
  
  ## go through any biasing functions and add linear terms to af corresponding to
  ## the quantile cutoffs
  for (i in seq_along(gaf@functions)) {
    func.i <- gaf@functions[[i]]
    if (is(func.i, "AdditiveFunction") && func.i@operator == operatorsMap[["*"]]) {
      for (j in seq_along(func.i@terms)) {
        term.j <- func.i@terms[[j]]
        if (term.j@name %not_in% me.names) {
          af.main@terms[[length(af.main@terms) + 1L]] <-
            new("UnivariateAdditiveTerm", name = term.j@name, functions = list(newContinuousBaseFunction(C, NULL, "linear")))
          af.main@weights <- c(af.main@weights, af.main@weights[1L])
          me.names <- c(me.names, term.j@name)
        }
      }
    }
  }
  
  afs <- list()
  
  for (i in seq_along(gafTransformations)) {
    numFunctions <- rgpois(1L, afProbs[i,1L], afProbs[i,2L])
    afs.i <- list()
    
    if (numFunctions > 0L) for (j in seq_len(numFunctions)) {
      af <- generateAdditiveFunctionFromPrototype(C, df, afScales[i,], afDists[[i]], af.main, alignmentProb)
      if (is.null(af)) next
      
      af@transformation <- newTransformation(gafTransformations[[i]])
      
      afs.i[[length(afs.i) + 1L]] <- af
    }
    afs <- append(afs, afs.i)
  }
  
  new("GeneralizedAdditiveFunction", functions = afs, weights = as.list(rep(1, length(afs))))
}

generateGeneralizedAdditiveFunctionFromPrototypeAndForceOverlapAlignment <-
  function(C, df, gafTransformations, afProbs, afScales, afDists, gaf, alignmentProb)
{
  if (!is.matrix(afProbs))  afProbs <- matrix(afProbs, 1)
  if (!is.matrix(afScales)) afScales <- matrix(afScales, 1)
  if (!is.list(afDists))    afDists <- list(afDists)
  
  ## find principal linear additive function to duplicate effects
  af.main <- NULL
  for (i in seq_along(gaf@functions)) {
    func.i <- gaf@functions[[i]]
    if (is(func.i, "AdditiveFunction") && func.i@operator == operatorsMap[["+"]] && length(func.i@terms) > 0L &&
        length(func.i@terms[[1L]]@functions) > 0L && func.i@terms[[1L]]@functions[[1L]]@type != baseFunctionsMap[["constant"]]) {
      af.main <- func.i
      break
    }
  }
  if (is.null(af.main)) stop("cannot find principal additive function")
  
  ## see if there is a biasing function, if so add and interact the terms
  af.bias <- NULL
  for (i in seq_along(gaf@functions)) {
    func.i <- gaf@functions[[i]]
    if (is(func.i, "AdditiveFunction") && func.i@operator == operatorsMap[["*"]]) {
      af.bias <- func.i
      break
    }
  }
  
  addedBiasTerms <- is.null(af.bias)
   
  afs <- list()
  
  for (i in seq_along(gafTransformations)) {
    numFunctions <- rgpois(1L, afProbs[i,1L], afProbs[i,2L])
    afs.i <- list()
    
    if (numFunctions > 0L) for (j in seq_len(numFunctions)) {
      af <- generateAdditiveFunctionFromPrototype(C, df, afScales[i,], afDists[[i]], af.main, alignmentProb)
      if (is.null(af)) next
      
      af@transformation <- newTransformation(gafTransformations[[i]])
      
      if (!addedBiasTerms && is(af, "AdditiveFunction") &&
          af@transformation@type == transformationsMap[["identity"]] &&
          af@operator == operatorsMap[["+"]])
      {
        numTerms <- length(af.bias@terms)
        biasTerms <- vector("list", 2L^numTerms - 1L)
        biasIndex <- 1L
        for (k in seq_len(numTerms)) {
          combos <- combn(numTerms, k)
          
          if (nrow(combos) == 1L) {
            for (l in seq_len(ncol(combos))) {
              term.l <- af.bias@terms[[l]]
              biasTerms[[biasIndex]] <- generateUnivariateTerm(C, term.l@name, df[[term.l@name]], afDists[[i]]@baseFunctionDist)
              biasIndex <- biasIndex + 1L
            }
          } else {
            for (l in seq_len(ncol(combos))) {
              ## indexing should work since we will have already generated these
              terms.l <- biasTerms[combos[,l]]
              cols.l  <- df[sapply(terms.l, function(term) term@name)] 
              biasTerms[[biasIndex]] <- generateInteractionTerm(C, terms.l, cols.l, afDists[[i]]@interactionRetentionProbs)
              biasIndex <- biasIndex + 1L
            }
          }
        }
        af@terms <- append(af@terms, biasTerms)
        af@weights <- getWeightsForTerms(af@terms)
      }
      
      afs.i[[length(afs.i) + 1L]] <- af
    }
    afs <- append(afs, afs.i)
  }
  
  new("GeneralizedAdditiveFunction", functions = afs, weights = as.list(rep(1, length(afs))))
}

getWeightsForTerms <- function(terms)
{
  effectiveNumTerms <- 0
  for (i in seq_along(terms))
    effectiveNumTerms <- effectiveNumTerms + countEffectiveNumTerms(terms[[i]])
  
  rep_len(1 / sqrt(effectiveNumTerms), length(terms))
}

createInteractionFunctionFromList <- function(afs, transformation = "power", pars = 1 / length(afs))
{
  if (length(afs) == 2) {
    result <- new("InteractionFunction",
                  transformation = newTransformation(transformation, pars),
                  func1 = afs[[1]],
                  func2 = afs[[2]])
  } else {
    result <- new("InteractionFunction",
                  transformation = newTransformation(transformation, pars),
                  func1 = afs[[1]],
                  func2 = createInteractionFunctionFromList(tail(afs, -1), "identity", numeric()))
  }
  result
}

if (FALSE) generateBiasingFunction <- function(df, afProbs, afScales, afDists)
{
  if (!is.matrix(afProbs)) afProbs <- matrix(afProbs, 1)
  if (!is.matrix(afScales)) afScales <- matrix(afScales, 1)
  if (!is.list(afDists)) afDists <- list(afDists)
  
  afs <- list()
  for (i in seq_along(afDists)) {
    numFunctions <- rgpois(1L, afProbs[i,1L], afProbs[i,2L])
    afs.i <- list()
    
    if (numFunctions > 0L) for (j in seq_len(numFunctions)) {
      af <- generateAdditiveFunction(df, afScales[i,], afDists[[i]])
      if (length(af@terms) == 0L) next
      
      af@transformation <- newTransformation("sigmoid")
      afs.i[[length(afs.i) + 1L]] <- af
    }
    afs <- append(afs, afs.i)
  }
  if (length(afs) > 1L) {
    createInteractionFunctionFromList(afs)
  } else if (length(afs) == 1L) {
    afs[[1L]]
  } else {
    new("AdditiveFunction", transformation = newTransformation("identity"),
        inputScale = 1.0, outputScale = 1.0, terms = list(), weights = numeric())
  }
}

if (FALSE) rescaleAdditiveFunction <- function(af) {
  if (is(af, "AdditiveFunction")) {
    for (i in seq_along(af@terms))
      af@terms[[i]] <- rescaleInteractionTerm(af@terms[[i]], getTermSize(af@terms[[i]]))
  } else {
    af@func1 <- rescaleAdditiveFunction(af@func1)
    af@func2 <- rescaleAdditiveFunction(af@func2)
  }
  af
}

absAdditiveFunction <- function(af) {
  if (is(af, "AdditiveFunction")) {
    for (i in seq_along(af@terms)) {
      af@terms[[i]] <- absInteractionTerm(af@terms[[i]])
    }
  } else {
    af@func1 <- absAdditiveFunction(af@func1)
    af@func2 <- absAdditiveFunction(af@func2)
  }
  af
}

if (FALSE) generateBiasingFunction2 <- function(df, afProbs, afScales, afDists)
{
  if (!is.matrix(afProbs)) afProbs <- matrix(afProbs, 1)
  if (!is.matrix(afScales)) afScales <- matrix(afScales, 1)
  if (!is.list(afDists)) afDists <- list(afDists)
  
  afs <- list()
  for (i in seq_along(afDists)) {
    numFunctions <- rgpois(1L, afProbs[i,1L], afProbs[i,2L])
    afs.i <- list()
    
    if (numFunctions > 0L) for (j in seq_len(numFunctions)) {
      af <- generateAdditiveFunction(df, afScales[i,], afDists[[i]])
      if (length(af@terms) == 0L) next
      
      af <- absAdditiveFunction(af)
      af <- rescaleAdditiveFunction(af)
      #if (is(af, "AdditiveFunction")) {
      #  for (j in seq_along(af@terms)) {
      #    af@terms[[j]] <- absInteractionTerm(af@terms[[j]])
      #    af@terms[[j]] <- rescaleInteractionTerm(af@terms[[j]], getTermSize(af@terms[[j]]))
      #  }
      #} else {
      
      af@transformation <- newTransformation("identity")
      afs.i[[length(afs.i) + 1L]] <- af
    }
    afs <- append(afs, afs.i)
  }
  if (length(afs) > 1L) {
    createInteractionFunctionFromList(afs, "identity", numeric())
  } else if (length(afs) == 1L) {
    afs[[1L]]
  } else {
    new("AdditiveFunction", transformation = newTransformation("identity"),
                 inputScale = 1.0, outputScale = 1.0, terms = list(), weights = numeric())
  }
}

generateBiasingFunction3 <- function(C, df, scales, dist)
{
  p <- ncol(df)
  
  termInclusionProb <- if (is.function(dist@termInclusionProb)) dist@termInclusionProb(p) else dist@termInclusionProb
  
  includedTerms <- sampleIndices(p, termInclusionProb)
  
  numTerms <- length(includedTerms)
  termNames <- colnames(df)[includedTerms]
  
  terms <- vector("list", numTerms)
  for (i in seq_len(numTerms)) {
    terms[[i]] <- generateUnivariateTerm(C, termNames[i], df[[termNames[i]]], dist@baseFunctionDist)
  }
  
  numTerms <- length(terms)
  weights <- rep_len(1, numTerms)
  
  new("AdditiveFunction", inputScale = scales[1], outputScale = scales[2], terms = terms, weights = weights,
               operator = operatorsMap[["*"]])
}

addAdditiveFunctionToGeneralizedAdditiveFunction <- function(gaf, af, weight = 1)
{
  gaf@functions <- append(gaf@functions, af)
  gaf@weights <- append(gaf@weights, weight)
  gaf
}

addTreatmentToGeneralizedAdditiveFunction <- function(C, gaf, df, z, dist)
{
  interactionInclusionsProbs <- dist@interactionInclusionProbs
  interactionRetentionProbs  <- dist@interactionRetentionProbs
  
  for (i in seq_along(gaf@functions)) {
    func.i <- gaf@functions[[i]]
    
    if (!isConstantAdditiveFunction(func.i) && is(func.i, "AdditiveFunction") && 
        func.i@operator == operatorsMap[["+"]])
    {
      if (func.i@transformation@type == transformationsMap[["identity"]]) {
        gaf@functions[[i]] <- addTreatmentToAdditiveFunction(C, func.i, df, z, interactionInclusionsProbs, interactionRetentionProbs)
      } else {
        treatmentTerm <- generateUnivariateTerm(C, ".z", z, c(step.discrete = 1))
        if ("addTreatmentToGAF" %in% C$RUN_BUGGED)
          treatmentTerm@functions[[1L]] <- shuffleBaseline(treatmentTerm@functions[[1L]])
        
        af <-
          new("AdditiveFunction", transformation = newTransformation("identity"), 
              inputScale = 1, outputScale = 1, terms = list(treatmentTerm), weights = 1,
              operator = operatorsMap[["+"]])
          
        gaf@functions[[i]] <-
          new("InteractionFunction", transformation = newTransformation("identity"),
              func1 = gaf@functions[[i]], func2 = af)
      }
    }
  }
  gaf
}

addTreatmentToInteractionFunction <- function(C, af, df, z, interactionInclusionProbs, interactionRetentionProbs)
{
  af@func1 <- if (is(af@func1, "InteractionFunction"))
    addTreatmentToInteractionFunction(C, af@func1, df, z, interactionInclusionProbs, interactionRetentionProbs)
  else
    addTreatmentToAdditiveFunction(C, af@func1, df, z, interactionInclusionProbs, interactionRetentionProbs)
  
  af@func2 <- if (is(af@func2, "InteractionFunction"))
    addTreatmentToInteractionFunction(C, af@func2, df, z, interactionInclusionProbs, interactionRetentionProbs)
  else
    addTreatmentToAdditiveFunction(C, af@func2, df, z, interactionInclusionProbs, interactionRetentionProbs)
  
  af
}

if (FALSE) addTreatmentMainEffectToGeneralizedAdditiveFunction <- function(C, gaf, z, mean, scale)
{
  treatmentTerm <- generateUnivariateTerm(".z", z, c(step.discrete = 1))
  # treatmentTerm@functions[[1L]] <- shuffleBaseline(treatmentTerm@functions[[1L]])
  treatmentTerm@functions[[1L]]@pars[2L] <- rt(1L, C$BF_DF, mean) * scale
  environment(treatmentTerm@functions[[1L]]@f)$pars <- treatmentTerm@functions[[1L]]@pars
  
  termAdded <- FALSE
  for (i in seq_along(gaf@functions)) {
    if (isConstantAdditiveFunction(gaf@functions[[i]])) next
    if (is(gaf@functions[[i]], "InteractionFunction")) next
    if (gaf@functions[[i]]@transformation@type != transformationsMap[["identity"]]) next
    
    gaf@functions[[i]]@terms <- append(treatmentTerm, gaf@functions[[i]]@terms)
    gaf@functions[[i]]@weights <- append(gaf@functions[[i]]@weights[1L], gaf@functions[[i]]@weights)
    termAdded <- TRUE
  }
  if (i == length(gaf@functions) && termAdded == FALSE) {
    af <- new("AdditiveFunction", transformation = newTransformation("identity"), 
              inputScale = 1, outputScale = 1, terms = list(treatmentTerm), weights = 1,
              operator = operatorsMap[["+"]])
    gaf@functions <- append(af, gaf@functions)
    gaf@weights <- append(gaf@weights[[1L]], gaf@weights)
  }
  gaf
}

if (FALSE) addTreatmentMainEffectToGeneralizedAdditiveFunction2 <- function(gaf, z, value)
{
  treatmentTerm <- generateUnivariateTerm(".z", z, c(step.discrete = 1))
  pars <- c(0, value)
  attr(pars, "map") <- c("ctl", "trt")
  treatmentTerm@functions[[1L]]@pars <- pars
  environment(treatmentTerm@functions[[1L]]@f)$pars <- pars
  
  termAdded <- FALSE
  for (i in seq_along(gaf@functions)) {
    if (isConstantAdditiveFunction(gaf@functions[[i]])) next
    if (is(gaf@functions[[i]], "InteractionFunction")) next
    if (gaf@functions[[i]]@transformation@type != transformationsMap[["identity"]]) next
    
    gaf@functions[[i]]@terms <- append(treatmentTerm, gaf@functions[[i]]@terms)
    gaf@functions[[i]]@weights <- append(1 / gaf@functions[[i]]@outputScale, gaf@functions[[i]]@weights)
    termAdded <- TRUE
    break
  }
  if (i == length(gaf@functions) && termAdded == FALSE) {
    af <- new("AdditiveFunction", transformation = newTransformation("identity"), 
              inputScale = 1, outputScale = 1, terms = list(treatmentTerm), weights = 1,
              operator = operatorsMap[["+"]])
    gaf@functions <- append(af, gaf@functions)
    gaf@weights <- append(1, gaf@weights)
  }
  gaf
}

addTreatmentMainEffectToGeneralizedAdditiveFunction3 <- function(C, gaf, z, c0, c1)
{
  treatmentTerm <- generateUnivariateTerm(C, ".z", z, c(step.discrete = 1))
  pars <- c(0, c1 - c0)
  attr(pars, "map") <- c("ctl", "trt")
  treatmentTerm@functions[[1L]]@pars <- pars
  environment(treatmentTerm@functions[[1L]]@f)$pars <- pars
  
  termAdded <- FALSE
  constantChanged <- FALSE
  for (i in seq_along(gaf@functions)) {
    if (isConstantAdditiveFunction(gaf@functions[[i]])) {
      constFunc <- gaf@functions[[1L]]@terms[[1L]]@functions[[1L]]
      gaf@functions[[1L]]@terms[[1L]]@functions[[1L]] <-
        setConstantForFunction(constFunc, constFunc@pars + c0)
      constantChanged <- TRUE
      next
    }
    if (is(gaf@functions[[i]], "InteractionFunction")) next
    if (gaf@functions[[i]]@transformation@type != transformationsMap[["identity"]]) next
    
    gaf@functions[[i]]@terms <- append(treatmentTerm, gaf@functions[[i]]@terms)
    gaf@functions[[i]]@weights <- append(1 / gaf@functions[[i]]@outputScale, gaf@functions[[i]]@weights)
    termAdded <- TRUE
    break
  }
  if (i == length(gaf@functions) && constantChanged == FALSE) {
    gaf <- addConstantTermToGeneralizedAdditiveFunction(gaf, c0)
    i <- i + 1L
  }
  if (i == length(gaf@functions) && termAdded == FALSE) {
    af <- new("AdditiveFunction", transformation = newTransformation("identity"), 
              inputScale = 1, outputScale = 1, terms = list(treatmentTerm), weights = 1,
              operator = operatorsMap[["+"]])
    gaf@functions <- append(af, gaf@functions)
    gaf@weights <- append(1, gaf@weights)
  }
  gaf
}


addTreatmentToAdditiveFunction <- function(C, af, df, z, interactionInclusionProbs, interactionRetentionProbs)
{
  if (!is.null(names(interactionInclusionProbs))) {
    index <- which(names(interactionInclusionProbs) == transformationsMap[af@transformation@type])
    
    interactionInclusionProbs <- interactionInclusionProbs[[index]]
    interactionRetentionProbs <- interactionRetentionProbs[[index]]
  }
  
  treatmentTerms <- list()
  
  baseTermIndex <- 1L
  if (!is.null(interactionInclusionProbs) && length(interactionInclusionProbs) > 0) for (i in seq_along(interactionInclusionProbs))
  {
    marginalInteractionProb <- interactionInclusionProbs[[i]]
    
    numTerms <- 0L
    for (termIndex in seq.int(baseTermIndex, length(af@terms))) {
      if (getTermSize(af@terms[[termIndex]]) == i) numTerms <- numTerms + 1L
      else break
    }
    if (numTerms == 0L) next
    
    if (is.list(marginalInteractionProb) && !is.null(names(marginalInteractionProb[[2]]))) {
      globalTermProbs <- marginalInteractionProb[[2L]]
      termProbs <- rep_len(0, numTerms)
      
      for (j in seq.int(baseTermIndex, baseTermIndex + numTerms - 1L)) {
        termNames <- getTermNames(af@terms[[j]])
        if (any(termNames %in% names(globalTermProbs))) {
          termProbs[j] <- mean(unname(globalTermProbs[match(termNames, names(globalTermProbs))]))
        }
      }
      marginalInteractionProb[[2L]] <- termProbs / sum(termProbs)
    }
    
    interactionIndices <- baseTermIndex + sampleIndices(numTerms, marginalInteractionProb) - 1L
    
    for (j in seq_along(interactionIndices)) {
      interactionTerm <- interactTermWithTreatment(C, af@terms[[interactionIndices[j]]], df, z, interactionRetentionProbs)
      if (!is.null(interactionTerm))
        treatmentTerms <- append(treatmentTerms, interactionTerm)
    }
    
    baseTermIndex <- termIndex + 1L
    
    if (baseTermIndex > length(af@terms)) break
  }
  
  af@terms <- append(af@terms, treatmentTerms)
  
  #effectiveNumTerms <- 0
  #for (i in seq_along(af@terms))
  #  effectiveNumTerms <- effectiveNumTerms + countEffectiveNumTerms(af@terms[[i]])
  
  #af@weights <- rep_len(1 / sqrt(effectiveNumTerms), length(af@terms))
  #af@weights <- c(af@weights, rt(length(treatmentTerms), 1) / 4.5)
  
  af@weights <- getWeightsForTerms(af@terms)
  
  af
}

evaluateAdditiveFunctionToDataframe <- function(af, x, mainEffectsOnly = FALSE)
{
  if (!is.data.frame(x) && is.null(names(x))) stop("argument to evaluate must be a data frame")
  
  if (is(af, "InteractionFunction")) {
    df1 <- evaluateAdditiveFunctionToDataframe(af@func1, x, mainEffectsOnly)
    df2 <- evaluateAdditiveFunctionToDataframe(af@func2, x, mainEffectsOnly)
    
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
    return(result)
  }
  
  for (i in seq_along(x))
    if (!is.discrete(x[[i]])) x[[i]] <- af@inputScale * x[[i]]
  
  if (af@operator == operatorsMap[["+"]]) {
    result <- list()
    for (i in seq_along(af@terms)) {
      df.i <- evaluateAdditiveTermToDataframe(af@terms[[i]], x)
      df.i <- lapply(df.i, function(col) col * af@outputScale * af@weights[i])
      result <- append(result, df.i)
    }
  } else {
    result <- rep(1, NROW(x))
    if (length(af@terms) == 0L)
      result <- 0
    else for (i in seq_along(af@terms)) {
      result <- result * af@weights[i] * evaluate(af@terms[[i]], x)
    }
    result <- list(af@outputScale * result)
    #names(result) <- paste0(sapply(af@terms, function(term.i) term.i@name), collapse = " * ")
    names(result) <- format(af)
  }
  if (af@transformation@type == transformationsMap[["exp"]]) {
    numCols <- length(result)
    
    result.exp <- list(C = rep_len(1, nrow(x)))
    result.exp <- append(result.exp, result)
    
    result.sq <- vector("list", choose(numCols, 2))
    names.sq <- character(choose(numCols, 2))
    index <- 1L
    for (i in seq_len(numCols)) {
      result.sq[[index]] <- 0.5 * result[[i]]^2
      names.sq[index] <- paste0(names(result)[i], "^2")
      index <- index + 1L
      
      if (i < numCols - 1L) for (j in seq.int(i + 1L, numCols)) {
        result.sq[[index]] <- result[[i]] * result[[j]]
        names.sq[index] <- paste0(names(result)[i], " * ", names(result)[j])
        index <- index + 1L
      }
    }
    names(result.sq) <- names.sq
    result <- result.exp
  }
  result
}

evaluateGeneralizedAdditiveFunctionToDataframe <- function(gaf, x, mainEffectsOnly = FALSE)
{
  if (!is.data.frame(x) && is.null(names(x))) stop("argument to evaluate must be a data frame")
  
  result <- list()
  for (i in seq_along(gaf@functions)) {
    af.result <- evaluateAdditiveFunctionToDataframe(gaf@functions[[i]], x)
    w <- if (is.function(gaf@weights[[i]])) gaf@weights[[i]](x) else gaf@weights[[i]]
    af.result <- lapply(af.result, function(col) col * w)
    result <- append(result, af.result)
  }
  attr(result, "class") <- "data.frame"
  attr(result, "row.names") <- seq_len(NROW(x))
  result
}

getTermNamesForGeneralizedAdditiveFunction <- function(gaf)
{
  result <- character()
  for (i in seq_along(gaf@functions)) {
    result <- c(result, getTermNamesForAdditiveFunction(gaf@functions[[i]]))
  }
  unique(result)
}

getTermNamesForAdditiveFunction <- function(af)
{
  if (is(af, "InteractionFunction"))
    return(c(getTermNamesForGeneralizedAdditiveFunction(af@func1), getTermNamesForGeneralizedAdditiveFunction(af@func2)))
  
  result <- character()
  for (i in seq_along(af@terms))
    result <- c(result, getTermNames(af@terms[[i]]))
  
  result
}

getMainEffectsForGeneralizedAdditiveFunction <- function(gaf)
{
  result <- character()
  for (i in seq_along(gaf@functions)) {
    result <- c(result, getMainEffectsForAdditiveFunction(gaf@functions[[i]]))
  }
  unique(result)
}

getMainEffectsForAdditiveFunction <- function(af)
{
  if (is(af, "InteractionFunction") || af@operator != operatorsMap[["+"]] ||
      af@transformation@type != transformationsMap[["identity"]]) return(character())
  
  result <- character()
  for (i in seq_along(af@terms)) {
    term.i <- if (isTreatmentInteractionTerm(af@terms[[i]])) af@terms[[i]]@term1 else af@terms[[i]]
    if (is(term.i, "UnivariateAdditiveTerm") && term.i@name != "constant" && term.i@name != ".z") result <- c(result, term.i@name)
  }
  result
}

scaleGeneralizedAdditiveFunction <- function(gaf, mean, scale)
{
  scaleFound <- FALSE
  for (i in seq_along(gaf@functions)) {
    af <- gaf@functions[[i]]
    if (isConstantAdditiveFunction(af)) {
      newConst <- scale * gaf@functions[[i]]@terms[[1L]]@functions[[1L]]@pars + (1 - scale) * mean
      gaf@functions[[i]]@terms[[1L]]@functions[[1L]]@pars <- newConst
      environment(gaf@functions[[i]]@terms[[1L]]@functions[[1L]]@f)$pars <- newConst
      scaleFound <- TRUE
    } else {
      gaf@functions[[i]]@outputScale <- scale * gaf@functions[[i]]@outputScale
    }
  }
  if (!scaleFound) gaf <- addConstantTermToGeneralizedAdditiveFunction(gaf, mean * (1 - scale))
  gaf
}

countInteractionsForGeneralizedAdditiveFunction <- function(gaf, countTreatment = FALSE)
{
  result <- 0L
  for (i in seq_along(gaf@functions)) {
    func.i <- gaf@functions[[i]]
    if (!is(func.i, "AdditiveFunction") || func.i@operator != operatorsMap[["+"]] ||
        func.i@transformation@type != transformationsMap[["identity"]]) next
    
    for (j in seq_along(func.i@terms)) {
      term.j <- func.i@terms[[j]]
      if (is(term.j, "UnivariateAdditiveTerm")) next
      
      if (countTreatment == FALSE) {
        if (!isTreatmentInteractionTerm(term.j)) result <- result + 1L
      } else {
        if (isTreatmentInteractionTerm(term.j)) result <- result + 1L
      }
    }
  }
  result
}

countFunctionalFormsForGeneralizedAdditiveFunction <- function(gaf, functions)
{
  result <- 0L
  for (i in seq_along(gaf@functions)) {
    func.i <- gaf@functions[[i]]
    if (!is(func.i, "AdditiveFunction") || func.i@operator != operatorsMap[["+"]] ||
        func.i@transformation@type != transformationsMap[["identity"]]) next
    
    for (j in seq_along(func.i@terms)) {
      term.j <- func.i@terms[[j]]
      if (!is(term.j, "UnivariateAdditiveTerm") || term.j@name == ".z") next
      
      for (k in seq_along(term.j@functions)) {
        func.k <- term.j@functions[[k]]
        if (baseFunctions[func.k@type] %in% functions) result <- result + 1L
      }
    }
  }
  result
}

countPurelyLinearTermsForGeneralizedAdditiveFunction <- function(gaf)
{
  result <- 0L
  for (i in seq_along(gaf@functions)) {
    func.i <- gaf@functions[[i]]
    if (!is(func.i, "AdditiveFunction") || func.i@operator != operatorsMap[["+"]] ||
        func.i@transformation@type != transformationsMap[["identity"]]) next
    
    for (j in seq_along(func.i@terms)) {
      term.j <- func.i@terms[[j]]
      if (!is(term.j, "UnivariateAdditiveTerm") || term.j@name == ".z") next
      
      if (length(term.j@functions) > 1L) next
      
      if (term.j@functions[[1L]]@type == baseFunctionsMap["linear"]) result <- result + 1L
    }
  }
  result
}

countNumberOfMainEffectsInGeneralizedAdditiveFunction <- function(gaf){
  result <- 0L
  for (i in seq_along(gaf@functions)) {
    func.i <- gaf@functions[[i]]
    if (!is(func.i, "AdditiveFunction") || func.i@operator != operatorsMap[["+"]] ||
        func.i@transformation@type != transformationsMap[["identity"]]) next
    
    for (j in seq_along(func.i@terms)) {
      term.j <- func.i@terms[[j]]
      if (!is(term.j, "UnivariateAdditiveTerm") || term.j@name == ".z") next
      
      result <- result + 1L
    }
  }
  result

}

