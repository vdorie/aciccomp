## definitions required for subsequent 2016 functions
setClass("BaseFunction",
         slot = list(type = "integer",
                     pars = "numeric",
                     f    = "function"))
setClass("Transformation",
         slot = list(type = "integer",
                     pars = "numeric"))

setClass("AdditiveTerm")
setClass("UnivariateAdditiveTerm", contains = "AdditiveTerm",
         slots = list(
           name = "character",
           functions = "list"))
setClass("InteractionAdditiveTerm", contains = "AdditiveTerm",
         slots = list(
           term1 = "AdditiveTerm",
           term2 = "AdditiveTerm"))


setClass("GAFFunction",
         slots = list(
           transformation = "Transformation",
           isScale        = "logical"),
         prototype = list(
           transformation = new("Transformation", type = 1L, pars = numeric())))
setClass("AdditiveFunction", contains = "GAFFunction",
         slots = list(
           inputScale  = "numeric",
           outputScale = "numeric",
           terms     = "list",
           weights   = "numeric",
           operator = "integer"))

setClass("InteractionFunction", contains = "GAFFunction",
         slots = list(
           func1 = "GAFFunction",
           func2 = "GAFFunction"))

setClass("GeneralizedAdditiveFunction",
         slots = list(functions = "list",
                      weights = "list"))

setClassUnion("FunctionOrNumericOrList", c("function", "numeric", "list"))
setClassUnion("NumericOrListOrNull", c("numeric", "list", "NULL"))
setClass("FunctionDistribution",
         slots = list(termInclusionProb = "FunctionOrNumericOrList",
                      baseFunctionDist = "list",
                      interactionInclusionProbs = "NumericOrListOrNull",
                      interactionRetentionProbs = "NumericOrListOrNull"))

