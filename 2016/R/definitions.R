baseFunctions <-
  c("constant", "linear", "quadratic", "cubic", "step.constant", "step.linear", "step.discrete", "sigmoid", "step.quantile")
baseFunctionsMap <- as.list(seq_along(baseFunctions))
names(baseFunctionsMap) <- baseFunctions

transformations <-
  c("identity", "exp", "power", "sigmoid")
transformationsMap <- as.list(seq_along(transformations))
names(transformationsMap) <- transformations

operators <-
  c("+", "*")
operatorsMap <- as.list(seq_along(operators))
names(operatorsMap) <- operators

