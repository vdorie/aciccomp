.constants <- list(
  RSP_INPUT_SCALE  = 1 / 2.6,
  RSP_OUTPUT_SHAPE_1 = 32 * 3,
  RSP_OUTPUT_RATE  = 3,
  RSP_OUTPUT_SHAPE_2 = 5,
  TRT_INPUT_SCALE  = 1 / 2.6,
  TRT_OUTPUT_SCALE = 1.8,
  TRT_BIAS_SCALE   = -5.5,
  RSP_SIGMA_Y = 1,
  
  ## base function contants
  BF_CONSTANT_SCALE = 1 / 3,
  BF_LINEAR_SCALE   = 2 / 3,
  BF_QUADRATIC_SHAPE_1 = 1.15,
  BF_QUADRATIC_SHAPE_2 = 64,
  BF_QUADRATIC_RATE    = 32,
  BF_QUADRATIC_SCALE   = 1.275,
  BF_CUBIC_SHAPE = 0.5,
  BF_CUBIC_RATE = 32,
  BF_CUBIC_SCALE   = 0.475 * 3.3,
  BF_CONTINUOUS_SCALE = 8,
  
  BF_STEP_SHAPE = 1.15,
  BF_STEP_CONSTANT_SCALE = 2 / 3,
  BF_STEP_LINEAR_SCALE = 5.7 / 3,
  
  BF_SIGMOID_SHAPE_1 = 15 * 16,
  BF_SIGMOID_RATE_1 = 16,
  BF_SIGMOID_SHAPE_2 = 0.8 * 64,
  BF_SIGMOID_RATE_2 = 64,
  
  BF_QUANTILE_SHAPE_1 = 8,
  BF_QUANTILE_SHAPE_2 = 3.428,
  
  BF_TWEAK_SIGN_PROB = 0.75,
  BF_TWEAK_NORMAL_SCALE = 1 / 10,
  BF_TWEAK_GAMMA_SHAPE = 2 * 32,
  BF_TWEAK_GAMMA_RATE = 32,
  
  TRT_BF_DF = 5,
  RSP_BF_DF = 3,
  
  TRT_LINEAR_SCALE_SHAPE_1 = 32 * 3,
  TRT_LINEAR_SCALE_SHAPE_2 = 0.4,
  TRT_LINEAR_SCALE_RATE = 3,
  
  RSP_EXP_SCALE_SHAPE = 256 / 3,
  RSP_EXP_SCALE_RATE  = 256,
  
  RSP_EXP_WEIGHT_SHAPE = 2 * 128,
  RSP_EXP_WEIGHT_RATE  = 128,
  
  RSP_TE_MEAN      = 0.75,
  RSP_TE_SCALE     = 1 / 6,
  RSP_TE_DF        = 3,
  
  SPARSE_COVARIATE_WEIGHT = 1,
  CONTINUOUS_COVARIATE_WEIGHT = 6,
  DEFAULT_COVARIATE_WEIGHT = 1.5
)

linear.colName        <- paste0("linear_", c(0, 1))
quadratic.colName     <- paste0("quadratic_", c(0, 1))
cubic.colName         <- paste0("cubic_", c(0, 1))
step.constant.colName <- paste0("step.constant_", c(0, 1))
step.linear.colName   <- paste0("step.linear_", c(0, 1))

.constants$TRT_BASELINE_SHIFT <- function(x)
  (1 / 30) / (abs((0.5 - x) / 0.15)^1.75 + 1)

.constants$BASE_FUNCTION_DIST_LIN <- list(linear = 1)

.constants$BASE_FUNCTION_DIST_POLY <- list(
  linear = 1,
  quadratic = structure(c(0, 1 / 3), names = linear.colName),
  cubic = matrix(c(0, 0.10,
                   0, 0.25), 2, byrow = TRUE,
                 dimnames = list(quadratic.colName, linear.colName)))

.constants$BASE_FUNCTION_DIST_STEP <- list(
  linear = 1,
  step.constant = structure(c(0, 1 / 3), names = linear.colName),
  step.linear = matrix(c(0, 0.10,
                         0, 0.25), 2, byrow = TRUE,
                       dimnames = list(step.constant.colName, linear.colName)))

.constants$BASE_FUNCTION_DIST_EXP <- list(
  linear = 1,
  quadratic = structure(c(0, 0.25), names = linear.colName),
  step.constant = matrix(c(0, 0,
                           0, 0.25), 2, byrow = TRUE,
                         dimnames = list(quadratic.colName, linear.colName)))


.constants$dist.lin <- new("FunctionDistribution",
  termInclusionProb = list(c(6.5, 1.4), "placeholder"),
  baseFunctionDist = .constants$BASE_FUNCTION_DIST_LIN,
  interactionInclusionProbs = NULL,
  interactionRetentionProbs =
  c(linear = 1, step.discrete = 1,
    step.constant = 0.75, step.linear = 0.5,
    quadratic = 0.75, cubic = 0.5))

.constants$dist.int <- new("FunctionDistribution",
  termInclusionProb = list(c(6.5, 1.4), "placeholder"),
  baseFunctionDist = .constants$BASE_FUNCTION_DIST_LIN,
  interactionInclusionProbs = c(3 / 28, 2 / 56),
  interactionRetentionProbs = c(linear = 1, step.discrete = 1))

.constants$dist.pure.poly <- new("FunctionDistribution",
  termInclusionProb = list(c(6.5, 1.4), "placeholder"),
  baseFunctionDist = .constants$BASE_FUNCTION_DIST_POLY,
  interactionInclusionProbs = NULL,
  interactionRetentionProbs = NULL)

.constants$dist.poly <- new("FunctionDistribution",
  termInclusionProb = list(c(6.5, 1.4), "placeholder"),
  #termInclusionProb = list(c(10, 4), "placeholder"),
  baseFunctionDist = .constants$BASE_FUNCTION_DIST_POLY,
  interactionInclusionProbs = c(3 / 28, 2 / 56),
  #interactionInclusionProbs = c(8 / 45, 10 / 120),
  interactionRetentionProbs =
  c(linear = 1, step.discrete = 1,
    step.constant = 0.75, step.linear = 0.5,
    quadratic = 0.75, cubic = 0.5))

.constants$dist.step <- new("FunctionDistribution",
  termInclusionProb = list(c(6.5, 1.4), "placeholder"),
  baseFunctionDist = .constants$BASE_FUNCTION_DIST_STEP,
  interactionInclusionProbs = c(3 / 28, 2 / 56),
  interactionRetentionProbs =
  c(linear = 1, step.discrete = 1,
    step.constant = 0.75, step.linear = 0.5,
    quadratic = 0.75, cubic = 0.5))

.constants$dist.exp <- new("FunctionDistribution",
  termInclusionProb = list(c(4, 1), "placeholder"),
  baseFunctionDist = .constants$BASE_FUNCTION_DIST_EXP,
  interactionInclusionProbs = list(c(1.35, 0.15)),
  interactionRetentionProbs = c(linear = 1, step.discrete = 1, step.constant = 1,
  step.linear = 1, quadratic = 1, cubic = 1))

.constants$dist.bias1 <- new("FunctionDistribution",
  termInclusionProb = c(2.5, 0.30),
  baseFunctionDist = list(step.quantile = 1),
  interactionInclusionProbs = NULL,
  interactionRetentionProbs = NULL)

.constants$dist.bias2 <- new("FunctionDistribution",
  termInclusionProb = c(3.5, 0.45),
  baseFunctionDist = list(step.quantile = 1),
  interactionInclusionProbs = NULL,
  interactionRetentionProbs = NULL)

.constants$dist.hetero.med <- new("FunctionDistribution",
  termInclusionProb = numeric(0L),
  baseFunctionDist  = list(),
  interactionInclusionProbs = list(c(3, 0.4), 1 / 10, 1 / 100),
  interactionRetentionProbs = c(linear = 1, step.discrete = 1, step.constant = 0.75,
                                step.linear = 0.5, quadratic = 0.75, cubic = 0.5))

.constants$dist.hetero.high <- new("FunctionDistribution",
  termInclusionProb = numeric(0L),
  baseFunctionDist  = list(),
  interactionInclusionProbs = list(c(6, 1), 3, 2),
  interactionRetentionProbs = c(linear = 1, step.discrete = 1, step.constant = 0.75,
                                step.linear = 0.5, quadratic = 0.75, cubic = 0.5))



rm(linear.colName, quadratic.colName, cubic.colName, step.constant.colName, step.linear.colName)

constants_2016 <- function(...)
{
  result <- .constants
  args <- list(...)
  
  if (length(args) > 0L) {
    if (any(names(args) %not_in% names(result)))
      stop("unrecognized constants: ", paste0(evalx(names(args), x[x %not_in% names(result)]), collapse = ", "))
    
    result[names(args)] <- args
  }
  
  result
}

