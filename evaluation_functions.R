evaluation_functions <- list(
  ## not actually bias, obviously
  bias = function(data, results, results.ind) {
    (mean(data$y[data$z == 1]) - mean(data$y[data$z == 0])) - results$est
  },
  ## and obviously not coverage
  coverage = function(data, results, results.ind) {
    te <- mean(data$y[data$z == 1]) - mean(data$y[data$z == 0])
    results$ci_lower <= te & te <= results$ci_upper
  })
