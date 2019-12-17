aciccomp 2016 Fitted Values
===========================

This directory contains of the fitted values for the 15 original methods and 9 post-competition additions in the paper "Automated versus do-it-yourself methods for causal inference: Lessons learned from a data analysis competition", Statist. Sci., Volume 34, Number 1 (2019), 43-68.

## Format

Each RData file contains an object named `results', which is itself a 7700 x matrix containing columns `sate` for the estimate, `sate.lower` for the lower 95% confidence bound, and `sate.upper` for the upper.

In order to avoid name collisions, it is recommended to load each object in a special environment. For example:

```R
loadEnv <- new.env(parent = baseenv())

load("sl_bart_tmle.RData", envir = loadEnv)
```

The rownames of the result matrix correspond to the parameter set and simulation number, iterating first through the 77 cases and then through the 100 replications. Maps between row numbers and parameter/simulation number are given by:

```R
getParametersForIter <- function(i.row) {
  c((i.row - 1L) %% nrow(aciccomp2016::parameters_2016) + 1L, (i.row - 1L) %/% nrow(aciccomp2016::parameters_2016) + 1L)
}

getIterForParameters <- function(i.par, i.sim)
{
  i.par + (i.sim - 1L) * nrow(aciccomp2016::parameters_2016)
}
```
