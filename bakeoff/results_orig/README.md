aciccomp 2016 Fitted Values
===========================

This directory contains of the fitted values for the 15 original methods and 9 post-competition additions in the paper "Automated versus do-it-yourself methods for causal inference: Lessons learned from a data analysis competition", Statist. Sci., Volume 34, Number 1 (2019), 43-68.

## Format

Each rds file is a 7700 x matrix containing columns `satt` for the estimate, `satt.lower` for the lower 95% confidence bound, and `satt.upper` for the upper. As an example of usage:

```R
> sl_bart_tmle_results <- readRDS("sl_bart_tmle.rds")

> dim(sl_bart_tmle_results)
[1] 7700    3
> head(sl_bart_tmle_results)
                                                       satt satt.lower satt.upper
linear_0.35_one-term_linear_0.75_high_001          2.710912   2.432793   2.989031
polynomial_0.35_one-term_exponential_0.75_none_001 4.643160   4.571654   4.714667
linear_0.35_one-term_linear_0.75_none_001          6.591090   6.524476   6.657703
polynomial_0.35_full_exponential_0.75_high_001     4.036371   3.848314   4.224428
linear_0.35_one-term_exponential_0.75_high_001     3.430520   3.125533   3.735508
polynomial_0.35_one-term_linear_0.75_high_001      2.995059   2.835973   3.154145
```

The rownames of the result matrix correspond to the parameter set and simulation number, iterating first through the 77 cases and then through the 100 replications. Maps between row numbers and parameter/simulation number are given by:

```R
getParametersForIter <- function(i.row) {
  c((i.row - 1L) %%  nrow(aciccomp2016::parameters_2016) + 1L,
    (i.row - 1L) %/% nrow(aciccomp2016::parameters_2016) + 1L)
}

getIterForParameters <- function(i.par, i.sim)
{
  i.par + (i.sim - 1L) * nrow(aciccomp2016::parameters_2016)
}
```
