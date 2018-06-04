\name{parameters_2017}
\docType{data}
\alias{parameters_2017}
\title{Parameters Data for the 2017 ACIC Competition}
\description{
  Data set containg the parameters used to generate data for the 2017 Atlantic Causal Inference Conference
  competition.
}
\usage{parameters_2017}
\format{A data frame describing 32 scenarios that vary across 4 features.
  \enumerate{
  \item \code{errors} - The type of measurent error/residual noise added to the data. Can be \code{"group_corr"},
                        \code{"heteroskedastic"}, \code{"iid"}, or \code{"nonadditive"}.
  \item \code{magnitude} - Size of the average treatment effect. Can be 0 for small and 1 for large.
  \item \code{noise}     - Size of the signal-to-noise ratio. Can be 0 for low signal and 1 for high.
  \item \code{confounding} - Strength of the relationship between covariates and treament model, which in turn
                             creates more or less confounding. Can be 0 for low confounding and 1 for high.
  }
}
\source{
  Original release.
}
\references{
  Hahn P.R., Dorie V. and Murray J.R. (2018) Atlantic Causal Inference Conference (ACIC) Data Analysis Challenge
  2017, online pdf \url{https://math.la.asu.edu/~prhahn/debrief.pdf}.
}
\keyword{datasets}

