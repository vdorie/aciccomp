\name{dgp_2017}
\title{Data Generating Process for the 2017 ACIC Competition}
\alias{dgp_2017}
\description{
  Applies the data generating process used in the Atlantic Causal Inference Competition of 2017.
}
\usage{
  dgp_2017(parameters, random.seed)
}
\arguments{
  \item{parameters}{A named list containing elements in the form of \code{\link{parameters_2017}},
                    a row of the same object, or an integer specifying which row of
                    \code{parameters_2017} is to be used; see that page for details.}
  \item{random.seed}{A list of arguments to be used in a call to \code{\link{set.seed}} or an
                     integer between 1 and 250 specifying the random seed associated with an
                     iteration from the competition.}
}
\details{
  Creates a causal inference problem by using the passed in \code{parameters} to generate a
  treatment assignment mechanism (probability of treatment for each individual), treatment effect,
  and finally observed data.
  
  Data are derived from the covariates in \code{\link{input_2017}}.
}
\value{
  A named list containing:
  
  \item{\code{z}}{Vector of treatment assignments.}
  \item{\code{y}}{Vector of observed response variables, $y(z)$.}
  \item{\code{alpha}}{Vector of individual level treatment effects, $y(1) - y(0)$.}
}
\author{
  P Richard Hahn \email{prhahn@asu.edu} and Vincent Dorie \email{vdorie@gmail.com}.
}
\references{
  Hahn P.R., Dorie V. and Murray J.R. (2018) Atlantic Causal Inference Conference (ACIC) Data Analysis Challenge
  2017, online pdf \url{https://math.la.asu.edu/~prhahn/debrief.pdf}.
}
\examples{
\dontrun{
sim <- dgp_2017(1, 1)
att <- mean(sim$alpha[sim$z == 1])
}
}
\keyword{simulation}
\keyword{causal inference}

