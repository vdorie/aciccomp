\name{dgp_2016}
\title{Data Generating Process for the 2016 ACIC Competition}
\alias{dgp_2016}
\description{
  Applies the data generating process used in the Atlantic Causal Inference Competition of 2016.
}
\usage{
  dgp_2016(x, parameters, random.seed,
           constants = constants_2016(),
           extraInfo = FALSE)
}
\arguments{
  \item{x}{Input data in the form of a data frame, most likely \code{\link{input_2016}}.}
  \item{parameters}{A named list containing elements in the form of \code{\link{parameters_2016}},
                    a row of the same object, or an integer specifying which row of
                    \code{parameters_2016} is to be used; see that page for details.}
  \item{random.seed}{A list of arguments to be used in a call to \code{\link{set.seed}} or an
                     integer between 1 and 100 specifying the random seed associated with an
                     iteration from the competition.}
  \item{constants}{A named list containing elements as returned by \code{\link{constants_2016}};
                   see there for details.}
  \item{extraInfo}{Boolean determining if additional information is to be returned, including the
                   treatment and control response surfaces, the transformed input data, and
                   whether or not a simulation would have been deemed interesting enough to
                   include in the competition.}
}
\details{
  Creates a causal inference problem by taking the input \code{x} and using the passed in
  \code{parameters} to generate a treatment assignment mechanism (probability of treatment for
  each individual), response surface (expected value under treatment and control), and finally
  observed data. The \code{parameters} provide high-level controls to adjust the result for
  causal inference features that may be of interest, while \code{constants} control at a lower
  level the parameters of generated functions.

  \subsection{Generalized Additive Functions}{
    The 2016 competition used a unique set of software that was internally described as
    ``Generalized Additive Functions'' (GAFs). A GAF consists of many small functions applied to
    various features/columns of the input that are added together or interacted with each other.
    The complete sum may then be passed through a link function to achieve a result in a
    transformed space. The small functions are randomly derived from a library of functions,
    so that the general features of the result can vary according to high level parameters.
    
    This package reproduces GAFs as they were used in the 2016 contest without the intention
    that they be further applied. It may be possible to use \code{dgp_2016} with different
    input data and changes to the \code{constants} should propogate, however these extensions
    will not be widely supported.
  }
}
\value{
  A named list containing:
  
  \item{\code{z}}{Vector of treatment assignments. If \code{extraInfo} is \code{FALSE}, \code{z}
                  contains 0s and 1s. If \code{TRUE}, \code{z} is a factor with levels \code{ctl}
                  and \code{trt}.}
  \item{\code{y}}{Vector of observed response variables, \eqn{y(z)}.}
  \item{\code{y.0}}{Vector of response variables under the control condition, \eqn{y(0)}.}
  \item{\code{y.1}}{Vector of response variables under the treatment condition, \eqn{y(1)}.}
  \item{\code{mu.0}}{Vector of expected response under the control condition, \eqn{E[Y(0)]}.}
  \item{\code{mu.1}}{Vector of expected response under the treatment condition, \eqn{E[Y(1)]}.}
  \item{\code{e}}{Vector of propensity scores, \eqn{P(Z = 1)}.}
  \item{\code{f.z}}{Optional - the GAF for the treatment assignment mechanism.}
  \item{\code{f.y}}{Optional - the GAF for the response surface.}
  \item{\code{x}}{Optional - the transformed input passed to \code{f.z} and \code{f.y}.}
  \item{\code{valid}}{Optional - boolean if the simulation would be rejected as "uninteresting".}
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\references{
  Dorie V., Hill J., Shalit U., Scott M. and Cervone D. (2017) Automated versus do-it-yourself methods for causal
  inference: Lessons learned from a data analysis competition, preprint arXiv \url{https://arxiv.org/abs/1707.02641}.
}
\examples{
\dontrun{
# to test a method
ate <- matrix(NA, 77, 100)
for (i in seq_len(77)) {
  for (j in seq_len(100)) {
    sim <- dgp_2016(input_2016, i, j)
    df <- input_2016
    df$y <- sim$y
    df$z <- sim$z
    fit <- lm(y ~ ., df)
    ate[i,j] <- coef(fit)["z"]
  }
}

## undocumented features, getting closest approximate linear model
sim <- dgp_2016(input_2016, 1, 1, extraInfo = TRUE)

e <- aciccomp:::evaluate(sim$f.z, sim$x)
x.z.approx <- aciccomp:::evaluateGeneralizedAdditiveFunctionToDataframe(sim$f.z, sim$x)

x.temp <- sim$x
x.temp$.z <- sim$z
x.y.approx <- aciccomp:::evaluateGeneralizedAdditiveFunctionToDataframe(sim$f.y, x.temp)
}
}
\keyword{simulation}
\keyword{causal inference}

