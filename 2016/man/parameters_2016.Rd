\name{parameters_2016}
\docType{data}
\alias{parameters_2016}
\title{Parameters Data for the 2016 ACIC Competition}
\description{
  Data set containg the parameters used to generate data for the 2016 Atlantic Causal Inference Conference
  competition.
}
\usage{parameters_2016}
\format{A data frame describing 77 scenarios that vary across 6 features.
  \enumerate{
  \item \code{model.trt} - Function distribution over the treatment assignment mechanism. Can be \code{"linear"},
                           \code{"polynomial"}, or \code{"step"}.
  \item \code{root.trt} - Baseline probability of receiving treatment.
  \item \code{overlap.trt} - Term that controls the addition of overlap-penalizing terms that forcibly exclude
                             observations from the treatment group by carving out hyper-rectangles of the
                             covariate space and assigning their treatment probability to 0. Can be \code{"full"}
                             for complete overlap, \code{"one-term"} for adding a single function as described
                             above, or \code{"two-term"} for adding two. Two-terms were not used in the
                             competition and is not thoroughly tested.
  \item \code{model.rsp} - Function distribution over the response surface. Can be \code{"linear"},
                           \code{"polynomial"}, \code{"step"}, or \code{"exponential"}.
  \item \code{alignment} - A numeric value that determines the degree to which terms from the treatment
                           assignment function appear in response surface function.
  \item \code{te.hetero} - A term that controls the degree of treatment effect heterogeneity. Can be \code{"none"}
                           for parallel surfaces, \code{"med"} or \code{"high"}. Higher heterogeneity is achieved
                           by selectively interacting terms from the response surface with a treatment indicator.
  }
}
\source{
  Original release.
}
\references{
  Dorie V., Hill J., Shalit U., Scott M. and Cervone D. (2017) Automated versus do-it-yourself methods for causal
  inference: Lessons learned from a data analysis competition, preprint arXiv \url{https://arxiv.org/abs/1707.02641}.
}
\keyword{datasets}

