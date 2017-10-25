\name{constants_2016}
\title{Constants Used in DGP for ACIC Competition 2016}
\alias{constants_2016}
\description{
  Returns or sets elements of a named list containing all of the constants required to run the
  data generating processes for the 2016 ACIC Competition.
}
\usage{
  constants_2016(...)
}
\arguments{
  \item{...}{Options from the list below.}
}
\details{
  Returns default values or sets them, as appropriate. Minimal error checking is performed.
}
\value{
  \item{\code{RSP_INPUT_SCALE}}{Scaling factor applied to covariates before evaluating the
                                response function.}
  \item{\code{RSP_OUTPUT_SHAPE_1}}{The first shape parameter in a beta-prime used to generate
                                   the output scale of the response function.}
  \item{\code{RSP_OUTPUT_RATE}}{The inverse scale parameter in a beta-prime used to generate
                                the output scale of the response function.}
  \item{\code{RSP_OUTPUT_SHAPE_2}}{The second shape parameter in a beta-prime used to generate
                                   the output scale of the response function.}
  \item{\code{TRT_INPUT_SCALE}}{Scaling factor applied to covariates before evaluating the
                                treatment assignment function.}
  \item{\code{TRT_OUTPUT_SCALE}}{Scaling factor applied to result of the treatment assignment
                                 function.}
  \item{\code{TRT_BIAS_SCALE}}{Approximate scale for treatment biasing functions when
                               \code{overlap} parameter is not \code{"full"}.}
  \item{\code{RSP_SIGMA_Y}}{Scale of noise added to response.}
  
  \item{\code{BF_CONSTANT_SCALE}}{Scale of constant base function parameter.}
  \item{\code{BF_LINEAR_SCALE}}{Scale of linear base function parameter.}
  \item{\code{BF_QUADRATIC_SHAPE_1}}{First shape parameter used to generate quadratic base
                                     function root parameter.}
  \item{\code{BF_QUADRATIC_SHAPE_2}}{Second shape parameter used to generate quadratic base
                                     function root parameter.}
  \item{\code{BF_QUADRATIC_RATE}}{Rate parameter used to generate quadratic base function
                                  root parameter.}
  \item{\code{BF_QUADRATIC_SCALE}}{Scale of quadratic base function parameter.}
  \item{\code{BF_CUBIC_SHAPE}}{Shape parameter used to generate cubic base function root
                               parameters.}
  \item{\code{BF_CUBIC_RATE}}{Rate parameter used to generate cubic base function root
                              parameters.}
  \item{\code{BF_CUBIC_SCALE}}{Scale of cubic base function parameter.}
  \item{\code{BF_CONTINUOUS_SCALE}}{Scale parameter shared by continuous base functions.}
   
  \item{\code{BF_STEP_SHAPE}}{Shape parameter for step base functions.}
  \item{\code{BF_STEP_CONSTANT_SCALE}}{Scale of step-wise constant base function parameter.}
  \item{\code{BF_STEP_LINEAR_SCALE}}{Scale of piece-wise linear base function parameter.}
  \item{\code{BF_SIGMOID_SHAPE_1}}{First shape parameter used to generate sigmoid base function
                                   parameters.}
  \item{\code{BF_SIGMOID_RATE_1}}{First rate parameter used to generate sigmoid base function
                                   parameters.}
  \item{\code{BF_SIGMOID_SHAPE_2}}{Second shape parameter used to generate sigmoid base function
                                   parameters.}
  \item{\code{BF_SIGMOID_RATE_2}}{Second rate parameter used to generate sigmoid base function
                                   parameters.}
  \item{\code{BF_QUANTILE_SHAPE_1}}{First shape parameter used to generate quantile base
                                    function cutoff.}
  \item{\code{BF_QUANTILE_SHAPE_2}}{Second shape parameter used to generate quantile base
                                    function cutoff.}
  
  \item{\code{BF_TWEAK_SIGN_PROB}}{Probability of changing sign when copy/modifying base
                                   function.}
  \item{\code{BF_TWEAK_NORMAL_SCALE}}{Scale of normal noise added to unconstrained base function
                                      parameters when copy/modifying.}
  \item{\code{BF_TWEAK_GAMMA_SHAPE}}{Shape parameter of positive noise added to constrained base
                                     function parameters when copying/modifying.}
  \item{\code{BF_TWEAK_GAMMA_RATE}}{Rate parameter of positive noise added to constrained base
                                    function parameters when copying/modifying.}
  
  \item{\code{TRT_BF_DF}}{Base function degrees of freedom when generating treatment assignment
                          mechanism.}
  \item{\code{RSP_BF_DF}}{Base function degrees of freedom when generating response surface.}

  \item{\code{TRT_LINEAR_SCALE_SHAPE_1}}{First scale parameter used to generate overall scale
                                         of treatment assignment mechanism.}
  \item{\code{TRT_LINEAR_SCALE_SHAPE_2}}{Second scale parameter used to generate overall scale
                                         of treatment assignment mechanism.}
  \item{\code{TRT_LINEAR_SCALE_RATE}}{Rate parameter used to generate overall scale of treatment
                                      assignment mechanism.}
  
  \item{\code{RSP_EXP_SCALE_SHAPE}}{Shape parameter used when generating scale factor for
                                    exponential functions.}
  \item{\code{RSP_EXP_SCALE_RATE}}{Rate parameter used when generating scale factor for
                                   exponential functions.}
  \item{\code{RSP_EXP_WEIGHT_SHAPE}}{Shape parameter used when generating relative weight factor
                                     for exponential functions.}
  \item{\code{RSP_EXP_WEIGHT_RATE}}{Rate parameter used when generating relative weight factor
                                     for exponential functions.}

  \item{\code{RSP_TE_MEAN}}{Expected value for population average treatment effect.}
  \item{\code{RSP_TE_SCALE}}{Scale factor for population average treatment effect.}
  \item{\code{RSP_TE_DF}}{Degrees of freedom for population average treatment effect.}
  
  \item{\code{SPARSE_COVARIATE_WEIGHT}}{Weight of inclusion for sparse, discrete covariates.}
  \item{\code{CONTINUOUS_COVARIATE_WEIGHT}}{Weight of inclusion for continuous covariates.}
  \item{\code{DEFAULT_COVARIATE_WEIGHT}}{Default weight of inclusion for covariates.}
  
  \item{\code{TRT_BASELINE_SHIFT}}{Function used to derive a scale when generating a baseline
                                   treatment probability from \code{root.trt}} 

  \item{\code{BASE_FUNCTION_DIST_LIN}}{Base function distribution containing only linear
                                       functions.}
  \item{\code{BASE_FUNCTION_DIST_POLY}}{Base function distribution containing linear, quadratic,
                                        and cubic functions.}
  \item{\code{BASE_FUNCTION_DIST_STEP}}{Base function distribution containing linear, step-wise
                                        constant, and piece-wise linear functions.}
  \item{\code{BASE_FUNCTION_DIST_EXP}}{Base function distribution containing third order
                                       polynomials to be used in exponential functions.}
  \item{\code{dist.lin}}{Function distribution for purely linear treatment or response.}
  \item{\code{dist.int}}{Function distribution with linear terms and interactions.}
  \item{\code{dist.pure.poly}}{Function distribution with quadratic terms and no interactions.}
  \item{\code{dist.poly}}{Function distribution with cubic terms and interactions.}
  \item{\code{dist.step}}{Function distribution with linear terms, step-wise constant terms,
                          and interactions.}
  \item{\code{dist.exp}}{Function distribution with quadratic terms and interactions appropriate
                         for use with exponential link functions.}
  \item{\code{dist.bias1}}{Function distribution over treatment assignment biasing functions.}
  \item{\code{dist.bias2}}{Function distribution over treatment assignment biasing functions.}
  \item{\code{dist.hetero.med}}{Function distribution specifying interaction retention
                                probabilities for medium degrees of treatment effect
                                heterogeneity.}
  \item{\code{dist.hetero.high}}{Function distribution specifying interaction retention
                                 probabilities for high degrees of treatment effect
                                 heterogeneity.}
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\keyword{simulation}
\keyword{causal inference}

