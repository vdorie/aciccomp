newTransformation <- function(name, pars = NULL)
{
  if (is.null(pars)) {
    pars <- switch(name,
      identity = numeric(),
      exp      = numeric(),
      power    = 0.5,
      sigmoid  = c(rgamma(1L, 5 * 64, 64), sample(c(-1, 1), 1) * rgamma(1L, 0.5 * 128, 128)))
  }
  methods::new("Transformation", type = transformationsMap[[name]], pars = pars)
}

applyTransformation <- function(transformation, x)
{
  pars <- transformation@pars
  
  switch(transformations[transformation@type],
         identity = x,
         exp      = exp(x),
         power    = sign(x) * abs(x)^pars,
         sigmoid  =
           if (pars[2L] < 0)
             pnorm(-pars[1L] * (x - pars[2L]))
           else
             pnorm( pars[1L] * (x - pars[2L]))
         )
}

getPrefixForTransformation <- function(transformation, symbolic, digits = getOption("digits"))
{
  switch(transformations[transformation@type],
         identity = if (symbolic) "" else "(",
         exp      = "exp(",
         power    = "power(",
         sigmoid  = paste0("pnorm(",
           if (symbolic)
             "a"
           else 
             paste0(if (transformation@pars[2L] < 0) "-" else "", signif(transformation@pars[1L], digits)),
           " * (")
  )
}

getSuffixForTransformation <- function(transformation, symbolic, digits = getOption("digits"))
{
  switch(transformations[transformation@type],
         identity = if (symbolic) "" else ")",
         exp      = ")",
         power    = if (symbolic) ", b)" else paste0(", ", signif(transformation@pars, digits), ")"),
         sigmoid  = paste0(" - ", if (symbolic) "b" else signif(transformation@pars[2L], digits), "))")
  )
}

