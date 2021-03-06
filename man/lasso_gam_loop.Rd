% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{lasso_gam_loop}
\alias{lasso_gam_loop}
\title{Internal Function}
\usage{
lasso_gam_loop(
  data,
  response.name,
  families,
  formulae,
  num.iter,
  tolerance,
  offset.name,
  weights,
  seed
)
}
\arguments{
\item{data}{The data with value for all the linear and smooth predictors}

\item{response.name}{The name of the response variable. Vector of two if
\code{family = "binomial"}}

\item{families}{List of two families as returned by \code{find_family}}

\item{formulae}{List of formulae as returned by \code{formula_setup}}

\item{num.iter}{Number of iterations for the gamlasso loop}

\item{tolerance}{Tolerance for covergence of the gamlasso loop}

\item{offset.name}{The name of the offset variable. \code{NULL} (default) if not provided}

\item{weights}{Vector with values of the weights variable if it exists.
\code{NULL} otherwise.}

\item{seed}{The random seed can be specified for reproducibility. This is used
for fitting the gam and lasso models, or fixed before each loop of gamlasso.}
}
\description{
Undocumented function. Do not use directly
}
