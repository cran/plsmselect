###########
# WRAPPER #
###########

#' Fitting a gamlasso model
#'
#' This function will fit a gamlasso model with the given penalties. For some
#' special cases using \code{\link[mgcv]{gam}} or \code{\link[glmnet]{glmnet}}
#' might be more efficient and/or flexible
#'
#' @param formula A formula describing the model to be fitted
#' @param response The name of the response variable. Could be two variables
#'   in case of a general binomial fit (see details below)
#' @param linear.terms The names of the variables to be used as linear predictors
#' @param smooth.terms The names of the variables to be used as smoothers
#' @param data The data with which to fit the model
#' @param family The family describing the error distribution and link function
#'   to be used in the model. A character string which can only be
#'   \code{"gaussian"} (default), \code{"binomial"}, \code{"poisson"} or
#'   \code{"cox"}. For \code{family = "binomial"}, \code{response} can be
#'   a vector of two and for \code{family="cox"}, \code{weights} must
#'   be provided (see details below).
#' @param linear.penalty The penalty used on the linear predictors. A character
#'   string which can be \code{"none"} (default), \code{"l1"} or \code{"l2"}. If
#'   \code{"l1"} is used then we use the gam and lasso loop. Otherwise only a
#'   gam model is fitted (with penalities on parametric terms if
#'   \code{linear.penalty = "l2"} ).
#' @param smooth.penalty The penalty used on the smoothers. A character
#'   string which can be \code{"l1"} or \code{"l2"} (default). \code{"l2"} refers
#'   to the inherent second order penalty smoothers have for controlling their
#'   shape, so \code{"none"} is not an option. For \code{"l1"} basis is specified
#'   by \code{bs='ts'}, else \code{bs='tp'} is used. (see \code{\link[mgcv]{gam}}
#'   for details on basis types)
#' @param num.knots Number of knots for each smoothers. Can be a single integer
#'   (recycled for each smoother variable) or a vector of integers the same length
#'   as the number of smoothers.
#' @param offset The name of the offset variable. \code{NULL} (default) if not provided
#' @param weights The name of the weights variable. \code{NULL} (default) if not
#'   provided. See details below.
#' @param interactions logical. Should interactions be included as covariates.
#'   If \code{TRUE} then the smoothers are fitted with \code{\link[mgcv]{ti}}
#'   instead of \code{\link[mgcv]{s}} so that the added effects of the interactions
#'   can be quantified separately.
#' @param seed The random seed can be specified for reproducibility. This is used
#'   for fitting the gam and lasso models, or fixed before each loop of gamlasso.
#' @param num.iter Number of iterations for the gamlasso loop
#' @param tolerance Tolerance for covergence of the gamlasso loop
#' @param prompts logical. Should \code{gamlassoChecks} provide interactive
#'   user prompts for corrective action when needed.
#' @param verbose logical. Should there be "progress reports" printed to the
#'   console while fitting the model.
#' @param ... Additional arguments
#'
#' @usage NULL
#'
#' @details \code{gamlasso} allows for specifying models in two ways:
#'   1) with the the formula approach, and 2) with the term specification approach.
#'
#'   The formula approach is appropriate for when the user wants an L1-penalty on the
#'   linear terms of the model, in which case the user is required to specify the linear terms
#'   in a model matrix named "X" appended to the input data frame. A typical formula specification
#'   would be "\code{y ~ X + s(z) + ...}" where "\code{X}" corresponds to the model-matrix of
#'   linear terms subject to an L1-penalty, while everything to the right of "\code{X}" is
#'   considered part of the gam formula (i.e. all smooth terms). In light of the above formula,
#'   gamlasso iterates (until convergence) between the following two lines of pseudo code:
#'
#'   \itemize{
#'     \item \code{model.cv.glmnet <- cv.glmnet(y=y, x=X, offset="model.gam fitted values")}
#'     \item \code{model.gam <- gam(y ~ s(z) + ..., offset="model.cv.glmnet fitted values")}
#'   }
#'
#'   The term specification approach can fit the same type of models as the formula approach
#'   (i.e. models with L1-penalty on the linear terms). However, it is more flexible in terms
#'   of penalty-structure and can be useful if the user has big data sets with lots of variables
#'   making the formula specification cumbersome. In the term specification approach
#'   the user simply specifies the names of the data columns corresponding to the
#'   \code{response}, \code{linear.terms} and \code{smooth.terms} and then specifies
#'   whether to put a \code{linear.penalty="l1"}, \code{"l2"} or \code{"none"}
#'   (on \code{linear.terms}) and whether to put a \code{smooth.penalty="l1"} or
#'   \code{"l2"} (on \code{smooth.terms}).
#'
#'   While fitting a binomial model for binary responses (0/1) include the response
#'   variable before "~" if using the formula approach or when using the term-
#'   specification approach the \code{response} argument will be a single variable name.
#'   In general if the responses are success/failure counts then the formula should
#'   start with something similar to \code{cbind(success,failure) ~ ...} and for
#'   using the term-specification approach the \code{response} argument should be a
#'   vector of length two giving the success and failure variable names.
#'
#'   If \code{family="cox"} then the \code{weights} argument must be provided
#'   and should correspond to a status variable (1-censor). For other models
#'   it should correspond to a custom weights variables to be used for the
#'   weighted log-likelihood, for example the total counts for fitting a
#'   binomial model. (weights for families other than "cox" currently not
#'   implemented)
#'
#'   Both the formula and term-specification approaches can fit interaction models as
#'   well. There are three kinds of interactions - those between two linear predictors,
#'   between two smooth predictors and between linear and smooth predictors. For the
#'   formula approach the first type of interaction must be included as additional
#'   columns in the "\code{X}" matrix and the other two types must be mentioned in the
#'   smooth terms part of the formula. For the term-specification approach the argument
#'   \code{interaction} must be \code{TRUE} in which case all the pairwise
#'   interactions are used as predictors and variable selection is done on all of them.
#'
#' @note The default values of \code{num.iter} and \code{tolerance} are
#'   essentially arbitrary. Also for each step when we check for convergence
#'   between the new and old predictions by the gam and lasso predictions,
#'   we use the following distance metric
#'   \deqn{ d(x,y) = \frac{1}{length(x)} \sum_{i=1}^{length(x)} (x_i - y_i)^2 }
#'
#' @return If the arguments fail the basic checking by \code{gamlassoChecks}
#'   then returns \code{NULL}. Else the function calls \code{gamlassoFit} which
#'   returns a list of two models, \code{gam} and \code{cv.glmnet}.
#'   Either of these could be \code{NULL} but if both are non-null then
#'   \code{convergence}, a matrix of values determining the convergence
#'   of the gamlasso loop is also returned.
#'   \code{gamlassoFit} also returns \code{inherit}, a list of select
#'   arguments used to fit the \code{gamlasso} model and some more values needed
#'   for prediction.
#'
#' @export
#'
#' @seealso \code{\link[mgcv]{gam}}, \code{\link[glmnet]{glmnet}}
#'
#' @examples
#' library(plsmselect)
#'
#' data(simData)
#'
#' ## Fit gaussian gamlasso model using the formula approach:
#' ## (L1-penalty both on model matrix (X) and smooth terms (bs="ts"))
#' simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]
#'
#' gfit = gamlasso(Yg ~ X +
#'                    s(z1, k=5, bs="ts") +
#'                    s(z2, k=5, bs="ts") +
#'                    s(z3, k=5, bs="ts") +
#'                    s(z4, k=5, bs="ts"),
#'                    data = simData,
#'                    seed=1)
#'
#' \donttest{## Equivalently with term specification approach:
#' gfit = gamlasso(response="Yg",
#'                   linear.terms=paste0("x",1:10),
#'                   smooth.terms=paste0("z",1:4),
#'                   data=simData,
#'                   linear.penalty = "l1",
#'                   smooth.penalty = "l1",
#'                   num.knots = 5,
#'                   seed=1)
#' }
#' ## The two main components of gfit are
#' ## gfit$cv.glmnet (LASSO component) and gfit$gam (GAM components):
#'
#' ## Extract lasso estimates of linear terms:
#' coef(gfit$cv.glmnet, s="lambda.min")
#'
#' ## Plot the estimates of the smooth effects:
#' plot(gfit$gam, pages=1)
#'
#' # See ?summary.gamlasso for an example fitting a binomial response model
#' # See ?predict.gamlasso for an example fitting a poisson response model
#' # See ?cumbasehaz for an example fitting a survival response model
gamlasso = function(...) {
  UseMethod("gamlasso")
}

#' @rdname gamlasso
#'
#' @export
gamlasso.formula = function(formula, data, family = "gaussian",
                            linear.penalty = "l1", smooth.penalty = "l2", num.knots = 5,
                            offset = NULL, weights = NULL, interactions = F,
                            seed = .Random.seed[1], num.iter = 100, tolerance = 1e-4, ...)
{
  linear.penalty = switch(linear.penalty,
                          "none" = 0,
                          "l1" = 1,
                          "l2" = 2)
  smooth.penalty = switch(smooth.penalty,
                          "l1" = 1,
                          "l2" = 2)

  ## There are no checks currently

  output = gamlassoFit(data = data, formula = formula,
                       family = family, linear.penalty = linear.penalty,
                       smooth.penalty = smooth.penalty, offset.name = offset,
                       weights.name = weights, num.knots = num.knots,
                       num.iter = num.iter, interactions = interactions,
                       tolerance = tolerance, seed = seed, ...)

  return(output)

}

#' @rdname gamlasso
#'
#' @export
gamlasso.default = function(response, linear.terms, smooth.terms, data, family = "gaussian",
                            linear.penalty = "l1", smooth.penalty = "l2", num.knots = 5,
                            offset = NULL, weights = NULL, interactions = F,
                            seed = .Random.seed[1], num.iter = 100, tolerance = 1e-4,
                            prompts = F, verbose = T, ...)
{
  # M = data.frame(alph = c("none","l1","l2"), num = 0:2)
  # linear.penalty = M$num[M$alph == linear.penalty]

  linear.penalty = switch(linear.penalty,
                          "none" = 0,
                          "l1" = 1,
                          "l2" = 2)
  smooth.penalty = switch(smooth.penalty,
                          "l1" = 1,
                          "l2" = 2)

  ### CHECKS AND BALANCES
  allcheck = FALSE
  check_clean = gamlassoChecks(data = data, response.name = response,
                                linear.name = linear.terms, smooth.name = smooth.terms,
                                family = family, linear.penalty = linear.penalty,
                                smooth.penalty = smooth.penalty, offset.name = offset,
                                weights.name = weights, num.knots = num.knots, num.iter,
                                tolerance, seed, prompts = prompts)
  allcheck = check_clean$check
  if(!allcheck)
  {
    if(verbose) message("\n\nSomething went wrong. Please try again. Exiting ... \n")
    return(NULL)
  }

  if(verbose) message("\nPreliminary checks passed. ")
  ## DO STUFF ACCORDING TO THE RESULTS OF THE PROMPTS

  train.data = check_clean$clean$train.data
  linear.terms = check_clean$clean$linear.name
  smooth.terms = check_clean$clean$smooth.name
  linear.penalty = check_clean$clean$linear.penalty
  smooth.penalty = check_clean$clean$smooth.penalty
  num.knots = check_clean$clean$num.knots

  output = gamlassoFit(data = train.data, response.name = response,
                        linear.name = linear.terms, smooth.name = smooth.terms,
                        family = family, linear.penalty = linear.penalty,
                        smooth.penalty = smooth.penalty, offset.name = offset,
                        weights.name = weights, num.knots = num.knots,
                        num.iter = num.iter, interactions = interactions,
                        tolerance = tolerance, seed = seed, verbose = verbose)
  return(output)

}

###################
## MAIN FUNCTION ##
###################

#' The function fitting a gamlasso model
#'
#' This function is the workhorse for fitting a gamlasso model. Not recommended
#' to call directly. It is slightly more efficient than \code{gamlasso.default} since
#' it doesn't perform any quality checks. Only use if the data has been cleaned
#' and no errors are expected to occur.
#'
#' @inheritParams gamlassoChecks
#' @inheritParams gamlasso.formula
#' @param interactions logical. Should interactions be included.
#'
#' @return See \code{\link{gamlasso}}
#' @export
#' @importFrom dplyr "%>%"
#' @importFrom stats as.formula model.matrix
#'
#' @examples ## Not recommended to use directly. Please see examples of gamlasso
gamlassoFit = function(data, formula = NULL, response.name = NULL,
                        linear.name = NULL, smooth.name = NULL,
                        family = "gaussian", linear.penalty = 0, smooth.penalty = 2,
                        offset.name = NULL, weights.name = NULL, num.knots = 5,
                        num.iter = 100, interactions = F, tolerance = 1e-4,
                        seed = .Random.seed[1], verbose = TRUE)
{
  formulae = formula_setup(formula, response.name, linear.name, smooth.name, family,
                           smooth.penalty, num.knots, offset.name, interactions)
  families = find_family(family)
  if(is.null(response.name)) response.name = detect_one_two(formulae$response)

  train.offset = train.weights = NULL
  if(!is.null(offset.name)) train.offset = data %>% dplyr::pull(offset.name)
  if(!is.null(weights.name)) train.weights = data %>% dplyr::pull(weights.name)

  onlysmooth = is.null(formulae$linear)
  onlylinear = is.null(formulae$smooth)

  model_gam = model_lasso = model_convmat = NULL
  if(onlysmooth)
  {
    eqn = as.formula(paste(formulae$response, "~", formulae$smooth, formulae$offset))
    set.seed(seed)
    if(verbose) message("Fitting a gam model with only smoothers ... \n")
    model_gam = mgcv::gam(eqn, family = families$gam_family, data = data,
                          weights = train.weights)
  }
  else
  {
    if(linear.penalty!=1)
    {
      if(onlylinear)
      {
        if(verbose) message("Fitting a gam model with only parametric terms ... \n")
        eqn = as.formula(paste(formulae$response, "~", formulae$linear, formulae$offset))
      }
      else
      {
        if(verbose) message("Fitting a gam model with parametric and smoothing terms ... \n")
        eqn = as.formula(paste(formulae$response, "~", formulae$linear,
                               "+", formulae$smooth, formulae$offset))
      }

      PP = NULL
      if(linear.penalty == 2)
      {
        xmat = model.matrix(as.formula(paste("~", formulae$linear)),
                            data = data)[,-1]
        dimP = dim(xmat)[2]
        PP = list(X=list(rank=dimP,0.5*diag(dimP)))
      }

      set.seed(seed)
      model_gam = mgcv::gam(eqn, family = families$gam_family, data = data,
                            weights = train.weights, paraPen = PP)
    }
    else
    {
      if(onlylinear)
      {
        xmat = model.matrix(as.formula(paste("~", formulae$linear)),
                            data = data)
        y = lasso_response(response.name, families$lasso_family,
                           data, train.weights)
        set.seed(seed)

        if(verbose) message("Fitting penalised glm ... \n")
        cv.lasso.model = glmnet::cv.glmnet(x=xmat, y=y, offset=train.offset,
                                           family=families$lasso_family, maxit=1000)
        lasso.model = glmnet::glmnet(x=xmat, y=y, offset=train.offset,
                                     family=families$lasso_family, maxit=1000)
        model_lasso = list(cv.model = cv.lasso.model, model = lasso.model)
      }
      else
      {
        if(verbose) message("Fitting a lasso and gam loop ... \n")
        models = lasso_gam_loop(data, response.name, families,
                                formulae, num.iter, tolerance,
                                offset.name, train.weights, seed)
        model_lasso = models$lasso
        model_gam = models$gam
        model_convmat = models$convmat
        if(models$convergence == FALSE)
        {
          warning(paste("Maximum number of iterations reached for",
                        "the gamlasso loop without convergence.",
                        "Try increasing the num.iter argument.") )
          # AND PRINT 7 NUMBERS SHOWING HOW CLOSE STUFF WAS TO THE END
          # M = models$convmat
          # M = M[,ncol(M) - (1:0)]
          # print(c(M[1:4,2],M[5,],M[6,2])) # add description
        }
      }
    }
  }

  names = list(response = response.name, linear = linear.name,
               smooth = smooth.name, offset = offset.name,
               weights = weights.name)
  inherit = list(train.data = data, fit.names = names,
                 family = family, formulae = formulae)

  output = list(gam = model_gam, cv.glmnet = model_lasso$cv.model,
                convergence = model_convmat, inherit = inherit)
  class(output) = "gamlasso"
  return(output)
}
