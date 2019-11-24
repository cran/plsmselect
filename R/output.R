#########
# PRINT #
#########

#' Print a gamlasso object
#'
#' The default print method for a \code{gamlasso} object
#'
#' @param x fitted model object of the class \code{gamlasso} as
#'   produced by \code{gamlasso}
#' @param ... Other arguments
#'
#' @details Outputs a list of two. \code{lasso} prints the lasso
#'   model (the same output as \code{print(object$cv.glmnet$glmnet.fit)}) if
#'   it is non-null and \code{gam} prints the gam model (the same output
#'   as \code{print(object$gam)}) if it is non-null.
#' @export
#'
#' @seealso \code{\link{gamlasso}}, \code{\link{summary.gamlasso}},
#'   \code{\link[mgcv]{print.gam}}, \code{\link[glmnet]{print.glmnet}}.
#'
#' @examples ## Please see the examples in ?gamlasso
print.gamlasso = function(x, ...)
{
  if(class(x) != "gamlasso")
    stop("Please enter a gamlasso object to print")

  print_lasso = print_gam = NULL
  if(!is.null(x$cv.glmnet))
  {
    print_lasso = print(x$cv.glmnet$glmnet.fit)
  }

  if(!is.null(x$gam))
  {
    print_gam = mgcv::print.gam(x$gam)
  }

  output = list(lasso = print_lasso, gam = print_gam)
  return(output)
}

###########
# SUMMARY #
###########

#' Summary for a gamlasso fit
#'
#' Default sumary method for a \code{gamlasso} object
#'
#' @param object fitted model object of the class \code{gamlasso} as
#'   produced by \code{gamlasso}
#' @param s Value of the lasso penalty parameter \code{lambda} at which
#'   predictions are required. Default is \code{"lambda.min"} but alternatively
#'   \code{"lambda.1se"} can be used.
#' @param ... Other arguments
#'
#' @details Outputs a list of two. \code{gam} prints a summary of the
#'   gam model (the same output as \code{summary(object$gam)}) if it
#'   is non-null. Objects of the class \code{cv.glmnet} do not have a
#'   default summary method, so the list item \code{lasso} produces the
#'   coefficients of the cross-vaidated lasso fit corresponding to the
#'   lowest value of the \eqn{\lambda} used ( the same output as
#'   \code{coef(object$cv.glmnet, s = "lambda.min")} if it is non-null).
#' @export
#'
#' @seealso \code{\link{gamlasso}}, \code{\link[mgcv]{summary.gam}},
#'   \code{\link[glmnet]{coef.cv.glmnet}}.
#'
#' @examples
#' library(plsmselect)
#'
#' data(simData)
#'
#' ## Fit binomial gamlasso model using the term specification
#' ## approach with binomial counts response
#' ## (L2-penalty on linear terms & L1-penalty on smooth terms)
#' bfit = gamlasso(c("success","failure"),
#'                 linear.terms=paste0("x",1:10),
#'                 smooth.terms=paste0("z",1:4),
#'                 data=simData,
#'                 family = "binomial",
#'                 linear.penalty = "l2",
#'                 smooth.penalty = "l1",
#'                 num.knots = 5,
#'                 seed=1)
#'
#' ## Since the above model has linear.penalty = "l2" it is
#' ## a pure GAM model (i.e. no LASSO component):
#' bfit$cv.glmnet
#'
#' ## Summary of model (here essentially the same as summary(bfit$gam)
#' ## because there is no LASSO component, i.e. linear.penalty="l2")
#' summary(bfit)
#'
#' \donttest{## We could use the formula approach below to fit the same model as above:
#' simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]
#' bfit = gamlasso(cbind(success,failure) ~ X + s(z1, bs="ts") +
#'                  s(z2, bs="ts") + s(z3, bs="ts") + s(z4, bs="ts"),
#'                 data = simData,
#'                 family = "binomial",
#'                 linear.penalty = "l2",
#'                 smooth.penalty = "l1",
#'                 seed=1)
#'
#' ## For a binary responses we only need one response variable in the formula
#' bfit2 = gamlasso(Yb ~ X + s(z1, bs="ts") + s(z2, bs="ts") + s(z3, bs="ts") + s(z4, bs="ts"),
#'                   data = simData,
#'                   family = "binomial",
#'                   seed=1)
#' }
#' # See ?gamlasso for an example fitting a gaussian response model
#' # See ?predict.gamlasso for an example fitting a poisson response model
#' # See ?cumbasehaz for an example fitting a survival response model
summary.gamlasso = function(object, s = "lambda.min", ...)
{
  if(class(object) != "gamlasso")
    stop("Please enter a gamlasso object to summarise")

  summary_lasso = summary_gam = NULL
  if(!is.null(object$cv.glmnet))
  {
    summary_lasso = stats::coef(object$cv.glmnet, s = s)
  }

  if(!is.null(object$gam))
  {
    summary_gam = mgcv::summary.gam(object$gam)
  }

  output = list(lasso = summary_lasso, gam = summary_gam)
  return(output)

}

##############
# PREDICTION #
##############

## type CAN BE "link" OR "response"

#' Prediction from a fitted gamlasso model
#'
#' Takes a fitted \code{gamlasso} object produced by \code{gamlasso} and
#'   returns predictions given a new set of values of the linear and smooth
#'   variables.
#'
#' @param object fitted model object of the class \code{gamlasso} as
#'   produced by \code{gamlasso}
#' @param newdata A data frame with the values of the linear and smooth
#'   variables for which predictions are to be made. If not provided then
#'   predictions corresponding to the original data used to fit \code{object}
#'   is returned. If provided then the variable names (column names) should
#'   match with the variable names used to fit \code{object}: the code throws
#'   an error if not.
#' @param type When this has the value \code{"link"} (default) then the
#'   linear predictor (with offset added if needed) is returned. When
#'   \code{type = "response"} predictions on the response scale is returned,
#'   depending on the family used while fitting \code{object}.
#' @param s Value of the lasso penalty parameter \code{lambda} at which
#'   predictions are required. Default is \code{"lambda.min"} but alternatively
#'   \code{"lambda.1se"} can be used.
#' @param new.event.times A vector of new event times to be used for predicting
#'   survival times when \code{type = "response"} for a gamlasso object fitted
#'   with \code{family = "cox"}
#' @param ... Other arguments
#'
#' @return Returns a vector of the same length as \code{nrow(newdata)} with
#'   the values of the linear predictor or on the response scale depending
#'   on \code{type}. For \code{type = "link"} the value is simply the elementwise
#'   sum of the predictions from the gam and lasso models in \code{object}.
#'   For \code{type = "response"} the values are on the response scale, for
#'   example exponential of the linear response is returned if
#'   \code{object$inherit$family = "poisson"}
#'
#' @export
#' @importFrom stats as.formula model.matrix
#'
#' @details Lasso models do not have standard errors so \code{predict.gamlasso}
#'   does not provide them either. The standard errors for the gam part of the
#'   model can be accesed by using \code{mgcv::predict.gam} with suitable options.
#'   Offsets are always included in the prediction if present in the original call
#'   to \code{gamlasso}. Also if \code{type} is anything other than \code{"link"}
#'   or \code{"response"} then the function throws an error.
#'
#' @seealso \code{\link{gamlasso}}, \code{\link[mgcv]{predict.gam}},
#'   \code{\link[glmnet]{predict.glmnet}}.
#'
#' @importFrom dplyr pull
#' @examples
#' library(plsmselect)
#'
#' data(simData)
#'
#' ## Fit poisson gamlasso model using the term specification approach:
#' ## (L2-penalty on linear terms & L2-penalty on smooth terms)
#' pfit = gamlasso(response="Yp",
#'                 linear.terms=paste0("x",1:10),
#'                 smooth.terms=paste0("z",1:4),
#'                 data=simData,
#'                 linear.penalty = "l2",
#'                 smooth.penalty = "l2",
#'                 family="poisson",
#'                 num.knots = 5,
#'                 seed=1)
#'
#' ## fitted values (of linear predictor):
#' fitted.values <- predict(pfit)
#'
#' ## predicted values on response scale:
#' pred.response <- predict(pfit, type="response", newdata=simData)
#'
#' \donttest{## For same model as above, but with L1-penalty on linear terms
#' ## i.e. L1-penalty on the model matrix (X) we can use formula approach:
#' simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]
#'
#' pfit = gamlasso(Yp ~ X +
#'                    s(z1, k=5) + # L2-penalty (bs="tp") is default (see ?mgcv::s)
#'                    s(z2, k=5) +
#'                    s(z3, k=5) +
#'                    s(z4, k=5),
#'                  family="poisson",
#'                  data = simData,
#'                  seed=1)
#' }
#' # See ?gamlasso for an example fitting a gaussian response model
#' # See ?summary.gamlasso for an example fitting a binomial response model
#' # See ?cumbasehaz for an example fitting a survival response model
predict.gamlasso = function(object, newdata = NULL, type = "link", s = "lambda.min",
                            new.event.times = NULL, ...)
{
  if(class(object) != "gamlasso")
    stop("Please enter a gamlasso object to predict")

  test.data = newdata
  if(is.null(newdata))
    test.data = object$inherit$train.data
  else
  {
    ## CHECKS ON VARIABLE NAMES IN test.data
    names = unlist(object$inherit$fit.names)
    if( !all(names %in% colnames(test.data) ) )
    {
      errmsg = paste("The variables",
                     paste(setdiff(names, colnames(test.data)), collapse = ", "),
                     "are not in newdata")
      stop(errmsg)
    }
    ## NOT ALL type IS IMPLEMENTED
    if(!(type %in% c("link", "response")))
    {
      stop("'type' must be either 'link' or 'response'.")
    }
  }

  test.offset = 0
  if(!is.null(object$inherit$fit.names$offset))
    test.offset = test.data %>% dplyr::pull(object$inherit$fit.names$offset)

  Xbeta.linear = Xbeta.smooth = 0

  if(is.null(object$gam))
  {
    test.xmat = model.matrix(as.formula(paste("~", object$inherit$formulae$linear,
                                              object$inherit$formulae$offset)),
                             data = test.data)
    Xbeta.linear = stats::predict(object$cv.glmnet,
                                  newx = test.xmat, s = s,
                                  type = "link", newoffset = test.offset)
  }
  else
  {
    Xbeta.smooth = mgcv::predict.gam(object$gam,
                                     newdata = test.data, type = "link")
    if(!is.null(object$cv.glmnet))
    {
      test.xmat = model.matrix(as.formula(paste("~", object$inherit$formulae$linear)),
                               data = test.data)[,-1]

      ## Don't put newoffset to be 0 because we would need to add it to
      ## Xbeta.final. But in case there is ony gam fitted then that has
      ## offset in the formula so it will be included in the predict
      ## already leading to double offset adding
      Xbeta.linear = stats::predict(object$cv.glmnet,
                                    newx = test.xmat, s = s,
                                    type = "link", newoffset = test.offset)
    }
  }

  Xbeta.final = as.numeric(Xbeta.linear) + as.numeric(Xbeta.smooth)
  if(type == "link")
    return(Xbeta.final)
  else if(type == "response")
  {
    ## IF FAMILY IS cox THEN RETURN SURVIVAL PROBABLITIES
    if(object$inherit$family == "cox")
    {
      lp = Xbeta.final
      if(!is.null(newdata))
        lp = predict.gamlasso(object, newdata = NULL, type = "link")
      event.time <- object$inherit$train.data %>%
        pull(object$inherit$fit.names$response)
      status <- object$inherit$train.data %>%
        pull(object$inherit$fit.names$weights)
      L0 = cbh(lp, event.time, status)

      S_mat = -outer(exp(Xbeta.final), L0(new.event.times)) %>% exp()
      return(S_mat)
    }

    ## ELSE DO STUFF BASED ON INVERSE LINK FUNCTION OF object$inherit$family
    linv = find_family(object$inherit$family)$gam_family$linkinv
    return(linv(Xbeta.final))
  }

}


###############
# DIAGNOSTICS #
###############

## Will be omitted for now. Might bring back later after careful testing

## #' Internal Function
## #'
## #' Undocumented function. Do not use directly
## #'
## #' @param M the matrix in the \code{Details} section of
## #'   \code{\link{gamlassoConvergence}}
## #'
## #' @importFrom graphics legend par plot points
## plotconvmat = function(M)
## {
##   old_par = par(mfrow = c(2,2))
##
##   plot(M[1,-1], xlab = "iteration (i)", type = 'l',
##        ylab = expression(paste("|| ", hat(beta)[i] - hat(beta)[i-1], " || for lasso")),
##        main = "Distance between succesive iterates \n of the lasso coefficients")
##
##   plot(M[3,-1], xlab = "iteration (i)", type = 'l',
##        ylab = expression(paste("|| ", hat(beta)[i] - hat(beta)[i-1], " || for gam")),
##        main = "Distance between succesive iterates \n of the gam coefficients")
##
##   plot(M[2,-1], xlab = "iteration (i)", type = 'l', ylim = c(-1,1), col = "red",
##        ylab = expression(paste("Corr(", hat(beta)[i], ",", hat(beta)[i-1], ")")),
##        main = "Correlation between succesive iterates")
##   points(M[4,-1], type = 'l', col = "blue")
##   legend("bottomright", legend = c("lasso", "gam"),
##          col = c("red", "blue"), lty = 1)
##
##   plot(0:(ncol(M)-1), M[5,], xlab = "iteration (i)", type = 'b', col = "red",
##        ylab = "No. of zeros", main = "Variable selection in the lasso iterates")
##   points(1:(ncol(M)-1), M[6,-1], type = 'b', col = "blue")
##   legend("bottomright", legend = c("individual", "common \nwith prev"),
##          col = c("red", "blue"), lty = 1, pch = 1, pt.bg = 'white')
##
##   suppressWarnings(par(old_par))
## }
##
##
## #' Gamlasso convergence diagnostics
## #'
## #' If the gamlasso object has both gam and lasso objects then this function
## #' checks the extent of convergence between the two models.
## #'
## #' @param object fitted model object of the class \code{gamlasso} as
## #'   produced by \code{gamlasso}
## #'
## #' @return The function produces diagnostics plots of convergence. It also
## #'   silently returns the matrix used in the plots, which can be used for
## #'   customised plotting.
## #'
## #' @details \code{gamlasso} returns a matrix of values determining the extent
## #'   of convergence of the gamlasso loop. This matrix has 6 rows with number of
## #'   columns equal to the number of iterations of the gam and lasso loop. The
## #'   first and third rows of the matrix has the distances between coefficients
## #'   of the lasso and gam models respectively from succesive iterations of the
## #'   loop ( the distance is the same as described in the \code{Details} section
## #'   in \code{\link{gamlasso}} ). The second and fourth rows contain the correlation
## #'   between the above discussed coefficients. The fifth row has the number of zero
## #'   coefficients of the lasso models and the sixth row has the number of common
## #'   zeros between succesive lasso coefficients.
## #'
## #' @seealso \code{\link{gamlasso}}, \code{\link[mgcv]{plot.gam}},
## #'   \code{\link[glmnet]{plot.glmnet}}.
## #'
## #' @examples ## In development
## gamlassoConvergence = function(object)
## {
##   # Plot the matrix object$convergence properly if exists.
##   if(is.null(object$convergence))
##   {
##     cat(paste("\nThere was no gam and lasso loop."))
##     return(invisible(NULL))
##   }
##   else
##   {
##     cat(paste("\nPlots for convergence of the gam and lasso loop"))
##
##     M = object$convergence
##     plotconvmat(M)
##     return(invisible(M))
##   }
##
##   # For the rest give a message to use plot.gam or plot.glmnet
##   plot_gam = "mgcv::plot.gam(.$gam)"
##   plot_lasso = "glmnet::plot.cv.glmnet(.$cv.glmnet)"
##
##   if(is.null(object$cv.glmnet))
##   {
##     cat(paste("\nUse", plot_gam, "for more plots."))
##   }
##   else
##   {
##     if(is.null(object$gam))
##     {
##       cat(paste("\nUse", plot_lasso, "for more plots."))
##     }
##     else
##     {
##       cat(paste("\nUse", plot_gam, "and", plot_lasso, "for more plots."))
##     }
##   }
## }
