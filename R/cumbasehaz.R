###################
# CUMULATIVE BASE #
# HAZARD FUNCTION #
###################

#' Cumulative Baseline Hazard of a gamlasso object
#'
#' This is only used when with family="cox"
#'
#' @param object fitted model object of the class \code{gamlasso} as
#'   produced by \code{gamlasso}
#'
#' @return This function returns the cumulative baseline hazard function of
#'   a \code{gamlasso} object if fitted using \code{family = "cox"}. More
#'   specifically, cumbasehaz(object) is the cumulative baseline hazard function
#'   corresponding to the linear predictor predict(object).
#' @export
#'
#' @seealso \code{\link{gamlasso}}
#'
#' @examples
#' library(plsmselect)
#'
#' data(simData)
#'
#' ## Fit Cox gamlasso model using the formula approach:
#' ## (L1-penalty both on X terms and smooth terms (bs="ts"))
#' simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]
#'
#' cfit = gamlasso(time ~ X +
#'                   s(z1, bs="ts", k=5) +
#'                   s(z2, bs="ts", k=5) +
#'                   s(z3, bs="ts", k=5) +
#'                   s(z4, bs="ts", k=5),
#'                 data = simData,
#'                 family = "cox",
#'                 weights="status",
#'                 seed=1)
#'
#' ## Obtain and plot predicted cumulative baseline hazard:
#' H0.pred <- cumbasehaz(cfit)
#'
#' time.seq <- seq(0, 60, by=1)
#' plot(time.seq, H0.pred(time.seq), type="l", ylab="Predicted Cumulative Baseline Hazard")
#'
#' ## Obtain predicted survial probabilities at month 1 and 2 (days 30 & 60):
#'
#' lp <- predict(cfit) # estimated linear predictor
#'
#' S.pred <- cbind(exp(-H0.pred(30)*exp(lp)), exp(-H0.pred(60)*exp(lp)))
#'
#' ## Obtain predicted survival at month 1 and 2 directly:
#' S.pred2 <- predict(cfit, type="response", new.event.times=c(30,60))
#'
#' ## Confirm that the two arrived at the same values:
#' all.equal(S.pred, S.pred2)
#'
#' # See ?gamlasso for an example fitting a gaussian response model
#' # See ?summary.gamlasso for an example fitting a binomial response model
#' # See ?predict.gamlasso for an example fitting a poisson response model
cumbasehaz = function(object)
{
  if(class(object) != "gamlasso")
    stop("Please enter a gamlasso object to print")

  if(object$inherit$family != "cox")
    stop("Cumulative base hazard function only makes sense for cox models")

  lp <- predict.gamlasso(object, newdata=object$inherit$train.data, type="link")
  event.time <- object$inherit$train.data %>%
    dplyr::pull(object$inherit$fit.names$response)
  status <- object$inherit$train.data %>%
    dplyr::pull(object$inherit$fit.names$weights)

  L0 = cbh(lp, event.time, status)
  return(L0)
}

## HELPER FUNCTION FOR cumbasehaz WHICH ACTUALLY MAKES THE FUNCTION

#' Internal Function
#'
#' Undocumented function. Do not use directly
#'
#' @param lp The linear predictor to be used as offset
#' @param event.time The event times
#' @param status Status indicating the complement of censoring
#'
#' @importFrom stats approxfun
cbh = function(lp, event.time, status)
{
  cox.model <- survival::coxph(survival::Surv(event.time, status) ~ offset(lp))
  sf <- survival::survfit(cox.model)
  H0 <- approxfun(sf$time, sf$cumhaz*exp( -mean(lp) ) )
  return(H0)
}



