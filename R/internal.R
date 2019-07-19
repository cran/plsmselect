## SECONDARY FUNCTIONS

#' Internal Function
#'
#' Undocumented function. Do not use directly
#'
#' @param fam Family in character form
find_family = function(fam)
{
  lasso_family = fam
  if(fam == "cox") gam_family = mgcv::cox.ph()
  else gam_family = eval(parse(text = paste(fam, "()", sep = "")))
  return(list(gam_family = gam_family, lasso_family = lasso_family))
}

#' Internal Function
#'
#' Undocumented function. Do not use directly
#'
#' @inheritParams gamlassoChecks
#' @param formula A formula to be parsed
#' @param interactions logical. Should interactions be included.
#'
#' @importFrom dplyr "%>%"
formula_setup = function(formula = NULL, response.name = NULL, linear.name = NULL,
                         smooth.name = NULL, family = NULL, smooth.penalty = NULL,
                         num.knots = NULL, offset.name = NULL, interactions = F)
{
  if(is.null(formula))
  {
    if(length(response.name)==1)
      r.formula = response.name
    else
    {
      if(family=="binomial")
        r.formula = paste("cbind(", paste(response.name, collapse = ","), ")", sep = "")
      else
        stop("Multiple responses allowed only for binomial models.")
    }

    if(is.null(linear.name)) {
      l.formula = NULL
    } else {
      l.formula = paste(linear.name, collapse = "+")
    }

    if(is.null(smooth.name)) {
      s.formula = NULL
    }
    else
    {
      if(smooth.penalty == 1)
        s.formula = paste(paste("s(",smooth.name,",k=",num.knots,",bs='ts')",sep=""),collapse="+")
      else
        s.formula = paste(paste("s(",smooth.name,",k=",num.knots,",bs='tp')",sep=""),collapse="+")
    }

    if(is.null(offset.name)) {
      o.formula = NULL
    } else {
      o.formula = paste("+offset(", offset.name, ")", sep = "")
    }

    if(interactions)
    {
      ll.int.formula = expand.grid(linear.name, linear.name) %>%
        apply(1, function(x) paste(x, collapse = ":")) %>%
        paste(collapse = "+")

      l.formula = paste(l.formula, ll.int.formula, collapse = "+")

      ## Number of knots might be a problem if diff for separate smooth.terms

      if(smooth.penalty==1)
      {
        s.formula = paste(paste("ti(",smooth.name,",k=",num.knots,",bs='ts')",
                                sep=""),collapse="+")
        ls.int.formula = expand.grid(smooth.name, linear.name) %>%
          apply(1, function(x) paste0("ti(",x[1],",by=",x[2],
                                      ",k=",num.knots,",ordered=T,bs='ts')")) %>%
          paste(collapse = "+")

        ss.int.formula = expand.grid(smooth.name, smooth.name) %>%
          apply(1, function(x) paste0("ti(",x[1],",",x[2],
                                      ",k=",num.knots,",bs='ts')")) %>%
          paste(collapse = "+")
      }
      else
      {
        s.formula = paste(paste("ti(",smooth.name,",k=",num.knots,",bs='tp')",
                                sep=""),collapse="+")
        ls.int.formula = expand.grid(linear.name, smooth.name) %>%
          apply(1, function(x) paste0("ti(",x[1],",by=",x[2],
                                      ",k=",num.knots,",ordered=T,bs='tp')")) %>%
          paste(collapse = "+")

        ss.int.formula = expand.grid(smooth.name, smooth.name) %>%
          apply(1, function(x) paste0("ti(",x[1],",",x[2],
                                      ",k=",num.knots,",bs='tp')")) %>%
          paste(collapse = "+")
      }

      s.formula = paste(s.formula, ls.int.formula, ss.int.formula, collapse = "+")
    }
  }
  else
  {
    fm = deparse(formula) %>% paste(collapse = "") %>%
      (function(x) {gsub("[[:space:]]", "", x)}) %>% strsplit("~") %>% unlist()
    # The anonymous function is only defined to circumvent the note that we get
    # for using . and which is not allowed in CRAN

    r.formula = fm[1]

    if(grepl("X+", fm[2])) {
      l.formula = "X"
    } else {
      l.formula = NULL
    }

    #fm1 = strsplit(fm[2], "X+", fixed = T) %>% unlist() %>% paste(collapse = "")
    fm1 = gsub("X+", "", fm[2], fixed = T)
    if(is.na(fm1) || nchar(fm1) == 0) {
      s.formula = NULL
    } else {
      s.formula = fm1
    }

    if(is.null(offset.name)) {
      o.formula = NULL
    } else {
      o.formula = paste("+offset(", offset.name, ")", sep = "")
    }
  }

  formulae = list(response = r.formula, linear = l.formula,
                 smooth = s.formula, offset = o.formula)
  return(formulae)
}

detect_one_two = function(response) {

  if(grepl("cbind", response)) {
    y = gsub("cbind(", "", response, fixed = T)
    y = gsub(")", "", y, fixed = T)
    y = strsplit(y, ",")
    y = unlist(y)
  } else {
    y = response
  }
  return(y)
}


lasso_response = function(response.name, family, data, train.weights)
{
  y = NULL
  if(family == "binomial" && length(response.name) == 2)
  {
    y = data %>% dplyr::select(response.name[2:1]) %>% as.matrix()
  }
  else if(family == "cox")
  {
    response = data %>% dplyr::pull(response.name)
    y = survival::Surv(response, train.weights)
  }
  else
  {
    y = data %>% dplyr::pull(response.name)
  }
  return(y)
}


#' Internal Function
#'
#' Undocumented function. Do not use directly
#'
#' @param x,y Vectors of the same length
meandist = function(x,y)
{
  d1 = sum((x-y)^2)
  d1 = sqrt(d1/length(x))
  return(as.numeric(d1))
}

#' Internal Function
#'
#' Undocumented function. Do not use directly
#'
#' @param x,y Vectors of the same length
nzeros = function(x,y=NULL)
{
  z1 = sum(x==0)
  z2 = NULL
  if(!is.null(y)) z2 = sum((x==0) & (y==0))
  return(c(z1,z2))
}

#' Internal Function
#'
#' Undocumented function. Do not use directly
#'
#' @param data The data with value for all the linear and smooth predictors
#' @param families List of two families as returned by \code{find_family}
#' @param formulae List of formulae as returned by \code{formula_setup}
#' @param weights Vector with values of the weights variable if it exists.
#'   \code{NULL} otherwise.
#' @inheritParams gamlassoChecks
#'
#' @importFrom dplyr one_of "%>%"
#' @importFrom stats as.formula model.matrix cor
lasso_gam_loop = function(data, response.name, families, formulae,
                          num.iter, tolerance, offset.name, weights, seed)
{
  off = train.offset = rep(0, nrow(data))
  if(!is.null(offset.name)) train.offset = data %>% dplyr::pull(offset.name)
  off = off + train.offset
  set.seed(seed)

  ## INITIAL GAM STEP
  eqn = as.formula(paste(formulae$response, "~", formulae$smooth))
  gam.model = mgcv::gam(eqn, family = families$gam_family, data = data,
                        weights = weights, offset = off)

  gam.pred = mgcv::predict.gam(gam.model, newdata = data) %>% as.numeric()
  gam.coefs = gam.model$coefficients %>% as.numeric()

  # gam prediction does not include offset if not specified in formula
  # We could just have gam.pred as the new offset but we would need to
  # change the formula in gam which leads to problems while predicting
  off = gam.pred + train.offset

  ## INITIAL LASSO STEP
  #xmat = model.matrix(as.formula(paste("~", formulae$linear)), data = data)[,-1]
  xmat = as.formula(paste("~", formulae$linear)) %>% model.matrix(data = data) %>%
    (function(x) {x[,-1]}) %>% as.matrix()
  # The anonymous function is only defined to circumvent the note that we get
  # for using . and which is not allowed in CRAN

  y = lasso_response(response.name, families$lasso_family, data, weights)

  ## For cox fits specifying intercept to be T or F or NULL all give warnings
  ## Although the output with intercept=F is the same as specifying nothing
  if(families$lasso_family == "cox")
  {
    cv.lasso.model = glmnet::cv.glmnet(x=xmat, y=y, offset=off,
                                       family=families$lasso_family, maxit=1000)
    lasso.model = glmnet::glmnet(x=xmat, y=y, offset=off,
                                 family=families$lasso_family, maxit=1000)
  }
  else
  {
    cv.lasso.model = glmnet::cv.glmnet(x=xmat, y=y, offset=off, intercept=F,
                                       family=families$lasso_family, maxit=1000)
    lasso.model = glmnet::glmnet(x=xmat, y=y, offset=off, intercept=F,
                                 family=families$lasso_family, maxit=1000)
  }


  lasso.pred = glmnet::predict.cv.glmnet(cv.lasso.model, newx=xmat, s="lambda.min",
                                         newoffset=0) %>% as.numeric()
  lasso.coefs = glmnet::coef.cv.glmnet(cv.lasso.model, s="lambda.min") %>% as.numeric()

  # glmnet prediction include offset
  # but suppressed by newoffset=0
  off = lasso.pred + train.offset

  convergence = FALSE

  convmat = matrix(0, nrow = 6, ncol = 1)
  convmat[5,1] = nzeros(lasso.coefs)

  for(i in 1:num.iter)
  {
    set.seed(seed)

    ## GAM STEP
    gam.model = mgcv::gam(eqn, family = families$gam_family, data = data,
                          weights = weights, offset = off)

    new.gam.pred = mgcv::predict.gam(gam.model, newdata = data) %>% as.numeric()
    new.gam.coefs = gam.model$coefficients %>% as.numeric()

    off = new.gam.pred + train.offset

    ## LASSO STEP
    ## For cox fits specifying intercept to be T or F or NULL all give warnings
    ## Although the output with intercept=F is the same as specifying nothing
    if(families$lasso_family == "cox")
    {
      cv.lasso.model = glmnet::cv.glmnet(x=xmat, y=y, offset=off,
                                         family=families$lasso_family, maxit=1000)
      lasso.model = glmnet::glmnet(x=xmat, y=y, offset=off,
                                   family=families$lasso_family, maxit=1000)
    }
    else
    {
      cv.lasso.model = glmnet::cv.glmnet(x=xmat, y=y, offset=off, intercept=F,
                                         family=families$lasso_family, maxit=1000)
      lasso.model = glmnet::glmnet(x=xmat, y=y, offset=off, intercept=F,
                                   family=families$lasso_family, maxit=1000)
    }

    new.lasso.pred = glmnet::predict.cv.glmnet(cv.lasso.model, newx=xmat,
                                               s="lambda.min", newoffset=0) %>% as.numeric()
    new.lasso.coefs = glmnet::coef.cv.glmnet(cv.lasso.model, s="lambda.min") %>% as.numeric()

    off = new.lasso.pred + train.offset

    ## CHECK CONVERGENCE
    if(meandist(lasso.pred, new.lasso.pred) < tolerance && meandist(gam.pred, new.gam.pred) < tolerance)
    {
      convergence = TRUE
      break
    }
    else
    {
      convnew = c(meandist(lasso.coefs, new.lasso.coefs),
                  cor(lasso.coefs, new.lasso.coefs),
                  meandist(gam.coefs, new.gam.coefs),
                  cor(gam.coefs, new.gam.coefs),
                  nzeros(lasso.coefs, new.lasso.coefs))
      convmat = cbind(convmat, convnew)

      lasso.pred = new.lasso.pred
      gam.pred = new.gam.pred
    }

  }

  return(list(lasso = list(cv.model = cv.lasso.model, model = lasso.model),
              gam = gam.model, convergence = convergence, convmat = convmat))
}



## Old gam and lasso loop where lasso was fitted first and then the gam
## Issue was that lasso fitted intercepts on top of intercept from gam
## We got rid of it by using intercept=F and reversing the order so we
## always start with some fitted intercept. If we didn't reverse we'd
## need another extra step just to fit a first intercept.
# lasso_gam_loop = function(data, response.name, families, formulae,
#                           num.iter, tolerance, offset.name, status, seed)
# {
#   off = train.offset = rep(0, nrow(data))
#   if(!is.null(offset.name)) train.offset = data %>% dplyr::pull(offset.name)
#   off = off + train.offset
#   set.seed(seed)
#
#   ## INITIAL LASSO STEP
#
#   xmat = model.matrix(as.formula(paste("~", formulae$linear)), data = data)[,-1]
#   response = data %>% dplyr::pull(response.name)
#   if(families$lasso_family != "cox")
#     y = response
#   else
#     y = survival::Surv(response, status)
#
#   cv.lasso.model = glmnet::cv.glmnet(x=xmat, y=y, offset=off, family=families$lasso_family, maxit=1000)
#   lasso.model = glmnet::glmnet(x=xmat, y=y, offset=off, family=families$lasso_family, maxit=1000)
#
#   lasso.pred = as.numeric(glmnet::predict.cv.glmnet(cv.lasso.model, newx=xmat, s="lambda.min", newoffset=0))
#   lasso.coefs = as.numeric(glmnet::coef.cv.glmnet(cv.lasso.model, s="lambda.min"))
#
#   # glmnet prediction include offset
#   off = lasso.pred + train.offset
#
#   ## INITIAL GAM STEP
#
#   eqn = as.formula(paste(response.name, "~", formulae$smooth))
#   gam.model = mgcv::gam(eqn, family = families$gam_family, data = data,
#                         weights = status, offset = off)
#
#   gam.pred = as.numeric(mgcv::predict.gam(gam.model, newdata = data))
#   gam.coefs = as.numeric(gam.model$coefficients)
#
#   # gam prediction does not include offset if not specified in formula
#   # We could just have gam.pred as the new offset but we would need to
#   # change the formula in gam which leads to problems while predicting
#   off = gam.pred + train.offset
#
#   convergence = FALSE
#
#   convmat = matrix(0, nrow = 6, ncol = 1)
#   convmat[5,1] = nzeros(lasso.coefs)
#
#   # lasso_dist = lasso_corr = vector()
#   # gam_dist = gam_corr = vector()
#   # lasso_zeros = matrix(nrow = 2, ncol = 0)
#   # lasso_zeros[1,1] = nzeros(lasso.coefs)
#
#   for(i in 1:num.iter)
#   {
#     set.seed(seed)
#     ## LASSO STEP
#
#     cv.lasso.model = glmnet::cv.glmnet(x=xmat, y=y, offset=off, family=families$lasso_family, maxit=1000)
#     lasso.model = glmnet::glmnet(x=xmat, y=y, offset=off, family=families$lasso_family, maxit=1000)
#
#     new.lasso.pred = as.numeric(glmnet::predict.cv.glmnet(cv.lasso.model, newx=xmat, s="lambda.min", newoffset=0))
#     new.lasso.coefs = as.numeric(glmnet::coef.cv.glmnet(cv.lasso.model, s="lambda.min"))
#
#     off = new.lasso.pred + train.offset
#
#     ## GAM STEP
#
#     gam.model = mgcv::gam(eqn, family = families$gam_family, data = data,
#                           weights = status, offset = off)
#
#     new.gam.pred = as.numeric(mgcv::predict.gam(gam.model, newdata = data))
#     new.gam.coefs = as.numeric(gam.model$coefficients)
#
#     off = new.gam.pred + train.offset
#
#     if(meandist(lasso.pred, new.lasso.pred) < tolerance && meandist(gam.pred, new.gam.pred) < tolerance)
#     {
#       convergence = TRUE
#       break
#     }
#     else
#     {
#       convnew = c(meandist(lasso.coefs, new.lasso.coefs),
#                   cor(lasso.coefs, new.lasso.coefs),
#                   meandist(gam.coefs, new.gam.coefs),
#                   cor(gam.coefs, new.gam.coefs),
#                   nzeros(lasso.coefs, new.lasso.coefs))
#       convmat = cbind(convmat, convnew)
#
#       lasso.pred = new.lasso.pred
#       gam.pred = new.gam.pred
#     }
#
#   }
#
#   return(list(lasso = list(cv.model = cv.lasso.model, model = lasso.model),
#               gam = gam.model, convergence = convergence, convmat = convmat))
# }


