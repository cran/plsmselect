## CHECKING FUNCTION

#' Internal Function
readconfirm = function()
{
  conf = readline(prompt = "Please enter y or n: ")
  if(!(conf %in% c("y","n")))
    return(readconfirm())
  return(as.character(conf))
}

#' Checking data before fitting gamlasso
#'
#' This function checks if the arguments entered for fitting a gamlasso model
#' are compatible with each other. Not recommended to call directly. Only use
#' if cleaning data prior to fitting \code{\link{gamlassoFit}}
#'
#' @param data The training data for fitting the model
#' @param response.name The name of the response variable. Vector of two if
#'   \code{family = "binomial"}
#' @param linear.name The names of the variables to be used as linear predictors
#' @param smooth.name The names of the variables to be used as smoothers
#' @param linear.penalty The penalty used on the linear predictors. Can be 0, 1 or 2
#' @param smooth.penalty The penalty used on the smoothers. Can be 1 or 2
#' @param offset.name The name of the offset variable. \code{NULL} (default) if not provided
#' @param weights.name The name of the weights variable. \code{NULL} (default)
#'   if not provided. See \code{Details} of \code{\link{gamlasso}}.
#' @inheritParams gamlasso
#'
#' @return \code{gamlassoChecks} produces a series of logical values:
#'   \code{allcheck} indicating if the arguments passed all the checks,
#'   \code{fit.smoothgam} indicating if there aren't any linear predictors and
#'   a model with only smoothers should be fitted, \code{fit.glmnet}
#'   is the counterpart for smooth predictors. It also returns the cleaned
#'   (if needed) arguments as a list named \code{cleandata} who's elements are:
#'
#'   \tabular{ll}{
#'     \code{train.data} \tab The training data with unnecessary columns deleted \cr
#'     \code{linear.name}, \code{smooth.name}, \code{num.knots} \tab The changed
#'      variable names and number of knots \cr
#'     \code{linear.penalty}, \code{smooth.penalty} \tab The changed penalties for linear and smooth
#'     terms. Reset to their default values only
#'     in the rare case of too few predictors
#'   }
#'
#' @export
#' @importFrom dplyr one_of "%>%"
#'
#' @note The arguments \code{offset.name}, \code{num.iter}, \code{tolerance}
#'   and \code{seed} are not currently not being used in testing.
#'
#' @examples ## Usage similar to gamlassoFit
gamlassoChecks = function(data, response.name, linear.name, smooth.name,
                           family, linear.penalty, smooth.penalty,
                           offset.name, weights.name, num.knots,
                           num.iter, tolerance, seed, prompts)
{
  if(prompts) cat("\nPreliminary checking ... \n")

  allcheck = TRUE
  p = rep("n", 18)
  fit.glmnet = fit.smoothgam = FALSE
  # changed in 13 and 9 resp

  ## CHECK IF data IS A DATA FRAME. IF NOT COERCE IT AND GIVE A WARNING
  if(!("data.frame" %in% class(data)))
  {
    msg = paste("Coercing data to a data.frame could lead to loss of",
                "variable names. Reenter data as data.frame")
    if(prompts)
    {
      cat(paste("\n\n", msg, "?"))
      p[1] = readconfirm()
    }
    else
    {
      warning(paste(msg, "recommended."))
    }
    if(p[1] == "n") data = as.data.frame(data)
  }


  ## ENSURE THAT family IS ENTERED AS A STRING
  if(class(family) != "character")
  {
    msg = paste("Argument 'family' should be of class 'character'.")
    if(prompts)
    {
      cat(paste("\n\n", msg))
      p[2] = "y"
    }
    else
    {
      stop(msg)
    }
  }

  ## CHECK THAT THE RESPONSE VARIABLE IS IN THE DATA
  if(!all(response.name %in% colnames(data)))
  {
    msg = paste(response.name, "is not a variable name in data.")
    if(prompts)
    {
      cat(paste("\n\n", msg))
      p[3] = "y"
    }
    else
    {
      stop(msg)
    }
  }

  ## CHECK THAT THE FEATURE NAMES ALSO EXIST IN THE DATA
  feature.names = c(linear.name, smooth.name)
  badnames = base::setdiff(feature.names, colnames(data))
  if(length(badnames) > 0)
  {
    msg = paste("The feature(s)", paste(badnames, collapse = ", "),
                "aren't in the dataset and will be ignored.")
    if(prompts)
    {
      cat(paste("\n\n", msg, "Do you believe these variables are",
                "important and wish to enter a new dataset?"))
      p[4] = readconfirm()
    }
    else
    {
      warning(msg)
    }
    if(p[4] == "n")
    {
      linear.name = base::setdiff(linear.name, badnames)
      smooth.name = base::setdiff(smooth.name, badnames)
    }
  }

  ## CHECK THAT THE OFFSET VARIABLE IS IN THE DATA
  if(!is.null(offset.name))
  {
    if(!(offset.name %in% colnames(data)))
    {
      msg = paste(offset.name, "is not a variable name in data.")
      if(prompts)
      {
        cat(paste("\n\n", msg))
        p[5] = "y"
      }
      else
      {
        stop(msg)
      }
    }
  }

  # CHECK THAT smooth.name AND linear.name HAVE NOTHING IN COMMON
  overlap = intersect(linear.name, smooth.name)
  if(length(overlap) > 0)
  {
    msg = paste("The variable(s)", paste(overlap, collapse = ", "),
                "overlap in the linear and smooth terms and will be ignored.")
    if(prompts)
    {
      cat(paste("\n\n", msg, "Do you believe these variables are important",
                "and wish to reenter them in their proper category?"))
      p[6] = readconfirm()
    }
    else
    {
      warning(msg)
    }
    if(p[6] == "n")
    {
      linear.name = base::setdiff(linear.name, overlap)
      smooth.name = base::setdiff(smooth.name, overlap)
    }
  }

  ## CHECK IF THE RELEVANT COLUMNS OF DATA HAVE CORRECT CLASS
  fit.names = c(feature.names, response.name, offset.name, weights.name)
  dataclass = data %>% dplyr::select(fit.names) %>% sapply(class)
  badcols = fit.names[!(dataclass %in% c("factor", "numeric", "integer", "matrix"))]
  if(length(badcols) > 0)
  {
    msg = paste("The variables", paste(badcols, collapse = ", "),
                "were not of class factor or numeric and will be ignored.")
    if(prompts)
    {
      cat(paste("\n\n", msg, "Do you believe these variables are important",
                "and wish to reenter them with their correct class?"))
      p[7] = readconfirm()
    }
    else
    {
      warning(msg)
    }
    if(p[7] == "n") data = data %>% dplyr::select(-one_of(badcols))
  }

  ## THE DATA MAY HAVE EXTRA COLUMNS. CHECK IF THEY ARE IMPORTANT
  extra.names = base::setdiff(colnames(data), fit.names)
  if(length(extra.names) > 0)
  {
    msg = paste("There are", length(extra.names), "more variable(s) in the",
                "data than provided in the response, linear, smooth,",
                "offset or weights names. These will be ignored.")
    if(prompts)
    {
      cat(paste("\n\n", msg, "Do you believe these variables are important",
                "and wish to restart the fit with new variable names?"))
      p[8] = readconfirm()
    }
    # else
    # {
    #   warning(msg)
    # }
    if(p[8] == "n") data = data %>% dplyr::select(-one_of(extra.names))
  }

  ## CHECK IF A VARIABLE IS USELESS, i.e., HAS ONLY 1 VALUE
  onevalue = data %>% lapply(unique) %>% sapply(length)
  onevalue = colnames(data)[onevalue==1]
  if(length(onevalue) > 0)
  {
    msg = paste("The variable(s)", paste(onevalue, collapse = ", "),
                "only take one value and will be ignored.")
    if(prompts)
    {
      cat(paste("\n\n", msg))
    }
    else
    {
      warning(msg)
    }
    data = data %>% dplyr::select(-one_of(onevalue))
  }

  ## IF THERE ARE NO LINEAR VARIABLES THEN DON'T NEED PARAMETRIC PART OF GAMLASSO
  if(length(linear.name)==0)
  {
    msg = paste("No variables were provided as linear predictors so",
                "fitting a generalised additive model")
    if(prompts)
    {
      cat(paste("\n\n", msg, "Do you wish to quit and fit",
                "using mgcv::gam for more flexibility?"))
      p[9] = readconfirm()
    }
    else
    {
      warning(msg)
    }
    if(p[9] == "n") fit.smoothgam = TRUE
  }

  ## CHECK THAT SMOOTH VARIABLES ARE ALL numeric
  smoothclass = data %>% dplyr::select(smooth.name) %>% sapply(class)
  badsmooth = which(smoothclass != "numeric")
  if(length(badsmooth) > 0)
  {
    msg = paste("The smoothing variable(s)", paste(badsmooth, collapse = ", "),
                "are not numeric and will be ignored.")
    if(prompts)
    {
      cat("\n\n", msg, "Do you believe these variables are",
          "important and wish to reenter them with correct class?")
      p[10] = readconfirm()
    }
    else
    {
      warning(msg)
    }
    if(p[10] == "n") smooth.name = base::setdiff(smooth.name, badsmooth)
  }

  ## CHECKING length(num.knots) is 1 or length(smooth.name)
  if(!(length(num.knots) %in% c(1,length(smooth.name)) ))
  {
    msg = paste("'num.knots' must be of length 1 or ",length(smooth.name),
                " not ",length(num.knots),".", sep = "")
    if(prompts)
    {
      cat(paste("\n\n", msg))
      p[11] = "y"
    }
    else
    {
      stop(msg)
    }
  }

  ## CHECK THAT SMOOTH VARIABLES HAVE ENOUGH UNIQUE VALUES
  if(p[11] == "n")
  {
    b = data %>% dplyr::select(smooth.name) %>%
      lapply(unique) %>% sapply(length)
    badsmooth2 = smooth.name[!(b > num.knots)]
    if(length(badsmooth2) > 0)
    {
      msg = paste("The variable(s)", paste(badsmooth2, collapse = ", "),
                  "do not have enough unique values to be used as",
                  "smoothers and will be ignored.")
      if(prompts)
      {
        cat(paste("\n\n", msg, "Do you believe these variables are",
                  "important and wish to use them as a linear predictor?"))
        p[12] = readconfirm()
      }
      else
      {
        warning(msg)
      }
      if(p[12] == "n")
      {
        w = which(smooth.name %in% badsmooth2)
        smooth.name = smooth.name[-w]
        num.knots = num.knots[-w]
      }
      if(p[12] == "y") linear.name = base::union(linear.name, badsmooth2)
    }
  }

  ## CHECK IF THERE ARE NO SMOOTH VARIABLES
  if(length(smooth.name)==0)
  {
    msg = paste("No smoothing variables found so fitting a",
                ifelse(linear.penalty==0, "", "penalised"),
                "generalised linear model.")
    if(prompts)
    {
      cat(paste("\n\n", msg, "Do you wish to quit and fit using",
                ifelse(linear.penalty==0, "stats::glm", "glmnet::glmnet"),
                "for more flexibility?"))
      p[13] = readconfirm()
    }
    else
    {
      warning(msg)
    }
    if(p[13] == "n") fit.glmnet = TRUE
  }

  ## IF THERE ARE TOO FEW VARIABLES SUGGEST NOT PUTTING IN PENALTIES
  if(length(linear.name) + length(smooth.name) <= 10)
  {
    if(linear.penalty != 0 || smooth.penalty !=2)
    {
      msg = paste("Variable selection with too few predictors",
                  "could be numerically unstable.")
      if(prompts)
      {
        cat(paste("\n\n", msg, "Fit the model without penalties?"))
        p[14] = readconfirm()
      }
      else
      {
        warning(msg)
      }
      if(p[14] == "y")
      {
        linear.penalty = 0
        smooth.penalty = 2
      }
    }
  }

  ## FOR FITTING poisson OR binomial CHECK FOR POSITIVE INTEGER RESPONSES
  ## while fitting general binomial, input could be proportions. Keep this
  ## check for poisson. For binomial we need a way more complicated
  ## check and clean
  if(family %in% c("poisson", "binomial"))
  {
    response = data %>% dplyr::select(response.name)
    if( !( all(response >= 0) && all(response == floor(response)) ) )
    {
      msg = paste("Only positive integer responses are allowed for",
                  "fitting a", family, "model.")
      if(prompts)
      {
        cat(paste("\n\n", msg))
        p[15] = "y"
      }
      else
      {
        stop(msg)
      }
    }
  }
  ## Right now we are restricting binomial responses to only integers
  ## later on we could expand to proportions. So right now we don't need
  ## weights but later we probably will.
  if(family == "binomial")
  {
    if(! (length(response.name) %in% 1:2) )
    {
      msg = paste("Please enter one or two response variables.")
      if(prompts)
      {
        cat(paste("\n\n", msg))
        p[16] = "y"
      }
      else
      {
        stop(msg)
      }
    }
  }



  ## FOR FITTING cox STATUS = 1-CENSOR MUST BE SUPPLIED
  ## AND SHOULD MATCH ONE OF THE VARIABLES IN THE DATA
  if(family == "cox")
  {
    if(is.null(weights.name))
    {
      msg = paste("Status variable not found for fitting a Cox model.")
      if(prompts)
      {
        cat(paste("\n\n", msg))
        p[17] = "y"
      }
      else
      {
        stop(msg)
      }
    }
    else if(!(weights.name %in% colnames(data)))
    {
      msg = paste(weights.name, "is not a variable name in data.")
      if(prompts)
      {
        cat(paste("\n\n", msg))
        p[18] = "y"
      }
      else
      {
        stop(msg)
      }
    }
  }

  if("y" %in% p[-c(12,14)]) allcheck = FALSE

  cleandata = list(train.data = data, num.knots = num.knots,
                   linear.name = linear.name, linear.penalty = linear.penalty,
                   smooth.name = smooth.name, smooth.penalty = smooth.penalty)

  return(list(fit.glmnet = fit.glmnet, fit.smoothgam = fit.smoothgam,
              promptans = p, check = allcheck, clean = cleandata) )

}
