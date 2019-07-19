## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(message=F, warning = F, cache=T, echo=T, eval=T,
                      fig.pos="H", fig.align="center", comment='#>')

## ---- echo=FALSE---------------------------------------------------------
tab <- matrix(c("\u2713","",rep("\u2713",4)),2,3)
rownames(tab) = c("$\\beta$","$f_j$'s")
colnames(tab) = c("None","$\\ell_1$","$\\ell_2$")
knitr::kable(tab, align = "rccc")

## ------------------------------------------------------------------------
library(tidyverse)
library(plsmselect)

data(simData)

## ---- echo=FALSE---------------------------------------------------------
knitr::kable(head(simData), format = "html", digits=2,
             caption="Table 1. First 6 samples of the simulated data set: simData.") %>% 
  kableExtra::kable_styling(bootstrap_options = "striped", font_size = 11.5)

## ----normalfit-----------------------------------------------------------
## Create model matrix X corresponding to linear terms
## (necessary for the formula option of gamlasso below)
simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]

## The formula approach
gfit = gamlasso(Yg ~ X +
                  s(z1, k=5, bs="ts") +
                  s(z2, k=5, bs="ts") +
                  s(z3, k=5, bs="ts") +
                  s(z4, k=5, bs="ts"),
                data = simData,
                seed = 1)

## ----normalfitalt, eval = F----------------------------------------------
#  ## The term specification approach
#  gfit = gamlasso(response = "Yg",
#                  linear.terms = paste0("x",1:10),
#                  smooth.terms = paste0("z",1:4),
#                  data = simData,
#                  linear.penalty = "l1",
#                  smooth.penalty = "l1",
#                  num.knots = 5,
#                  seed = 1)

## ------------------------------------------------------------------------
# mgcv::gam object:
class(gfit$gam)
# glmnet::cv.glmnet object
class(gfit$cv.glmnet)

## ------------------------------------------------------------------------
summary(gfit)

## ----plotgam, fig.width=6, fig.height=6----------------------------------
## Plot the estimates of the smooth effects:
plot(gfit$gam, pages=1)

## ----fittedvsobserved, fig.width=6, fig.height=6-------------------------
## Plot fitted versus observed values:
plot(simData$Yg, predict(gfit), xlab = "Observed values", ylab = "Fitted Values")

## ----poifit--------------------------------------------------------------
## Create model matrix X corresponding to linear terms
## (necessary for the formula option of gamlasso below)
simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]

## Poisson response. Formula approach.
pfit = gamlasso(Yp ~ X + 
                  s(z1, bs="ts", k=5) + 
                  s(z2, bs="ts", k=5) + 
                  s(z3, bs="ts", k=5) + 
                  s(z4, bs="ts", k=5),
                data = simData,
                family = "poisson",
                seed = 1)

## ----poifitalt, eval=FALSE-----------------------------------------------
#  ## Poisson response. Term-specification approach.
#  pfit = gamlasso(response = "Yp",
#                  linear.terms = paste0("x",1:10),
#                  smooth.terms = paste0("z",1:4),
#                  data = simData,
#                  linear.penalty = "l1",
#                  smooth.penalty = "l1",
#                  family = "poisson",
#                  num.knots = 5,
#                  seed = 1)

## ------------------------------------------------------------------------
coef(pfit$cv.glmnet, s="lambda.min")

## ---- fig.width=6, fig.height=3------------------------------------------
par(mfrow=c(1,2))
plot(pfit$gam, select=1) # estimate of smooth term z1
plot(pfit$gam, select=2) # estimate of smooth term z2

## ---- fig.width=6, fig.height=6------------------------------------------
plot(predict(pfit, type="response"), exp(simData$lp), xlab="predicted count", ylab="true expected count")

## ----binfit, eval=FALSE--------------------------------------------------
#  ## Create model matrix X corresponding to linear terms
#  ## (necessary for the formula option of gamlasso below)
#  simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]
#  
#  ## Bernoulli trials response
#  bfit = gamlasso(Yb ~ X +
#                    s(z1, bs="ts", k=5) +
#                    s(z2, bs="ts", k=5) +
#                    s(z3, bs="ts", k=5) +
#                    s(z4, bs="ts", k=5),
#                  data = simData,
#                  family = "binomial",
#                  seed = 1)

## ----binfitalt, eval = F-------------------------------------------------
#  ## The term specification approach
#  bfit = gamlasso(response = "Yb",
#                  linear.terms = paste0("x",1:10),
#                  smooth.terms = paste0("z",1:4),
#                  data = simData,
#                  family="binomial",
#                  linear.penalty = "l1",
#                  smooth.penalty = "l1",
#                  num.knots = 5,
#                  seed = 1)

## ---- eval=FALSE---------------------------------------------------------
#  summary(bfit)
#  plot(bfit$gam, pages=1)
#  pred.prob <- predict(bfit, type="response")
#  true.prob <- exp(simData$lp)/(1+exp(simData$lp))
#  plot(pred.prob, true.prob, xlab="predicted probability", ylab="true probability")

## ----binomfit, eval=FALSE------------------------------------------------
#  ## Create model matrix X corresponding to linear terms
#  ## (necessary for the formula option of gamlasso below)
#  simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]
#  
#  ## Binomial counts response. Formula approach.
#  bfit2 = gamlasso(cbind(success,failure) ~ X +
#                     s(z1, bs="ts", k=5) +
#                     s(z2, bs="ts", k=5) +
#                     s(z3, bs="ts", k=5) +
#                     s(z4, bs="ts", k=5),
#                   data = simData,
#                   family = "binomial",
#                   seed = 1)

## ----binomfitalt, eval = F-----------------------------------------------
#  ## Binomial counts response. Term specification approach
#  bfit2 = gamlasso(c("success","failure"),
#                   linear.terms=paste0("x",1:10),
#                   smooth.terms=paste0("z",1:4),
#                   data=simData,
#                   family = "binomial",
#                   linear.penalty = "l1",
#                   smooth.penalty = "l1",
#                   num.knots = 5,
#                   seed=1)

## ---- eval=FALSE---------------------------------------------------------
#  summary(bfit2)
#  plot(bfit2$gam, pages=1)
#  pred.prob <- predict(bfit2, type="response")
#  true.prob <- exp(simData$lp)/(1+exp(simData$lp))
#  plot(pred.prob, true.prob, xlab="predicted probability", ylab="true probability")

## ----coxfit--------------------------------------------------------------
## Create model matrix X corresponding to linear terms
## (necessary for the formula option of gamlasso below)
simData$X = model.matrix(~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10, data=simData)[,-1]

# Censored time-to-event response. Formula approach.
cfit = gamlasso(time ~ X +
                  s(z1, bs="ts", k=5) +
                  s(z2, bs="ts", k=5) +
                  s(z3, bs="ts", k=5) +
                  s(z4, bs="ts", k=5),
                data = simData,
                family = "cox",
                weights = "status",
                seed = 1)

## ----coxfitalt, eval=FALSE-----------------------------------------------
#  # Censored time-to-event response. Term specification approach.
#  cfit = gamlasso(response = "time",
#                  linear.terms = paste0("x",1:10),
#                  smooth.terms = paste0("z",1:4),
#                  data = simData,
#                  linear.penalty = "l1",
#                  smooth.penalty = "l1",
#                  family = "cox",
#                  weights="status",
#                  num.knots = 5,
#                  seed = 1)

## ----coxdiagnos, fig.width=6, fig.height=6-------------------------------
## Obtain and plot predicted cumulative baseline hazard:
H0.pred <- cumbasehaz(cfit)

time.seq <- seq(0, 60, by=1)
plot(time.seq, H0.pred(time.seq), type="l", xlab = "Time", ylab="",
     main = "Predicted Cumulative \nBaseline Hazard")

## ----coxpredict, fig.width=6, fig.height=6-------------------------------
## Obtain predicted survival at days 1,2,3,...,60:
S.pred <- predict(cfit, type="response", new.event.times=1:60)

## Plot the survival curve for sample (subject) 17:
plot(1:60, S.pred[17,], xlab="time (in days)", ylab="Survival probability", type="l")

## ----simdata, eval=FALSE-------------------------------------------------
#  generate.inputdata <- function(N) {
#  
#    id <- paste0("i",1:N)
#    ## Define linear predictors
#    linear.pred <- matrix(c(rbinom(N*3,size=1,prob=0.2),
#                            rbinom(N*5,size=1,prob=0.5),
#                            rbinom(N*2,size=1,prob=0.9)),
#                          nrow=N)
#    colnames(linear.pred) = paste0("x", 1:ncol(linear.pred))
#  
#    ## Define smooth predictors
#    smooth.pred = as.data.frame(matrix(runif(N*4),nrow=N))
#    colnames(smooth.pred) = paste0("z", 1:ncol(smooth.pred))
#  
#    ## Coalesce the predictors
#    dat = cbind.data.frame(id, linear.pred, smooth.pred)
#  
#    return(dat)
#  
#  }
#  
#  ## Simulate N input data observations:
#  N <- 100
#  simData2 <- generate.inputdata(N)

## ---- eval=FALSE---------------------------------------------------------
#  ## "True" coefficients of linear terms (last 7 are zero):
#  beta <- c(1, -0.6, 0.5, rep(0,7))
#  
#  ## "True" nonlinear smooth terms:
#  f1 <- function(x) x^3
#  f2 <- function(x) sin(x*pi)
#  
#  ## Calculate "True" linear predictor (lp) based on above data (simData2)
#  simData2$lp <- as.numeric(as.matrix(simData2[,paste0("x",1:10)]) %*% beta + f1(simData2$z1) + f2(simData2$z2))

## ---- eval=FALSE---------------------------------------------------------
#  ## Simulate gaussian response with mean lp:
#  simData2$Yg = simData2$lp + rnorm(N, sd = 0.1)
#  
#  ## Simulate bernoulli trials with success probability exp(lp)/(1+exp(lp))
#  simData2$Yb = map_int(exp(simData2$lp), ~ rbinom(1, 1, p = ( ./(1+.) ) ) )
#  
#  ## Simulate Poisson response with log(mean) = lp
#  simData2$Yp = map_int(exp(simData2$lp), ~rpois(1, .) )

## ---- eval=FALSE---------------------------------------------------------
#  ## Simulate binomial counts with success probability exp(lp)/(1+exp(lp))
#  sizes = sample(10, nrow(simData2), replace = TRUE)
#  # Simulate success:
#  simData2$success = map2_int(exp(simData2$lp), sizes, ~rbinom(1, .y, p = ( .x/(1+.x) )))
#  # Calculate failure:
#  simData2$failure = sizes - simData2$success

## ---- eval=FALSE---------------------------------------------------------
#  ## Function that simulates N samples of censored event times (time, status)
#  ## Event times ~ Weibull(lambda0, rho) with linear predictor lp
#  ## rho=1 corresponds to exponential distribution
#  simulWeib <- function(N, lambda0, rho, lp)
#  {
#  
#    # Censoring times ~ Exponential(lambdaC)
#    lambdaC=0.0005 # very mild censoring
#  
#    # Simulate Weibull latent event times
#    v <- runif(n=N)
#    Tlat <- (- log(v) / (lambda0 * exp(lp)))^(1 / rho)
#  
#    # Simulate censoring times
#    C <- rexp(n=N, rate=lambdaC)
#  
#    # Calculate follow-up times and event indicators
#    time <- pmin(Tlat, C)
#    status <- as.numeric(Tlat <= C)
#  
#    data.frame(time=time, status=status)
#  
#  }
#  
#  ## Set parameters of Weibull baseline hazard h0(t) = lambda*rho*t^(rho-1):
#  lambda0 <- 0.01; rho <- 1;
#  
#  ## Simulate (times, status) from above censored Weibull model and cbind to our data:
#  simData2 <- cbind.data.frame(simData2, simulWeib(N, lambda0, rho, simData2$lp))

