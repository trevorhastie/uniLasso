#' Fit a cross-validated univariate guided lasso model, followed by a lasso polish.
#'
#' This function has two stages. In the first we fit a
#' univariate-guided sparse regression `uniLasso` model using
#' cross-validation to select the lasso penalty
#' parameter (using \code{s = "lambda.min"}). In the second stage, we use the predictions from this chosen
#' model as an offset, and fit a cross-validated unrestricted lasso
#' model. For squared error loss, this means we post-fit a lasso model
#' to the residuals. Conveniently, it returns an object that inherits
#' from \code{cv.glmnet}, in which the two models are _stitched_
#' together. What this means is that the chosen coefficients from the
#' first model are added to the coefficients from the second, and other related components are updated as well.
#' This means at predict time we do not have to fiddle with offsets.  All
#' of the methods for \code{cv.glmnet} can be applied, such as
#' `predict`, `plot`, `coef`, `print`, and `assess.glmnet`.
#'
#'
#' @param x Input matrix, of dimension \code{nobs x nvars}; each row is
#' an observation vector.
##' @param y Response variable. Quantitative for \code{family = "gaussian"} or
#' \code{family = "poisson"} (non-negative counts). For \code{family="binomial"},
#' should be a numeric vector consisting of 0s and 1s. For \code{family="cox"},
#' y should be a two-column matrix with columns named 'time' and 'status'.
#' The latter is a binary variable, with '1' indicating death, and '0'
#' indicating right-censored.
#'
#' @param family one of "gaussian","binomial" or "cox". Currently only these families are implemented.
#' In the future others will be added.
#' @param weights optional vector of non-negative weights, default is NULL which results in all weights = 1.
#' @param \ldots additional arguments passed to \code{cv.uniLasso} and
#'     \code{cv.glmnet}. Note: by defaults \code{cv.uniLasso()} uses
#'     `standardize=FALSE`, and `cv.glmnet()` uses
#'     `standardize=TRUE`. These are both the sensible defaults for
#'     this function. Users can supply `standardize=FALSE` (via the
#'     \ldots argument) which will overide the `cv.glmnet()`
#'     default. Users should avoid using `standardize=TRUE`, since
#'     this will affect the first stage model as well, where this is
#'     not a suitable choice.
#'
#' @return An object of class \code{"polish.unilasso"} that inherits
#'     from class \code{"cv.glmnet"}. The \code{"glmnet.fit"} is the
#'     _stitched_ second-stage model, from which predictions are
#'     made. An additional component named \code{"cv.uniLasso"} is the
#'     first stage model.
#' @examples
#'
#' # Gaussian data, p=1000, n=300, SNR=1  "medium SNR"
#' # use the built-in simulate function to create Gaussian data
#' set.seed(101)
#' data <- simulate_uniLasso("medium-SNR")
#' attach(data) # has components "x","y","xtest","ytest","mutest","sigma"
#' pfit <- polish.uniLasso(x,y)
#' plot(pfit)
#' pred <- predict(pfit, newx = xtest,  s = "lambda.min") # ie predict from a "cv.glmnet" object.
#' mean((ytest-pred)^2) # test error
#' print(pfit)
#' print(pfit$glmnet.fit)
#' plot(pfit$glmnet.fit) # coefficient plot of the second stage
#' plot(pfit$cv.uniLasso) # cv.glmnet plot of the first stage
#' plot(pfit$cv.uniLasso$glmnet.fit) # coefficient plot of the first stage
#'
#'
#' # Binomial response
#'
#' yb =as.numeric(y>0)
#' pfitb = polish.uniLasso(x, yb, family="binomial")
#' predict(pfitb, xtest[1:10,], type="response") # predict at default s = "lambda.1se"
#' plot(pfitb)
#' plot(pfitb$glmnet.fit) # plot second stage lasso coefficient path
#' plot(pfitb$cv.uniLasso) # plot first stage cv.uniLasso results

#' # Cox response
#'
#' set.seed(10101)
#' N = 1000
#' p = 30
#' nzc = p/3
#' x = matrix(rnorm(N * p), N, p)
#' beta = rnorm(nzc)
#' fx = x[, seq(nzc)] %*% beta/3
#' hx = exp(fx)
#' ty = rexp(N, hx)
#' tcens = rbinom(n = N, prob = 0.3, size = 1)  # censoring indicator
#' y = cbind(time = ty, status = 1 - tcens)  # y=Surv(ty,1-tcens) with library(survival)
#' pfitc = polish.uniLasso(x, y, family = "cox")
#' plot(pfitc)
#' plot(pfitc$cv.uniLasso)
#'
#' @export

polish.uniLasso <- function(x,y,family=c("gaussian","binomial","cox"),weights=NULL,...){
    this.call = match.call()
    family=match.arg(family)
    fit0 = cv.uniLasso(x,y,family=family,weights=weights,...)
    pred0 = predict(fit0,s="lambda.min",newx=x)
    fit1 = cv.glmnet(x,y,offset=pred0,weights=weights,family=family,...)
    ## Now we stitch them together
    glmnet.fit1 = fit1$glmnet.fit
    glmnet.fit0= fit0$glmnet.fit
    index=fit0$index[1,] # indexes lambda.min
    beta1=glmnet.fit1$beta
    nc=ncol(beta1)
    beta0=glmnet.fit0$beta[,index,drop=FALSE]
    glmnet.fit1$beta <- beta1 + beta0[,rep(1,nc)]
    if(!is.null(glmnet.fit1$a0)){# Interecpt is present
        glmnet.fit1$a0=glmnet.fit1$a0+glmnet.fit0$a0[index]
    }
    glmnet.fit1$df = apply(glmnet.fit1$beta!=0,2,sum)
    glmnet.fit1$offset=FALSE
    fit1$glmnet.fit=glmnet.fit1
    ## end of stitch
    fit1$cv.uniLasso = fit0
    fit1$call=this.call
    class(fit1)=c("polish.uniLasso",class(fit1))
    fit1
    }

