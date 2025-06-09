#' Fit a cross-validated univariate guided lasso model.
#'
#' Fit a univariate-guided sparse regression `uniLasso` model using cross-validation to select the second stage lasso penalty parameter. Conveniently,  it returns an object that inherits from \code{cv.glmnet}, so that
#' all of the methods for \code{cv.glmnet} can be applied, such as `predict`, `plot`, `coef`, `print`,
#'  and `assess.glmnet`.
#'
#' Fits a two stage lasso model and selects the penalty parameter by cross validation. First stage replaces each feature by the univariate fit for that feature.
#'   Second stage fits a (positive) lasso using the first stage features. Hence the second stage selects and
#'   modifies the coefficients of the first stage model, similar to the adaptive lasso. Leads to potentially sparser models.
#'
#'  For "binomial" family `y` is a binary response.
#'  For "cox" family, `y` should be a `Surv` object for right censored data,
#'  or a matrix with columns labeled 'time' and 'status'
#'  Although `glmnet` has more flexible options say for binary responses, and for `cox`
#'  responses, these are not yet implemented (but are possible and will appear in future versions).
#'  Likewise, other \code{glm} families are possible as well, but not yet implemented.
#'
#'  This is a one-visit function that returns a \code{'cv.glmnet'} object.
#'  You can make predictions from the whole path, at 'lambda.min' etc just like you can for
#'  a `'cv.glmnet object'`.
#'
#'  `loo = TRUE` means it uses the prevalidated loo fits (approximate loo or 'alo' for binomial and cox) for each univariate model as features to avoid overfitting in the second stage. The coefficients are then multiplied into the original univariate coefficients to get the final model.
#'
#'  `loo = FALSE` means it uses the univariate fitted predictor,  and hence it is a form of adaptive lasso, but tends to overfit.
#'  `lower.limits = 0` means `uniLasso` constrains the sign of the coefs in the second round to be that of the univariate fits.
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
#' @param family one of "gaussian","binomial" or "cox". Currently only these families are implemented. In the future others
#' will be added.
#' @param weights optional vector of non-negative weights, default is NULL which results in all weights = 1.
#'
#' @param loo TRUE (the default) means that uniLasso uses the prevalidated loo fits (approximate loo or 'alo' for "binomial" and "cox") for each univariate model as features to avoid overfitting.
#' \code{loo=FALSE} means it uses the univariate fitted predictor.
#' @param lower.limits = 0 (default) means that uniLasso  constrains the sign of the coefs in the second round to be that of the univariate fits.
#' @param standardize input argument to  glmnet for final non-negative lasso fit. Strongly recommend \code{standardize=FALSE} (default) since the univariate fit determines the right scale for each variable.
#' @param info Users can supply results of \code{uniInfo} on external datasets rather than compute them on the same data used to fit the model. If this is supplied, its \code{$betas} are used. Default is NULL.
#' @param loob.nit Number of Newton iterations for GLM or Cox in computing univariate linear predictors. Default is 2.
#' @param loob.eps A small number used in regularizing the Hessian for the Cox model. Default is 0.0001.
#' @param \ldots additional arguments passed to \code{cv.glmnet}.
#' @return An object that inherits from class \code{"cv.glmnet"}. There is one additional parameter returned, which is `info` and has two components.
#' They are  \code{beta0} and \code{beta}, the intercepts and slopes for the usual (non-LOO) univariate fits from stage 1.
#'
#' @examples
#' # cv.uniLasso examples
#' # Default usage with Gaussian data
#'
#' sigma =3
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)*sigma
#' xtest=matrix(rnorm(n * p), n, p)
#' ytest <- xtest %*% beta + rnorm(n)*sigma
#'
#' cvfit <- cv.uniLasso(x, y)
#' plot(cvfit)
#' predict(cvfit,xtest[1:10,], s="lambda.min") # predict at some test data points
#'
#' # Two-stage variation where we carve off a small dataset for computing the univariate coefs.
#'
#' cset=1:20
#' info = uniInfo(x[cset,],y[cset])
#' cvfit_two_stage <- cv.uniLasso(x[-cset,], y[-cset], info = info)
#' plot(cvfit_two_stage)
#'
#' # Binomial response cv.uniLasso
#'
#' yb =as.numeric(y>0)
#' cvfitb = cv.uniLasso(x, yb, family="binomial")
#' predict(cvfitb, xtest[1:10,], type="response") # predict at default s = "lambda.1se"
#'
#'
#' # cv.uniLasso with same positivity constraints, but starting `beta`
#' # from univariate fits on the same data. With loo=FALSE, does not tend to do as well,
#' # probably due to overfitting.
#'
#'  cvfit_pos_adapt <- cv.uniLasso(x, y, loo = FALSE)
#'  plot(cvfit_pos_adapt)
#'
#' # cv.uniLasso with no constraints, but starting `beta` from univariate fits.
#' # This is a version of the adaptive lasso, which tends to overfit, and loses interpretability.
#'
#'  cvfit_adapt <- cv.uniLasso(x, y, loo = FALSE, lower.limits = -Inf)
#'  plot(cvfit_adapt)
#'
#' # Cox response cv.uniLasso
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
#' cvfitc = cv.uniLasso(x, y, family = "cox")
#' plot(cvfitc)
#'
#' @export


cv.uniLasso <- function(x,y,family=c("gaussian","binomial","cox"),weights=NULL,
                      loo=TRUE,
                      lower.limits=0,
                      standardize=FALSE,
                      info=NULL,
                      loob.nit=2,
                      loob.eps=0.0001,
                      ...){
    this.call = match.call()
    family=match.arg(family)
    if(is.null(info)){ # user did not supply info
        info = uniInfo(x,y,family,weights,loob.nit,loob.eps,loo)
    }
    else {
        if(!is.null(info$F))warning("You supplied info with a loo 'F' component; we ignore that, and use '$beta' and'$beta0' instead.")
        loo=FALSE # we cannot trust the supplied info to give the right number of rows
        }
    if(loo)
        xp=info$F
    else {
        ones=rep(1,nrow(x))
        xp=x*outer(ones,info$beta)+outer(ones,info$beta0)
    }
    dimnames(xp)=dimnames(x)
    fit = cv.glmnet(xp,y,weights=weights,
                    lower.limits=lower.limits,
                    family=family,standardize=standardize,...)
    gfit=fit$glmnet.fit
    offset=gfit$offset
    if(offset)gfit$offset=FALSE # temporarily disable offset for intercept calculation
    a0=drop(predict(gfit,info$beta0))
    if(offset)gfit$offset=TRUE
    gfit$beta=gfit$beta*outer(info$beta,rep(1,length(a0)))
    gfit$a0=a0
    fit$glmnet.fit=gfit
    fit$info = info[c("beta0","beta")]
    class(fit)=c("cv.uniLasso",class(fit))
    fit$call = this.call
    fit
    }




