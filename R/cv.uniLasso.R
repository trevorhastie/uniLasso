#' Fit a uniLasso model.
#' This function fits a univariate-guided sparse regression (lasso)  using cross-validation to select some tuning parameters. Conveniently,  it returns an objects that inherits from \code{cv.glmnet}, so that
#' all of the methods for \code{cv.glmnet} can be applied, such as predict, plot, coef, print,
#'  and assess.glmnet
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
#'
#' @param loo TRUE means (default) that `cv.uniLasso` uses the prevalidated loo fits (alo for binomial and cox) for each univariate model as features  to avoid overfitting.  \code{loo=FALSE} means it uses the univariate fitted predictor (the default \code{loo=TRUE}.
#'   FALSE means it uses the univariate fitted predictor.
#' @param lower.limits = 0 (default) means that uniLasso  constrains the sign of the coefs in the second round to be that of the univariate fits
#' @param standardize input argument to  glmnet for final non-negative lasso fit. Strongly recommend \code{standardize=FALSE} (default) since the univariate fit determins the right scale for each variable.
#' @param info Users can supply results of \code{uniInfo} on external datasets rather than compute them on the same data used to fit the model. If this is supplied, its \code{$betas} are used. Default is NULL.
#' @param loob.nit Number of Newton iterations for GLM or Cox in computing univariate linear predictors. Default is 4.
#' @param loob.eps A small number used in regularizing the Hessian for the Cox model. Default is 0.0001.
#' @param \ldots additional arguments passed to \code{cv.glmnet}.
#' @return An object that inherits from class \code{"cv.glmnet"}.
#'  @param info This has two components, \code{beta0} and \code{beta}, - the intercepts and slopes for the usual (non-LOO) univariate fits from stage 1.
#'
#' Fits a two stage lasso model. First stage replaces each feature by the univariate fit for that feature.
#'   Second stage fits a (positive) lasso using the first stage features. Hence the second stage selects and
#'   modifies the coefficents of the first stage model, similar to the adaptive lasso. Leads to potentially sparser models.
#'
#'  for "binomial" family y is a binary response
#'  for "cox" family y should be a Surv object for right censored data
#'  or a matrix with columns labeled 'time' and 'status'
#'  Although glmnet has more flexible options say for binary responses, and for cox
#'  responses, these are not yet implemented but are possible.
#'  Likewise, other glm families are possible as well, but not yet implemented.
#'
#'  This is a one-visit function that returns a \code{'cv.glmnet'} object.
#'  You can make predictions from the whole path, at 'lambda.min' etc just like you can for
#'  a `'cv.glmnet object'`.
#'
#'  loo = TRUE means it uses the prevalidated loo fits (alo for binomial and cox)for each univariate model as features
#'  loo = FALSE means it uses the univariate fitted predictor and hence it is a form of adaptive lasso.
#'  lower.limits = 0 means it constrains the sign of the coefs in the second round to be that of the univariate fits.
#'
#' @examples
#' # Gaussian model
#' set.seed(1)
#' sigma = 3
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)*sigma
#' xtest=matrix(rnorm(n * p), n, p)
#' ytest <- xtest %*% beta + rnorm(n)*sigma
#'
#' # Default usage
#'
#' fit <- cv.uniLasso(x, y)
#' plot(fit)  #plot cv curve
#' yhat<- predict(fit,xtest) #predict on test data
#' beta=coef(fit,s="lambda.min") #extract coefs of final model at value of lambda minimizing CV error.
#'
#' # Adpative lasso with same positivity constraints, but starting `beta` from univariate fits.
#'
#'  fit_pos_adapt <- cv.uniLasso(x, y, loo = FALSE)
#'  plot(fit_pos_adapt)
#'
#' # Adpative lasso with no constraints, but starting `beta` from univariate fits.
#'
#'  fit_adapt <- cv.uniLasso(x, y, loo = FALSE, lower.limits = -Inf)
#'  plot(fit_adapt)
#'
#' # Adaptive lasso where we carve off a small dataset for computing the univariate coefs.
#'
#' cset=1:20
#' info = uniInfo(x[cset,],y[cset])
#' fit_two_stage <- cv.uniLasso(x[-cset,], y[-cset], info = info)
#' plot(fit_two_stage)
#' @export


cv.uniLasso <- function(x,y,family=c("gaussian","binomial","cox"),
                      loo=TRUE,
                      lower.limits=0,
                      standardize=FALSE,
                      info=NULL,
                      loob.nit=4,
                      loob.eps=0.0001,
                      ...){

    family=match.arg(family)
    require(glmnet)
    if(is.null(info)){ # user did not supply info
        info = uniInfo(x,y,family,loob.nit,loob.eps,loo)
    }
    else {
        loo=FALSE # we cannot trust the supplied info to give the right number of rows
        }
    if(loo)
        xp=info$F
    else {
        ones=rep(1,nrow(x))
        xp=x*outer(ones,info$beta)+outer(ones,info$beta0)
        }
    fit = cv.glmnet(xp,y,
                    lower.limits=lower.limits,
                    family=family,standardize=standardize,...)
    gfit=fit$glmnet.fit
    a0=drop(predict(gfit,info$beta0))
    gfit$beta=gfit$beta*outer(info$beta,rep(1,length(a0)))
    gfit$a0=a0
    fit$glmnet.fit=gfit
    fit$info = info[c("beta0","beta")]
    class(fit)=c("cv.uniLasso",class(fit))
    fit
    }




