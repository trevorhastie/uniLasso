#' Fit a univariate guided regression model
#'
#' Fit a univariate-guided regression, by a two-stage procedure.
#' The first stage fits `p` separate univariate models to the response. The second stage fits a regression model, preserving the univariate signs.
#'
#' @param hard.zero if `TRUE` (default), the model fits the unpenalized regression. This may not be possible  when `p > n`. In this case `hard.zero = FALSE` is preferable, and the model is fit using the smallest value of `lambda` in the path.

#' @examples
#' # uniReg usage
#'
#' sigma =3
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)*sigma
#' xtest=matrix(rnorm(n * p), n, p)
#'
#' fit <- uniReg(x, y)
#' predict(fit,xtest[1:10,]) #predict on test data
#' coef(fit)
#' print(fit)
#'
#' fita <- uniReg(x, y, hard.zero = FALSE)
#' print(fita)
#'
#' fitb <- uniReg(x, y>0, family = "binomial")
#' coef(fitb)
#' print(fitb)
#'
#' @rdname uniLasso
#' @export uniReg


uniReg <- function(x,y,family=c("gaussian","binomial","cox"),
                   loo=TRUE,
                   lower.limits=0,
                   standardize=FALSE,
                   info=NULL,
                   loob.nit=2,
                   loob.eps=0.0001,
                   hard.zero = TRUE,
                      ...){
    this.call=match.call()
    family=match.arg(family)
    if(is.null(info)){ # user did not supply info
        info = uniInfo(x,y,family,loob.nit,loob.eps,loo)
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
    fit = glmnet(xp,y,
                    lower.limits=lower.limits,
                 family=family,standardize=standardize,...)
    lambda=fit$lambda
    if(hard.zero && min(lambda)>0){
        lambda = c(lambda,0)
        fit = glmnet(xp,y,lambda=lambda,
                    lower.limits=lower.limits,
                    family=family,standardize=standardize,...)
        }
    offset=fit$offset
    if(offset)fit$offset=FALSE # temporarily disable offset for intercept calculation
    a0=drop(predict(fit,info$beta0))
    if(offset)fit$offset=TRUE
    nlam=length(lambda)
    beta=(fit$beta*outer(info$beta,rep(1,length(a0))))[,nlam,drop=FALSE]
    colnames(beta)=NULL
    names(a0)=NULL
    fit$beta=beta
    fit$a0=a0[nlam]
    fit$lambda = fit$lambda[nlam]
    fit$dev.ratio=fit$dev.ratio[nlam]
    fit$df=fit$df[nlam]
    fit$info = info[c("beta0","beta")]
    fit$call = this.call
    class(fit)=c("uniReg","uniLasso",class(fit))
    fit
    }




