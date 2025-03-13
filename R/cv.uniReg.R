#' Fit a cross-validated univariate guided regression model.
#'
#' Fit a cross-validated univariate-guided sparse regression `uniLasso` model, with a focus on the end of the path which corresponds to the `uniReg` fit. Conveniently,  it returns an object that inherits from \code{cv.glmnet},and methods such as `predict`, `plot`, `coef`, `print` all gives sensible results.
#'
#' @param hard.zero if `TRUE` (default), the model includes `lambda=0` in the path. This may not be possible  when `p > n`. In this case `hard.zero = FALSE` is preferable, and the model is fit using the usual regularization path.

#' @examples
#' # cv.uniReg usage
#'
#' sigma =3
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)*sigma
#' xtest=matrix(rnorm(n * p), n, p)
#'
#'
#' fit <- cv.uniReg(x, y)
#' plot(fit)
#' coef(fit)
#' print(fit)
#' predict(fit,xtest[1:10,]) #predict on test data
#'
#' fita <- cv.uniReg(x, y, hard.zero = FALSE)
#' plot(fita)
#' print(fita)
#'
#' fitb <- cv.uniReg(x, y>0, family = "binomial")
#' plot(fitb)
#' print(fitb)
#'
#' @rdname cv.uniLasso
#' @export cv.uniReg



cv.uniReg <- function(x,y,family=c("gaussian","binomial","cox"),
                      loo=TRUE,
                      lower.limits=0,
                      standardize=FALSE,
                      info=NULL,
                      loob.nit=2,
                      loob.eps=0.0001,
                      hard.zero = TRUE,
                      ...){
    this.call = match.call()
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
    fit0 = glmnet(xp,y,
                    lower.limits=lower.limits,
                 family=family,standardize=standardize,...)
    lambda=fit0$lambda
    if(hard.zero && min(lambda)>0)lambda = c(lambda,0)
    fit = cv.glmnet(xp,y,
                    lower.limits=lower.limits,
                    family=family,standardize=standardize,lambda=lambda,...)
    gfit=fit$glmnet.fit
    offset=gfit$offset
    if(offset)gfit$offset=FALSE # temporarily disable offset for intercept calculation
    a0=drop(predict(gfit,info$beta0))
    if(offset)gfit$offset=TRUE
    gfit$beta=gfit$beta*outer(info$beta,rep(1,length(a0)))
    gfit$a0=a0
    fit$glmnet.fit=gfit
    fit$info = info[c("beta0","beta")]
    class(fit)=c("cv.uniReg","cv.uniLasso",class(fit))
    fit$call = this.call
    fit
    }



#' plot the cross-validation curve produced by cv.uniReg
#'
#' Plots the cross-validation `cv.uniLasso` curve (which is a `cv.glmnet` curve), and upper and lower standard deviation
#' curves, as a function of the \code{lambda} values used. It highlights the value at the end of the path, which is either `lambda = 0`, or else the smallest `lambda` if `hard.zero = FALSE`.
#'
#' A plot is produced, and nothing is returned.
#'
#' @param x fitted \code{"cv.uniReg"} object, which inherits from \code{"cv.uniLasso"}.
#' negative if \code{sign.lambda=-1}.
#' @param \dots Other graphical parameters to plot
#' @author Trevor Hastie and Rob Tibshirani\cr Maintainer:
#' Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{cv.uniLasso} and \code{glmnet:::cv.glmnet}.
#' @keywords models regression
#' @examples
#'
#' set.seed(1010)
#' n = 1000
#' p = 100
#' nzc = trunc(p/10)
#' x = matrix(rnorm(n * p), n, p)
#' beta = rnorm(nzc)
#' fx = (x[, seq(nzc)] %*% beta)
#' eps = rnorm(n) * 5
#' y = drop(fx + eps)
#' cvob0 = cv.uniReg(x, y)
#' plot(cvob0)
#' cvob = cv.uniReg(x, y, hard.zero = FALSE)
#' plot(cvob)
#'
#' @method plot cv.uniReg
#' @export
plot.cv.uniReg = function(x,...){
    lam=x$lambda
    nlam=length(lam)
    hard.zero=FALSE
    if(min(lam)==0){
        hard.zero=TRUE
        lam[nlam] = if(nlam>2)lam[nlam-1]^2/lam[nlam-2] else 0.001
        x$lambda <- lam
        }
    glmnet:::plot.cv.glmnet(x,...)
    points(log(lam[nlam]),x$cvm[nlam],pch=19,col="blue")
    leg = ifelse(hard.zero,expression(lambda==0),expression(lambda==smallest))
    legend("topleft",legend=leg,pch=19,col="blue")
}

#' print a cross-validated uniReg object
#'
#' Print a summary of the results of cross-validation for a uniReg model.
#'
#' A summary of the cross-validated uniReg fit is produced. This is an augmented summary of a cv.glmnet object, with an extra row corresponding to the smallest lambda in the path
#'
#' @param x fitted 'cv.uniReg' object
#' @param digits significant digits in printout
#' @param \dots additional print arguments
#' @author Trevor Hastie and Rob Tibshirani\cr Maintainer:
#' Trevor Hastie <hastie@@stanford.edu>
#' @examples
#'
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' fit1 = cv.uniReg(x, y)
#' print(fit1)
#' @method print cv.uniReg
#' @export
#' @export print.cv.uniReg
print.cv.uniReg <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall: ", deparse(x$call), "\n\n")

    optlams=c(x$lambda.min,x$lambda.1se)
    which=match(optlams,x$lambda)
    mat = with(x, cbind(optlams, which, cvm[which], cvsd[which], nzero[which]))
    dimnames(mat) = list(c("min", "1se"), c("Lambda", "Index","Measure",
                                            "SE", "Nonzero"))
    cat("Measure:", x$name,"\n\n")
    lambda=x$lambda
    nlam=length(lambda)

    mat = rbind(
        mat,
        "smallest"=with(x,c(lambda[nlam],nlam,cvm[nlam],cvsd[nlam],nzero[nlam])))
    mat=data.frame(mat,check.names=FALSE)
    class(mat)=c("anova",class(mat))
    print(mat,digits=digits)
}

#' make predictions from a "cv.uniReg" object.
#'
#' This function makes predictions from a cross-validated `uniReg` model, using
#' the stored \code{"glmnet.fit"} object, and by default the smallest value of `lambda` used.
#'
#' This function makes it easier to use the results of cross-validation to make
#' a prediction.
#'
#' @param object Fitted \code{"cv.uniReg"}.
#' @param newx Matrix of new values for \code{x} at which predictions are to be
#' made. Must be a matrix; can be sparse as in \code{Matrix} package. See
#' documentation for \code{predict.glmnet}.
#' @param s Value(s) of the penalty parameter \code{lambda} at which
#' predictions are required. Default is the value \code{s="zero"} which corresponds to the smallest value of `lambda` used.  Alternatively  \code{s="lambda.1se"} or \code{s="lambda.min"} can be used. If
#' \code{s} is numeric, it is taken as the value(s) of \code{lambda} to be
#' used.
#' @param \dots Not used. Other arguments to predict.
#' @return The object returned depends on the \dots{} argument which is passed
#' on to the \code{predict} method for \code{glmnet} objects.
#' @author Trevor Hastie and Rob Tibshirani\cr Maintainer:
#' Trevor Hastie <hastie@@stanford.edu>
#' @seealso \code{print}, and \code{coef} methods, and
#' \code{cv.uniReg}.
#' @keywords models regression
#' @examples
#'
#' x = matrix(rnorm(100 * 20), 100, 20)
#' y = rnorm(100)
#' cv.fit = cv.uniReg(x, y)
#' predict(cv.fit, newx = x[1:5, ])
#' coef(cv.fit)
#'
#' @method predict cv.uniReg
#' @export
predict.cv.uniReg=function(object,newx,s=c("zero","lambda.1se","lambda.min"),...){
  if(is.numeric(s))lambda=s
  else
      if(is.character(s)){
          s=match.arg(s)
          if(s!="zero")
              lambda=object[[s]]
          else
              lambda = min(object$lambda)
          names(lambda)=s
      }
    else stop("Invalid form for s")
  predict(object$glmnet.fit,newx,s=lambda,...)
}


#' @method coef cv.uniReg
#' @export
coef.cv.uniReg=function(object,s=c("zero","lambda.1se","lambda.min"),...){
  if(is.numeric(s))lambda=s
  else
      if(is.character(s)){
          s=match.arg(s)
          if(s!="zero")
              lambda=object[[s]]
          else
              lambda = min(object$lambda)
          names(lambda)=s
      }
    else stop("Invalid form for s")
  coef(object$glmnet.fit,s=lambda,...)
}
