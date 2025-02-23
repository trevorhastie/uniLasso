#' Compare the nonzero coefficients and univariate counterparts
#'
#' For a `cv.uniLasso` object, compare the CV-selected nonzero coefficints to their univariate counterparts. Also works for a `cv.glmne`t object.
#'
#' @param cv.object A `cv.uniLasso` or a `cv.glmnet` object.
#' @param info The result of a call to `uniInfo()`. If `cv.object` inherits from `cv.uniLasso`, the `$info` component will be used.
#' @param s the value of lambda to be used, with default `s="lambda.min"`. Alternatively, can be `s="lambda.1se".
#' @return a two columns matrix with one column being the non-zero coefficients from the `cv.object`, and the other being the corresponding univariate coefficients.
#'
#'  @examples
#'
#' sigma =3
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)*sigma
#'
#' cvfit <- cv.uniLasso(x, y)
#' uniCoef(cvfit)
#'
#' cvfit2 <- cv.glmnet(x,y)
#' uniCoef(cvfit2, info=cvfit$info)
#'
#' @export

uniCoef <- function(cv.object, info=NULL, s=c("lambda.min","lambda.1se"),...){
    type="uniLasso"
    s = match.arg(s)
    if(!inherits(cv.object,"cv.uniLasso"))type="Lasso"
    if(is.null(info))info = cv.object$info
    if(is.null(info))stop("Missing 'info'. Since you supplied a 'cv.glmnet' object, you must provide 'info', the result of a call to 'uniInfo'.")
    cfs = coef(cv.object,s=s,...)
    nz = cfs@i+1 #S4 class
    rns = row.names(cfs)
    if(rns[1]=="(Intercept)"){
        cfs=cfs[-1,]# drops to numeric
        rns=rns[-1]
        nz=nz[-1]-1
    }
    else cfs=cfs[,1]# drops to numeric
    mat=cbind(cfs[nz],info$beta[nz])
   dimnames(mat)=list(rns[nz],c(type,"Univariate"))
    mat
}
