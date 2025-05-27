#' Compute bootstrap confidence intervals for a univariate guided regression model
#'
#'
#' @param B Number of bootstrap samples. Default is 500.
#' @param alpha size of confidence interval.
#' @examples
#' # ci.uniReg usage
#'
#' sigma =3
#' set.seed(1)
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)*sigma
#' ci <- ci.uniReg(x, y, B=100)
#' print(ci)
#'
#' @rdname uniLasso
#'
#' @export ci.uniReg


ci.uniReg <- function(x,y,family=c("gaussian","binomial","cox"),weights=NULL,
                      B=500,
                      alpha=0.05,
                      ...){
    this.call=match.call()
    fit = uniReg(x,y,weights,...)
    beta = coef(fit)
    vn = row.names(beta)[-1]
    p=length(vn)
    betamat = matrix(0,p,B)
    w0 = weights
    n = nrow(x)
    if(is.null(w0))w0=rep(1,n)
    bootstrap_weights <- function(n) {
        tabulate(sample(1:n, n, replace = TRUE), nbins = n)
}
    for(b in 1:B){
        wb = bootstrap_weights(n)*w0
        bfit = uniReg(x,y,family=family,weights=wb,...)
        betamat[,b] = as.numeric(coef(bfit))[-1]
    }
    ci = apply(betamat, 1, quantile, probs = c(alpha/2,1-alpha/2))
    ci = t(ci)
    row.names(ci)=vn
    ci = cbind(Coef=beta[-1,],ci)
    allneg = apply(ci<0,1,all)
    allpos = apply(ci>0,1,all)
    sig = allneg|allpos
    sigstar=character(p)
    if(any(sig))sigstar[sig]="*"
    data.frame(ci,"Sig."=sigstar)
    }




