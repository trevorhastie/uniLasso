#' Create the univariate info for use in uniLasso
#'
#' Fit p separate univariate fits,  and if requested computes the loo fit matrix F.
#' It is called internally by \code{uniLasso}, or can be called externally on separate data and passed as input to \code{uniLasso}.
#' Currently this function can accommodate "gaussian", "binomial", and "Cox" families.
#'
#' @param X An n x p feature matrix
#' @param y A response object, depending on the family. For "gaussian" it is just a response vector. For "binomial" either a binary vector, a two level factor, or a two column non-negative matrix with rows summing to 1. For "cox" it is a Surv object (currently for right censored data).
#' @param family one of "gaussian","binomial" or "cox". Currently only these families are implemented. In the future others
#' will be added.
#' @param weights Vector of non-negative weights. Default is NULL, which results in all weights equal to 1.
#' @param nit Number of iterations if Newton steps are required (in "binomial" and "cox"). Default is 2. In principal more is better, but in some cases can run into convergence issues.
#' @param eps A small number to regularize the hessian for "cox"; default is 0.0001.
#' @param loo A logical, default=FALSE. If TRUE it computes the matrix of loo fits F.
#' @return an list with components \code{$beta} and \code{$beta0}, and if \code{loo=TRUE}, a n x p matrix \code{F} with the loo fits.
#' @examples
#' # Gaussian model
#' set.seed(1)
#' sigma=3
#' n <- 100; p <- 20
#' x <- matrix(rnorm(n * p), n, p)
#' beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
#' y <- x %*% beta + rnorm(n)*sigma
#'
#'
#' info = uniInfo(x,y)
#' names(info)
#'
#' yb = as.numeric(y>0)
#' info = uniInfo(x,yb, family = "binomial", loo = TRUE)
#' names(info)
#' @export

uniInfo = function(X,y,
                   family=c("gaussian","binomial","cox"),
                   weights=NULL,
                   nit = 2,
                   eps = 0.0001,
                   loo = FALSE){
### Check for constant columns in X
    s = apply(X,2,sd)
    zerosd = s==0
    if(all(zerosd))stop("All features are constant")
    if(any(zerosd))X = X[,!zerosd,drop=FALSE]
    if(is.null(weights))weights=rep(1,nrow(X))
    family=match.arg(family)
    info=switch(family,
                "gaussian" = loob_ols(X,y,w0=weights,loo=loo),
                "binomial" = loob_bin(X,y,w0=weights,nit=nit,loo=loo),
                "cox" = loob_cox(X,y,w0=weights,nit=nit,eps=eps,loo=loo),
                stop("Family not yet implemented")
                )
    if(any(zerosd))info = zerofix(info,zerosd)
    info
}
loob_ols = function(X,y,w0, loo=FALSE){
    ## LOO calculations for Gaussian model
    ## X is n x p model matrix, y is n-vector response
    ## w0 is weight vector, non-negative entries
    ## Returns the weighted univariate regression coefficients for each column of X
    ## if loo=TRUE, it also returns F the weighted prevalidated fit matrix (one at a time)
    n=length(y)
    w0 = n*w0/sum(w0)
    y=as.vector(y)
    p = ncol(X)
    W = outer(w0,rep(1,p))
    Z = outer(y,rep(1,p))
    wob = wlsu(X,W,Z)
    out = list()
    if(loo){
        Ws = sqrt(W)
        ones=rep(1,n)
        Xs=Ws*with(wob,scale(X,xbar,s))
        Ri = with(wob, n*(Ws*(Z-outer(ones,beta0))-Xs*outer(ones,beta))/(n-Ws^2-Xs^2))
        F = Z-Ri/Ws
        F[w0==0,] <- 0
        out$F=F
    }
    with(wob, c(out,list( beta=beta/s, beta0 = beta0-xbar*beta/s)))
    }
loob_bin = function(X,y,w0,nit=2, loo=FALSE){
    ## LOO calculations for Binomial model
    ## X is n x p model matrix
    ## y is either a two-column matrix, a binary vector, or a two level factor
    ## w0 is weight vector, non-negative entries
    ## Returns F the prevalidated fit matrix (one at a time)
    ## also the univariate coefficients
    n = nrow(X)
    p = ncol(X)
    ## check on y
    yout=multY(y)
    y=yout$y
    nc=dim(y)
    noo=nc[1]
    if(noo!=n)stop("X and y have different number of rows in call to uniInfo",call.=FALSE)
    nc=as.integer(nc[2])
    if(nc!=2)stop("More than two classes detected")
    y=y[,2]
    if(!is.null(yout$wts))w0=w0*yout$wts
    w0 = n*w0/sum(w0)
    W0 = outer(w0,rep(1,p))
    ## Initialization
    mus=(y+.5)/2
    w=mus*(1-mus)
    etas=log(mus/(1-mus))
    z = etas +(y-mus)/w
    W =outer(w,rep(1,p))
    Z=outer(z,rep(1,p))
   wob = wlsu(X,W*W0,Z)
    iter=1
    while(iter < nit){
        iter=iter+1
        Mus = 1/(1 + exp(-wob$Eta))
        W=Mus*(1-Mus)
        Z = wob$Eta + (outer(y,rep(1,p))-Mus)/W
        wob = wlsu(X,W*W0,Z)
    }
    out=list()
    if(loo){
        Ws=W*W0
        Ws = sqrt(scale(Ws,FALSE,colSums(Ws)/n))
        ones=rep(1,n)
        Xs=Ws*with(wob,scale(X,xbar,s))
        Ri = with(wob, n*(Ws*(Z-outer(ones,beta0))-Xs*outer(ones,beta))/(n-Ws^2-Xs^2))
        F = Z-Ri/Ws
        F[w0==0,] <- wob$Eta[w0==0,]
        isna=is.na(F)
        if(any(isna)) F[isna]=wob$Eta[isna]
        out$F=F
    }
    with(wob, c(out,list( beta=beta/s, beta0 = beta0-xbar*beta/s)))
}

loob_cox = function(X,y,w0,nit=4,eps=0.0001,loo=FALSE){
    ## LOO calculations for Cox PH survival model
    ## X is n x p model matrix, y is Surv  object as expected by glmnet
    ## Currently we do right censored, so y should have first column time and second column status
    ## In addition, we handle ties using Breslow method
    ## Returns F the prevalidated fit matrix (one at a time)
    ## also the univariate coefficients

    p = ncol(X)
    time <- y[, 1]
    d    <- y[, 2]
    n=length(time)
    w0 = n*w0/sum(w0)
    W0 = outer(w0,rep(1,p))

    ## Initialization
    Eta=matrix(0,n,p)
    gradob=coxgradu(Eta, time,d)
    o = gradob$o
    W = pmax(-gradob$diag_hessian,eps)
    Z = Eta + gradob$grad/W
    wlsu_ni = function(X,W,Z){
        beta= colSums(X*W*Z)/colSums(X*X*W)
        Eta = X*outer(rep(1,n),beta)
        list(beta=beta,Eta=Eta)
    }
    wob = wlsu_ni(X,W*W0,Z)
    iter=1
    while(iter < nit){
        iter=iter+1
        gradob = coxgradu(wob$Eta,time,d,o=o)
        W = pmax(-gradob$diag_hessian,eps)
        Z = wob$Eta + gradob$grad/W
        wob = wlsu_ni(X,W*W0,Z)
    }
    out=list()
    if(loo){
        X2w = X*X*W*W0
        X2w = scale(X2w,FALSE,colSums(X2w))
        Ri = (Z-X*outer(rep(1,n),wob$beta))/(1-X2w)
        F=Z-Ri
        F[w0==0,] <- wob$Eta[w0==0,]
        isna=is.na(F)
        if(any(isna)) F[isna]=wob$Eta[isna]
        out$F=F
    }
   c(out,list(beta=wob$beta, beta0 = rep(0,p)))
}


wlsu <- function(X,W,Z){
    n=nrow(X)
    totW = colSums(W)
    xbar = colSums(W*X)/totW
    Xm = X-outer(rep(1,n),xbar)
    s = sqrt(colSums(W*Xm^2)/totW)
    Xs = scale(Xm,FALSE,s)
    beta= colSums(Xs*W*Z)/totW
    beta0 = colSums(W*Z)/totW
    Eta = outer(rep(1,n),beta0) + Xs*outer(rep(1,n),beta)
    list(beta=beta,beta0=beta0,Eta=Eta,xbar=xbar,s=s)
}

coxgradu <- function(eta, time,d, w, o){
    ## eta is a n x p matrix of univariate cox fits
    ## o is an order vector, which is included in the result
    nobs = length(time)
    if (missing(w)) w=rep(1,nobs)
    w=w/sum(w)
    eta <- scale(eta, TRUE, FALSE)  # center eta so exponents are not too large

        # order exp(eta), time, d and w in ascending time order
        # for tied times, all deaths come before censored observations
        if(missing(o)) o <- order(time, d, decreasing = c(FALSE, TRUE))
        exp_eta <- exp(eta)[o,]
        time <- time[o]
        d <- d[o]
        w <- w[o]
    rskden <- apply((exp_eta*w)[nobs:1,],2,cumsum)[nobs:1,]
    ##reverse order inside;last guy is in all the risk sets

    ## for now we will rerun it each time
        ### See if there are dups in death times
        dups <- fid(time[d == 1],seq(length(d))[d == 1])
        dd <- d
        ww <- w

        ### next code replaces each sequence of tied death indicators by a new
        ### sequence where only the first is a 1 and the rest are zero. This
        ### makes the accounting in the following step work properly we also
        ### sums the weights in each of the tied death sets, and assign that
        ### weight to the first
        if(!is.null(ties<-dups$index_ties)){
            dd[unlist(ties)]=0
            dd[dups$index_first]=1
            wsum=sapply(ties,function(i,w)sum(w[i]),ww)
            tie1=sapply(ties,function(i)i[1])
            ww[tie1]=wsum
        }

        ### Get counts over risk sets at each death time
    rskcount=cumsum(dd)#this says how many of the risk sets each observation is in; 0 is none
### End of code that can be cached

        ### We now form partial sums of the 1/den just at the risk sets
        rskdeninv=apply((ww/rskden)[dd==1,], 2, cumsum)
        ### pad with a zero, so we can index it
        rskdeninv=rbind(0,rskdeninv)

        ### compute gradient for each obs
        grad <- w * (d - exp_eta * rskdeninv[rskcount+1,])
        grad[o,] <- grad

        # Now for the diag of Hessian
            rskdeninv2 <- apply((ww/(rskden^2))[dd==1,],2,cumsum)
            rskdeninv2 <- rbind(0, rskdeninv2)
            w_exp_eta <- w * exp_eta
            diag_hessian <- w_exp_eta^2 * rskdeninv2[rskcount+1,] - w_exp_eta * rskdeninv[rskcount+1,]
            diag_hessian[o,] <- diag_hessian
    list(grad=grad,diag_hessian=diag_hessian, o=o)
    }
fid <- function(x,index) {
    idup=duplicated(x)
    if(!any(idup)) list(index_first=index,index_ties=NULL)
    else {
        ndup=!idup
        xu=x[ndup]# first death times
        index_first=index[ndup]
        ities=match(x,xu)
        index_ties=split(index,ities)
        nties=sapply(index_ties,length)
        list(index_first=index_first,index_ties=index_ties[nties>1])
    }
}

zerofix <- function(info, zerosd){
    beta = beta0 = double(length(zerosd))
    beta[!zerosd] = info$beta
    beta0[!zerosd] = info$beta0
    infonew=list(beta=beta,beta0=beta0)
    F = info$F
    if(!is.null(F)){
        Fnew=matrix(0,nrow(F),length(zerosd))
        Fnew[,!zerosd] = F
        infonew$F=Fnew
    }
    infonew
}

multY <- function(y,minclass=0,warnclass=4){
    ## Utility function for binomial and multinomial
    ## This function takes y and produces a Y matrix
    ## If it already provided a y matrix, it ensures that it is of the right kind,
    ## and extracts a weight vector if appropriate.
    nc=dim(y)
if(is.null(nc)||nc[2]==1){
    ## Need to construct a y matrix, and include the weights
        y=as.factor(as.vector(y))
        ntab=table(y)
        nc=as.integer(length(ntab))
        if(nc<=1)stop("only 1 unique value in y")
        mincl=min(ntab)
        if(mincl<=minclass)stop(paste0("one binomial/multinomial class has ",minclass," or less observations; not allowed"))
        if(mincl<=warnclass)warning(paste0("one binomial/multinomial class has ",warnclass," or less observations; dangerous ground"))
        list(y=diag(nc)[as.numeric(y),],wts=NULL)
}else{
    wt = rowSums(y)
    if(any(wt!=1)){
        wtp=wt
        wtp[wt==0]=1
        y=y/wtp
    }
    else wt=NULL
    list(y=y,wts=wt)
}
}



