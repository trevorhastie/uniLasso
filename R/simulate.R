#' Simulate data for use in uniLasso and uniReg
#'
#' We use some standard examples in our uniLasso paper, and for convenience we provide generators for these datasets.
#' @param example which of the prepackaged examples to use. Choices are "low-SNR","medium-SNR","high-SNR","home-court","two-class","counter-example", as described in the uniLasso paper. The three SNRs used are 0.5 (low), 1.0 (medium) and 2.0 (high) (also used for home-court). The training sizes for the first four are 300, and test sizes 3000.
#' @param wide logical variable which determines if p>n (default, 1000) or not (100).
#' This function calls worker functions simulate_gaussian(), simulate_two-class(), and simulate_counterexample(), which are currently not documented.
#' @return a list with components "x", "y", "xtest", "ytest", "mutest", and "sigma", where "mutest" is the true test mean, and "ytest <- mutest + rnorm(nrow(xtest))*sigma."
#' @examples
#' dat = simulate_uniLasso("high-SNR")
#' fit = cv.uniLasso(dat$x, dat$y)
#' mse = mean( (predict(fit, dat$xtest)- dat$mutest)^2)
#' @export

simulate_uniLasso <- function(
                  example=c("low-SNR","medium-SNR","high-SNR","home-court","two-class","counter-example"),wide=TRUE){
    example=match.arg(example)
    p = ifelse(wide,1000,100)
    switch(example,
           "low-SNR" =    simulate_Gaussian(300,3000,p,snr=0.5),
           "medium-SNR" = simulate_Gaussian(300,3000,p,snr=1),
           "high-SNR" =   simulate_Gaussian(300,3000,p,snr=2),
           "home-court" = simulate_Gaussian(300,3000,p,snr=2,homecourt=TRUE),
           "two-class" =  simulate_twoclass(200,2000,wide),
           "counter-example" = simulate_counterexample(100, 1000)
           )
}

#' simulate Gaussian data
#'
#' A simulator that builds a training and test set with particular characteristics, as used in our "uniLasso" paper.
#' @param ntrain number of training examples.
#' @param ntest number of test examples.
#' @param p number of features.
#' @param snr desired SNR (signal-to-noise ratio).
#' @param rho for \code{homecourt=TRUE} 'rho' controls the autocorrelation between variables. Variables k units apart have correlation \code{rho^k}.
#' @param sparsity fraction of variables with nonzero coefficients.
#' @param homecourt logical; if \code{TRUE} then correlated features, with a special boost for large coefficients, mimicking the uniLasso two-stage algorithm.
#' @return a list with components "x", "y", "xtest", "ytest", "mutest", and "sigma", where "mutest" is the true test mean, and "ytest <- mutest + rnorm(ntest)*sigma."
#' @examples
#' dat = simulate_Gaussian(300,3000,p=500,snr=1.2)
#' fit = cv.uniLasso(dat$x, dat$y)
#' mse = mean( (predict(fit, dat$xtest)- dat$mutest)^2)
#' @export

simulate_Gaussian=function(ntrain = 300,
                           ntest = 3000,
                           p = 1000,
                           snr=1,
                           rho=0.8,
                           sparsity=0.1,
                           homecourt=FALSE){
    nonzero = round(sparsity * p)
    n = ntrain+ntest
    if(homecourt){
   ##Covariance matrix for AR(1)
        Sigma=rho^abs(outer(1:p,1:p,"-"))
        x = mvrnorm(n,mu=rep(0,p),Sigma)
        nonzero = pmin(nonzero,floor(p/2))
        nonzero_indices = seq(from = 1,to=nonzero)
        beta = double(p)
        beta[nonzero_indices] <- runif(nonzero, 0.5, 2)
        mu = x %*% beta
        sdmu = sd(mu)
        mu = mu/sd(mu)
        beta=beta/sdmu
        sigma = sqrt(1/snr)
        y = mu + sigma * rnorm(n)
        info = uniInfo(x,y)
        beta=beta*abs(info$beta)
    } else{
        u=rnorm(n)
        fac=sample(c(-1,1),size=p,replace=TRUE)
        x=matrix(rnorm(n*p),n,p)+scale(matrix(u,n,p),FALSE,fac)
        beta=c(rnorm(nonzero),rep(0,p-nonzero))
        }
    mu=x%*%beta
    sdmu = sd(mu)
    mu = mu/sd(mu)
    sigma = sqrt(1/snr)
    y = mu + sigma * rnorm(n)
    d = list(x=x,y=y,mu=mu,sigma=sigma)
    traintest_split(d, ntrain=ntrain)
}
#' simulate two class data
#'
#' @param ntrain number of training examples.
#' @param ntest number of test examples.
#' @param wide logical. If TRUE \code{p=500}, else \code{p=100}.
#' @return a list with components "x", "y", "xtest", "ytest", "mutest", and "sigma", where "mutest" is the true test mean, and "ytest <- mutest + rnorm(ntest)*sigma."
#' @examples
#' dat = simulate_twoclass(300,3000)
#' fit = cv.uniLasso(dat$x, dat$y, family="binomial")
#' misclass = mean( sign(predict(fit, dat$xtest,s="lambda.min"))== sign(dat$ytest-0.5))
#' @export
simulate_twoclass=function(ntrain,ntest,wide=TRUE){
    p=ifelse(wide,500,100)
    n = ntrain+ntest
        data = simulate_Gaussian(ntrain,ntest,p=p,homecourt=TRUE)
        y = sample(c(0,1),n,replace=TRUE)
        x=with(data,rbind(x,xtest))
        x[y==1,1:20]=x[y==1,1:20]+0.5
        data$x=x;data$y=y
        traintest_split(data, ntrain=ntrain)
    }
#' simulate counterexample data
#'
#' A particular counterexample where the first two features are strongly positively correlated,
#' yet they have coefficients of opposite sign in a multiple regression.
#' @param ntrain number of training examples.
#' @param ntest number of test examples.
#' @return a list with components "x", "y", "xtest", "ytest", "mutest", and "sigma", where "mutest" is the true test mean, and "ytest <- mutest + rnorm(ntest)*sigma."
#' @examples
#' dat = simulate_counterexample(300,3000)
#' fit = cv.uniLasso(dat$x, dat$y)
#' err = mean( (predict(fit, dat$xtest,s="lambda.min")- dat$mutest)^2)
#' @export
   simulate_counterexample = function(ntrain,ntest){
       p=20
       n = ntrain+ntest
        sigma=0.5
        x=matrix(rnorm(n*p),n,p)
        x[,2]=x[,2]+x[,1]
        beta=c(1,-.5,rep(0,p-2))
        mu=x%*%beta
        y=mu+rnorm(n)*sigma
       d = list(x=x,y=y,mu=mu,sigma=sigma)
       traintest_split(d,ntrain=ntrain)
   }
    traintest_split = function(d,ntrain){
        itrain=seq(ntrain)
        list(x=d$x[itrain,],y=d$y[itrain],xtest=d$x[-itrain,],ytest=d$y[-itrain],
             mutest=d$mu[-itrain],sigma=d$sigma)
    }


