#' Simulate data for use in uniLasso and uniReg
#' We use some standard examples in our uniLasso paper, and for convenience we provide generators for these datasets.
#' @param example which of the prepackaged examples to use. Choices are "low-SNR","medium-SNR","high-SNR","home-court","two-class","counter-example", as described in the uniLasso paper. The three SNRs used are 0.5 (low), 1.0 (medium) and 2.0 (high) (also used for home-court). The training sizes for the first four are 300, and test sizes 3000.
#' @param wide logical variable which determines if p>n (default, 1000) or not (100).
#' This function calls worker functions simulate_gaussian(), simulate_two-class(), and simulate_counterexample(), which are currently not documented.
#' @return a list with components "x", "y", "xtest", "ytest", "mutest", and "sigma", where "mutest" is the true test mean, and "ytest <- mutest + rnorm(3000)*sigma."
#' @examples
#' dat = simulate_uniLasso("high-SNR")
#' fit = cv.uniLasso(dat$x, dat$y)
#' mse = mean( (predict(fit, dat$xtest)- dat$mutest)^2)
#' @export

simulate_uniLasso <- function(
                  example=c("low-SNR","medium-SNR","high-SNR","home-court","two-class","counter-example"),wide=TRUE){
    example=match.arg(example)
    traintest_split = function(d,ntrain){
        itrain=seq(ntrain)
        list(x=d$x[itrain,],y=d$y[itrain],xtest=d$x[-itrain,],ytest=d$y[-itrain],
             mutest=d$mu[-itrain],sigma=d$sigma)
    }
    p = ifelse(wide,1000,100)
    switch(example,
           "low-SNR" = traintest_split(
               simulate_Gaussian(3300,p,snr=0.5),
               300),
           "medium-SNR" = traintest_split(
               simulate_Gaussian(3300,p,snr=1),
               300),
           "high-SNR" = traintest_split(
               simulate_Gaussian(3300,p,snr=2),
               300),
           "home-court" = traintest_split(
               simulate_Gaussian(3300,p,snr=2,homecourt=TRUE),
               300),
           "two-class" =  traintest_split(
               simulate_twoclass(2200,wide),
               200),
           "counter-example" = traintest_split(
               simulate_counterexample(1100),
               100)
           )
    }
simulate_Gaussian=function(n=300,
                           p = 1000,
                           snr=1,
                           rho=0.8,
                           sparsity=0.1,
                           homecourt=FALSE){
    nonzero = round(sparsity * p)
    if(homecourt){
   ##Covariance matrix for AR(1)
        Sigma=rho^abs(outer(1:p,1:p,"-"))
        x = mvrnorm(n,mu=rep(0,p),Sigma)
        nonzero = pmin(nonzero,floor(p/2))
        nonzero_indices = seq(from=1,by=2,length=nonzero)
        beta=double(p)
        beta[nonzero_indices] <- runif(nonzero, 0.5, 2.0)
    } else{
        u=rnorm(n)
        fac=sample(c(-1,1),size=p,replace=TRUE)
        x=matrix(rnorm(n*p),n,p)+scale(matrix(u,n,p),FALSE,fac)
        beta=c(rnorm(nonzero),rep(0,p-nonzero))
        }
    mu=x%*%beta
    sigma = drop(sqrt(var(mu)/snr))
    y=mu+sigma*rnorm(n)
    list(x=x,y=y,mu=mu,sigma=sigma)
}
simulate_twoclass=function(n,wide){
    p=ifelse(wide,500,100)
        data = simulate_Gaussian(n=n,p=p,homecourt=TRUE)
        y = sample(c(0,1),n,replace=TRUE)
        x=data$x
        x[y==1,1:20]=x[y==1,1:20]+0.5
        data$x=x;data$y=y
        data
    }
   simulate_counterexample = function(n){
        p=20
        sigma=0.5
        x=matrix(rnorm(n*p),n,p)
        x[,2]=x[,2]+x[,1]
        beta=c(1,-.5,rep(0,p-2))
        mu=x%*%beta
        y=mu+rnorm(n)*sigma
        list(x=x,y=y,mu=mu,sigma=sigma)
        }

