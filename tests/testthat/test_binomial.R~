context("Binomial")

sigma =3
set.seed(1)
n <- 100; p <- 20
x <- matrix(rnorm(n * p), n, p)
beta <- matrix(c(rep(2, 5), rep(0, 15)), ncol = 1)
y <- x %*% beta + rnorm(n)*sigma
xtest=matrix(rnorm(n * p), n, p)
ytest <- xtest %*% beta + rnorm(n)*sigma

# Binomial response uniLasso

yb =as.numeric(y>0)
fitb = uniLasso(x, yb, family="binomial")
predb = predict(fitb, xtest[1:10,], s=1, type="response")
cvfitb = cv.uniLasso(x, yb, family="binomial")
predcvb = predict(cvfitb, xtest[1:10,], type="response") # predict at default s = "lambda.1se"
fitub <- cv.uniReg(x, yb, family = "binomial")

objects = enlist(fitb,predb,cvfitb,predcvb,fitub)
###saveRDS(objects, "saved_results/test_binomial.RDS")

expected  <- readRDS("saved_results/test_binomial.RDS")
for (x in names(objects)) {
    cat(sprintf("Testing %s\n", x))
    expect_equal(objects[[x]], expected[[x]])
}

