# uniLasso (Univariate guided Lasso)
Fits a univariate-guided sparse regression (lasso), by a two-stage procedure. The first stage fits p separate univariate models to the response. The second stage gives more weight to the more important univariate features, and preserves their signs. Conveniently, it returns an objects that inherits from class `glmnet`, so that all of the methods for glmnet are available. See <doi:10.48550/arXiv.2501.18360> for details.

## Installation

To install the uniLasso R package directly from github, run the following in R:
```r
library(devtools)
install_github(repo="trevorhastie/uniLasso")
```
