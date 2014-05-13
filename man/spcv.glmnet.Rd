\name{spcv.glmnet}
\alias{spcv.glmnet}
\title{Spatial cross-validation for glmnet}
\usage{
spcv.glmnet(x, y, coords, dist_thres, longlat, weights, offset = NULL,
  lambda = NULL, type.measure = c("mse", "deviance", "class", "auc", "mae"),
  nfolds = 10, training_rows, grouped = TRUE, keep = FALSE,
  parallel = FALSE, ...)
}
\arguments{
  \item{x}{x matrix as in glmnet}

  \item{y}{response y as in glmnet.}

  \item{coords}{the spatial coordinates of the samples}

  \item{dist_thres}{the spatial distance threshold}

  \item{longlat}{boolean, if TRUE great circle distances
  are computed rather than euclidean distances}

  \item{...}{other variables that can be supplied to
  cv.glmnet}
}
\description{
Does leave-one-out cross-validation for glmnet, produces a
plot, and returns a value for lambda
}
\examples{
x = matrix(rnorm(100*20),100,20)
y = rnorm(100)
coords = expand.grid(1:10, 1:10)
nsp_cv = cv.glmnet(x,y)
sp_cv = spcv.glmnet(x, y, coords, 5, FALSE)
par(mfrow=c(1,2))
plot(nsp_cv)
plot(sp_cv)
}

