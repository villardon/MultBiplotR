PathwiseCoordinateOptimization <- function(y, X, beta=NULL, lambda=0.1,  tol=0.000001, maxiter=1000, show=FALSE){
  n = dim(X)[1]
  p = dim(X)[2]
  if (is.null(beta)) beta=as.matrix(ginv(t(X) %*% X) %*% t(X) %*% y)
  error=sum(beta^2)
  iter=0
  while ((error>tol) & (iter<maxiter)){
    iter=iter+1
    betaold=beta
    for (j in 1:p){
      r = y - X[,-j] %*% beta[-j,]
      b= sum(r * X[,j])/n
      beta[j,]=0
      if (b >0 & lambda < abs(b)) beta[j,]=b-lambda
      if (b <0 & lambda < abs(b)) beta[j,]=b+lambda
    }
    error=sum((beta-betaold)^2)
    if (show)
      print(c(iter,error))
  }
  return(beta)
}

