LASSORegression <- function(y, X, nlambdas=100, lambdas=NULL,  tol=0.00001, maxiter=100, show=FALSE){
  result=list()
  result$Means=apply(X, 2, mean)
  result$SD=apply(X, 2, sd)
  n = dim(X)[1]
  p = dim(X)[2]
  if (is.numeric(y)) {y= as.matrix(y)}
  if (is.data.frame(X)) X=DataFrame2Matrix4Regression(X, Intercept=FALSE)
  if (is.vector(X)) X=matrix(X, length(X),1)
  
  X=as.matrix(X)
  y=as.matrix(y)
  X=InitialTransform(X, transform = 5)$X 
  y=InitialTransform(y, transform = 5)$X 
  result$X=X
  result$y=y
  sct=sum(y^2)
  result$TSS=sct
  
  if (is.null(lambdas)){
    lambdamax= max((t(X) %*% y)/n)
    lambda.min.ratio = ifelse(n < p, 0.01, 1e-04)
    lambdamin = lambda.min.ratio * lambdamax
    lambdas=seq(lambdamin, lambdamax, length.out = nlambdas)
    lambdas=round(lambdas, digits=5)
  }

  result$Lambdas=lambdas
  nlambdas=length(lambdas)
  betas=matrix(0, p, nlambdas)
  beta=as.matrix(ginv(t(X) %*% X) %*% t(X) %*% y)
  for (j in 1:nlambdas){
    beta=PathwiseCoordinateOptimization(y, X, beta=beta, lambda=lambdas[j],  tol=tol, maxiter=maxiter, show=show)
    betas[,j]=beta
  }
  
  rownames(betas)=colnames(X)
  colnames(betas)=paste("l-",lambdas)
  result$Beta=betas
  fitted=X%*% betas
  resid=y %*% matrix(1, 1, nlambdas)-fitted
  sce=apply(fitted^2, 2, sum)
  scr=apply(resid^2, 2, sum)
  R2=round((1-scr/sct)*100, digits=2)
  result$ESS=sce
  result$RSS=scr
  result$R2=R2
  class(result)="LASSOReg"
  return(result)
  
}


