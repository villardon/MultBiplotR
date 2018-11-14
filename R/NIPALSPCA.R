NIPALSPCA <- function(X, dimens=2,  tol=0.000001, maxiter=1000){
  n = dim(X)[1]
  p = dim(X)[2]
  U= matrix(0, n, dimens)
  V= matrix(0, p, dimens)
  d=matrix(0, dimens, 1)
  E=X
  for (i in 1:dimens){
    u=as.matrix(E[,1])
    # u=u/sqrt(sum(u^2))
    error=1
    iter=0
    while ((error>tol) & (iter<maxiter)){
      iter=iter+1
      uold=u
      v=matrix(t(u) %*% E, ncol=1)/sum(u^2)
      v=v/sqrt(sum(v^2))
      u=E %*% v
      # u=u/sqrt(sum(u^2))
      error=sum((u-uold)^2)
    }
    d[i]=sqrt(sum(u^2))
    E=E-u %*% t(v)
    u=u/d[i]
    V[,i]=v
    U[,i]=u
  }
  U=X%*%V
  d=sqrt(apply(U^2,2,sum))
  result=list(d=as.vector(d), u=U, v=V)
  return(result)
}


