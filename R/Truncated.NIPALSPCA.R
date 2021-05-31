Truncated.NIPALSPCA <- function(X, dimens=2,  tol=0.000001, maxiter=1000, lambda=0.02){
  n = dim(X)[1]
  p = dim(X)[2]
  d=matrix(0, dimens, 1)
  SD=NIPALSPCA(X, dimens=dimens)
  U=SD$u
  V=SD$v
  E=X
  for (i in 1:dimens){
    u=as.matrix(SD$u[,i])* SD$d[i]
    v=as.matrix(SD$v[,i])
    error=1
    iter=0
    while ((error>tol) & (iter<maxiter)){
      iter=iter+1
      uold=u
      v=matrix(t(u) %*% E, ncol=1)/sum(u^2)
      
      v=v/sqrt(sum(v^2))
      v=v * (abs(v) >lambda)
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