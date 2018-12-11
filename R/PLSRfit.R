PLSRfit <- function(Y, X, S=2, tolerance=0.000005, maxiter=100, show=FALSE){
  I1=dim(X)[1]
  J=dim(X)[2]
  
  I2=dim(Y)[1]
  K=dim(Y)[2]
  I=I2
  
  T=matrix(0, I1, S)
  U=matrix(0, I1, S)
  W=matrix(0, J, S)
  C=matrix(0, K, S)
  P=matrix(0, J, S)
  Q=matrix(0, K, S)
  
  # Initial Step
  t0=matrix(1, I,1)
  c0=t(Y) %*% t0/ sum(t0^2)
  
  Y=Y-t0%*%t(c0)
  xb=matrix(apply(X,2,mean), ncol=1)
  X=X-t0%*%t(xb)
  
  for (i in 1:S){
    error=1
    iter=0
    u=matrix(X[,1],I1,1)
    while ((error>tolerance) & (iter<maxiter)){
      iter=iter+1
      w=(t(X) %*% u)/sum(u^2)
      w=w/sqrt(sum(w^2))
      t=X %*% w
      t=t/sum(t^2)
      c=t(Y) %*% t/ sum(t^2)
      newu= Y %*% c
      error=sum((u-newu)^2)
      u=newu
      if (show) print(c(i, iter, error))
    }
    T[,i]=t
    U[,i]=u
    W[,i]=w
    C[,i]=c
    p=t(X) %*% t/ sum(t^2)
    P[,i]=p
    q=t(Y) %*% u/ sum(u^2)
    Q[,i]=q
    X=X-t %*% t(p)
  }
  result=list(c0=c0, T=T, W=W, P=P, U=U, C=C, Q=Q)
  return(result)
}





