PLSR1BinFit <- function(Y, X, S=2, tolerance=0.000005, maxiter=100, show=FALSE, 
                        penalization=0, cte =TRUE, Algorithm=2, OptimMethod="CG"){
  
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
  
  freq=matrix(1,I,1)
  w=matrix(0,J,1)
  # El algoritmo 1 se basa en Bastien et al.
  
  if (Algorithm==1){ 
    t0=matrix(1, nrow=I, ncol=1)
    fit=RidgeBinaryLogistic(Y, t0, penalization=penalization, cte=FALSE)
    c0=fit$beta[1]
    
    for (j in 1:J){
      x=as.matrix(X[,j])
      fit=RidgeBinaryLogistic(Y, x, cte=TRUE)
      w[j]=fit$beta[2]
    }
    w=w/sqrt(sum(w^2))
    t=X %*% w
    T[,1]=t
    W[,1]=w
    p=t(X) %*% t/ sum(t^2)
    P[,1]=p
    X1=X-t %*% t(p)
    
    for (i in 2:S){
      for (j in 1:J){
        x=as.matrix(cbind(T[,1:(i-1)],X[,j]))
        fit=RidgeBinaryLogistic(Y, x, cte=TRUE)
        w[j]=fit$beta[(i+1)]
      }
      w=w/sqrt(sum(w^2))
      t=X1 %*% w
      T[,i]=t
      W[,i]=w
      p=t(X) %*% t/ sum(t^2)
      P[,i]=p
      X1=X1-t %*% t(p)
    }
    fit=RidgeBinaryLogistic(Y, T, penalization=penalization, cte=cte)
  }
  
  # El algoritmo 2 
  if (Algorithm==2){
      t0=matrix(1, nrow=I, ncol=1)
      fit=RidgeBinaryLogistic(Y, t0, penalization=penalization, cte=FALSE)
      q0=fit$beta[1]
      
      
    # We have to take the constant into account
    for (i in 1:S){
      error=1
      iter=0
      u= Y
      
      while ((error>tolerance) & (iter<maxiter)){
        iter=iter+1
        w=(t(X) %*% u)/sum(u^2)
        w=w/sqrt(sum(w^2))
        t=X %*% w
        TT=cbind(t0, t)
        
        q=fit$beta[2]
        Q=c(q0, q)
        CostBinPLS1B(q , Y=Y, U=TT, Q=Q)
        resbipB <- optimr(q, fn=CostBinPLS1B, gr=grBinPLS1B, method=OptimMethod, Y=Y, U=TT, Q=Q)
        q=resbipB$par
        Q=c(c0, q)
        
        grBinPLS1A(u, Y=Y, U=TT, Q=Q)
        parA <- optimr(u, fn=CostBinPLS1A, gr=grBinPLS1A, method=OptimMethod, Y=Y, U=TT, Q=Q)
        newu=matrix(parA$par, ncol=1)
        
        error=sum((u-newu)^2)
        u=newu
        if (show) print(c(i, iter, error))
      }
      T[,i]=t
      U[,i]=u
      W[,i]=w
      C[,i]=c
      Q[,i]=q
      p=t(X) %*% t/ sum(t^2)
      P[,i]=p
      X=X-t %*% t(p)
    }
    fit=RidgeBinaryLogistic(Y, T, penalization=penalization, cte=TRUE)
  }
  result=list(T=T, W=W, P=P, U=U, C=C, Q=Q, fit=fit)
  return(result)
}