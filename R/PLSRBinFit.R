PLSRBinFit <- function(Y, X, S=2, tolerance=0.000005, maxiter=100, show=FALSE,
                       penalization=0.1, cte =TRUE, OptimMethod="CG"){
  
  I1=dim(X)[1]
  J=dim(X)[2]
  
  I2=dim(Y)[1]
  K=dim(Y)[2]
  I=I2
  
  
  TT=matrix(0, I, S)
  P=matrix(0, J, S)
  Q=matrix(0, K, S)
  
  # Set the constants
  q0=matrix(0,nrow=K, ncol=1)
  for (j in 1:K)
    q0[j]=RidgeBinaryLogistic(y=Y[,j], matrix(1,I,1), penalization = 0)$beta
  Q=q0
  U=matrix(1,I1,1)
  # We suppose thet the X variables are at least centered
  
  for (k in 1:S){
    u=matrix(X[,1],I1,1)
    U=cbind(U,u)
    parQ=rep(1,K)/sqrt(K)
    Q=cbind(Q,parQ)
    t=u
    error=1
    iter=0
    while ((error>tolerance) & (iter<maxiter)){
      iter=iter+1
      told=t
      w=(t(X) %*% t)/sum(t^2)
      w=w/sqrt(sum(w^2))
      t=X %*% w
      t=scale(t)
      u=t
      U[,k+1]=u
      # Update Q
      resbipQ <- optimr(parQ, fn=JLogBiplotRegBRec, gr=grLogBiplotRegBRec, method=OptimMethod, X=Y, A=U, B=Q,lambda=penalization)
      parQ=resbipQ$par
      Q[,k+1]=parQ
      # Update U
      resbipU <- optimr(u, fn=JLogBiplotRegARec, gr=grLogBiplotRegARec, method=OptimMethod, X=Y, A=U, B=Q, lambda=penalization)
      u=resbipU$par
      U[,k+1]=u
      t=u
      error=sum((told-t)^2)
      if (show) cat("\n",round(iter), round(J, 3), round(error,7))
    }
    t=X %*% w
    p=X%*%w
    TT[,k]=t
    P[,k]=w
    X = X - t %*% t(w)
  }
  
  U=U[,-1]
  Q=Q[,-1]

  Lin= cbind(1, TT) %*% t(cbind(q0,Q))
  Expected=exp(Lin)/(1+exp(Lin))
  C=solve(t(TT)%*%TT)%*%t(TT)%*%U
  B=P%*%C%*%t(Q)
  Pred=matrix(as.numeric(Expected>0.5), nrow=I)
  Right=(Y==Pred)
  PercentCorrect=sum(Right)/(I*K)
  PercentCorrectCols=apply(Right, 2, sum)/I
  
  result=list(TT=TT, P=P, U=U, Q=Q, q0=q0, B=B, Linterm=Lin, Expected=Expected, Predictions=Pred, PercentCorrect=PercentCorrect,
              PercentCorrectCols=PercentCorrectCols)
  return(result)
}
