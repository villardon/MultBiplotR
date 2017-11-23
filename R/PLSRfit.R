PLSRfit <- function(Y, X, S=2,  InitTransform=5, grouping=NULL,  centerY=TRUE, scaleY=TRUE, tolerance=0.000005, maxiter=100, show=FALSE){
  
  if (is.data.frame(X)) X=as.matrix(X)
  
  
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", 
                              "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization", "Ranks")
  if (is.numeric(InitTransform)) 
    InitTransform = ContinuousDataTransform[InitTransform]
  
  
  result=list()
  I1=dim(X)[1]
  J=dim(X)[2]
  
  if (is.numeric(Y)) Y= matrix(Y, ncol=1)
  
  I2=dim(Y)[1]
  K=dim(Y)[2]
  inames=rownames(X)
  ynames=colnames(Y)
  xnames=colnames(X)
  dimnames=paste("Component", 1:S)
  result$Method="PLSR1"
  result$X=X
  result$Y=Y
  result$centerY=centerY
  result$scaleY=scaleY
  result$Initial_Transformation=InitTransform
  if (!(I1==I2)) stop('The number of rows of both matrices must be the same')
  else I=I1
  
  
  rownames(Y)<-rownames(X)
  colnames(Y)="Response"
  
  
  if (centerY & !scaleY){
    Y=TransformIni(Y,transform=4)
  }
  if (centerY & scaleY){
    Y=TransformIni(Y,transform=5)
  }
  
  Data = InitialTransform(X, transform = InitTransform, grouping=grouping)
  X = Data$X
  if (InitTransform=="Within groups standardization") result$Deviations = Data$ColStdDevs
  
  result$ScaledX=X
  result$ScaledY=Y
  uini=svd(X)
  T=matrix(0, I1, S)
  rownames(T)=inames
  colnames(T)=dimnames
  U=matrix(0, I1, S)
  rownames(U)=inames
  colnames(U)=dimnames
  W=matrix(0, J, S)
  rownames(W)=xnames
  colnames(W)=dimnames
  C=matrix(0, K, S)
  rownames(C)=ynames
  colnames(C)=dimnames
  P=matrix(0, J, S)
  Q=matrix(0, K, S)
  
  for (i in 1:S){
    error=1
    iter=0
    u=matrix(uini$u[,i],I1,1)
    while ((error>tolerance) & (iter<maxiter)){
      iter=iter+1
      w=(t(X) %*% u)/sum(u^2)
      w=w/sqrt(sum(w^2))
      t=X %*% w
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
  
  rownames(P)=xnames
  colnames(P)=dimnames
  rownames(Q)=ynames
  colnames(Q)=dimnames
  
  result$XScores=T
  result$XWeights=W
  result$XLoadings=P
  result$YScores=U
  result$YWeights=C
  result$YLoadings=Q
  result$XStructure=cor(result$X,T)
  result$YStructure=cor(result$Y,U)
  result$YXStructure=cor(result$Y,T)
  class(result)="PLSR"
  return(result)
}
