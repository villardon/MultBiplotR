PLSRfit <- function(Y, X, S=2, center=TRUE, scale=TRUE, tolerance=0.000005, maxiter=100, show=FALSE)
{
  result=list()
  I1=dim(X)[1]
  J=dim(X)[2]
  I2=dim(Y)[1]
  K=dim(Y)[2]
  inames=rownames(X)
  ynames=colnames(Y)
  xnames=colnames(X)
  dimnames=paste("Component", 1:S)
  result$Method="PLSR1"
  result$X=X
  result$Y=Y
  result$center=center
  result$scale=scale
  
  if (!(I1==I2)) stop('The number of rows of both matrices must be the same')
  else I=I1
  
  if (center & !scale){
    X=TransformIni(X,transform=4)
    Y=TransformIni(Y,transform=4)
    result$Initial_Transformation=4
  }
  if (center & scale){
    X=TransformIni(X,transform=5)
    Y=TransformIni(Y,transform=5)
    result$Initial_Transformation=5
  }
  
  result$ScaledX=X
  result$ScaledY=Y
  uini=svd(X)
  T=matrix(0, I, S)
  rownames(T)=inames
  colnames(T)=dimnames
  U=matrix(0, I, S)
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
    u=matrix(uini$u[,i],I,1)
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


biplot.PLSR <- function(plsr, ... ){
  X=plsr$X
  Y=plsr$Y
  
  I=dim(X)[1]
  J=dim(X)[2]
  K=dim(Y)[2]
  S=dim(plsr$XScores)[2]
  
  Biplot = list()
  Biplot$Title = " PLSR - Biplot"
  Biplot$Type = "PLSR" 
  Biplot$Initial_Transformation=plsr$Initial_Transformation
  Biplot$ncols=J
  Biplot$nrows=I
  Biplot$dim=S
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile)[2, ]
  Biplot$P75 = apply(X, 2, quantile)[4, ]
  Biplot$RowCoordinates = plsr$XScores
  Biplot$ColCoordinates = plsr$XLoadings
  Cont=CalculateContributions(plsr$ScaledX,plsr$XScores,  plsr$XLoadings )
  Biplot$RowContributions=Cont$RowContributions
  Biplot$ColContributions=Cont$ColContributions
  Biplot$Structure=Cont$Structure
  class(Biplot)="ContinuousBiplot"
  YBiplot=list()
  YBiplot$Title = " PLSR - Biplot (Y)"
  YBiplot$Type = "PLSR" 
  YBiplot$Initial_Transformation=plsr$Initial_Transformation
  YBiplot$ncols=K
  YBiplot$Means = apply(Y, 2, mean)
  YBiplot$Medians = apply(Y, 2, median)
  YBiplot$Deviations = apply(Y, 2, sd)
  YBiplot$Minima = apply(Y, 2, min)
  YBiplot$Maxima = apply(Y, 2, max)
  YBiplot$P25 = apply(Y, 2, quantile)[2, ]
  YBiplot$P75 = apply(Y, 2, quantile)[4, ]
  YBiplot$b0 = rep(0,K)
  
  YBiplot$ColCoordinates = plsr$YWeights
  Cont=CalculateContributions(plsr$ScaledY,plsr$XScores,  plsr$YWeights)
  YBiplot$ColContributions=Cont$ColContributions
  YBiplot$Structure=Cont$Structure
  class(YBiplot)="ContSupVarsBiplot"
  Biplot$ContSupVarsBiplot = YBiplot
  
  return(Biplot)
}