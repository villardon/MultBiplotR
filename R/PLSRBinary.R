PLSR1BinFit <- function(Y, X, S=2, center=TRUE, scale=TRUE, tolerance=0.000005,
                        maxiter=100, show=FALSE, penalization=0)
{
  
  if (!CheckBinaryVector(Y)) stop("The response must be binary (0 or 1)")
  result=list()
  I1=dim(X)[1]
  J=dim(X)[2]
  
  if (is.numeric(Y)) {Y= as.matrix(Y)
                      rownames(Y)<-rownames(X)
                      colnames(Y)="Response"}
  I2=dim(Y)[1]
  K=dim(Y)[2]
  inames=rownames(X)
  ynames=colnames(Y)
  xnames=colnames(X)
  dimnames=paste("Component", 1:S)
  
  
  result$Method="PLSR1 Binary"
  result$X=X
  result$Y=Y
  result$center=center
  result$scale=scale
  
  if (!(I1==I2)) stop('The number of rows of both matrices must be the same')
  else I=I1
  
  if (center & !scale){
    X=TransformIni(X,transform=4)
    result$Initial_Transformation=4
  }
  if (center & scale){
    X=TransformIni(X,transform=5)
    result$Initial_Transformation=5
  }
  
  result$ScaledX=X
  result$ScaledY=Y

  T=matrix(0, I, S)
  rownames(T)=inames
  colnames(T)=dimnames
  W=matrix(0, J, S)
  rownames(W)=xnames
  colnames(W)=dimnames
  C=matrix(0, K, S)
  rownames(C)=ynames
  colnames(C)=dimnames
  P=matrix(0, J, S)
  Q=matrix(0, K, S)
  freq=matrix(1,I,1)
  w=matrix(0,J,1)
  
  for (j in 1:J){
    x=as.matrix(X[,j])
    colnames(x)=xnames[j]
    fit=RidgeBinaryLogistic(Y, x)
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
      colnames(x)=c(dimnames[1:(i-1)] ,xnames[j])
      fit=RidgeBinaryLogistic(Y, x)
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
  
  fit=RidgeBinaryLogistic(Y, T)
  
  rownames(P)=xnames
  colnames(P)=dimnames
  
  result$XScores=T
  result$XWeights=W
  result$XLoadings=P
  result$YWeights=fit$beta
  result$XStructure=cor(result$X,T)
  class(result)="PLSR1Bin"
  return(result)
}


biplot.PLSR1BIN <- function(plsr, ... ){
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