PLSR1BinFit <- function(Y, X, S=2, InitTransform=5, grouping=NULL, tolerance=0.000005,
                        maxiter=100, show=FALSE, penalization=0)
{
  
  if (is.data.frame(X)) X=as.matrix(X)
  if (!CheckBinaryVector(Y)) stop("The response must be binary (0 or 1)")
  
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
  
  if (is.numeric(Y)) {Y= as.matrix(Y, I1,1)
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
  result$Initial_Transformation=InitTransform

  
  if (!(I1==I2)) stop('The number of rows of both matrices must be the same')
  else I=I1
  
  Data = InitialTransform(X, transform = InitTransform, grouping=grouping)
  X = Data$X
  if (InitTransform=="Within groups standardization") result$Deviations = Data$ColStdDevs
  result$ScaledX=X
  result$ScaledY=Y
  result$tolerance=tolerance
  result$maxiter=maxiter
  result$penalization=penalization
  
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
  
  fit=RidgeBinaryLogistic(Y, T, penalization=penalization)
  
  rownames(P)=xnames
  colnames(P)=dimnames

  result$XScores=T
  result$XWeights=W
  result$XLoadings=P
  result$YWeights=fit$beta
  result$XStructure=cor(result$X,T)
  result$BinaryFit=fit
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
  Biplot$alpha=0
  Biplot$Initial_Transformation=plsr$Initial_Transformation
  Biplot$ncols=J
  Biplot$nrows=I
  Biplot$dim=S
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  if (plsr$Initial_Transformation == "Within groups standardization")  Biplot$Deviations = plsr$Deviations
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile)[2, ]
  Biplot$P75 = apply(X, 2, quantile)[4, ]
  
  a=plsr$XScores
  b=plsr$XLoadings
  sca = sum(a^2)
  scb = sum(b^2)
  sca = sca/I
  scb = scb/J
  scf = sqrt(sqrt(scb/sca))
  a = a * scf
  b = b/scf
  
  Biplot$RowCoordinates = a
  Biplot$ColCoordinates = b
  
  Cont=CalculateContributions(plsr$ScaledX,plsr$XScores,  plsr$XLoadings )
  Biplot$Inertia=Cont$Fit*100
  Biplot$RowContributions=Cont$RowContributions
  StResponse=cor(plsr$BinaryFit$linterm, plsr$XScores)
  rownames(StResponse)="Response"
  Biplot$Structure=Cont$Structure
  Biplot$ColContributions=Cont$ColContributions
  Biplot$SupStructure=StResponse
  Biplot$SupColContributions=StResponse^2
  
  class(Biplot)="ContinuousBiplot"
  Biplot=AddBinVars2Biplot(Biplot, plsr$Y, penalization=plsr$penalization, tolerance = plsr$tolerance, maxiter = plsr$maxiter)
  return(Biplot)
}