PLSR1Bin <- function(Y, X, S=2, InitTransform=5, grouping=NULL, tolerance=0.000005,
                     maxiter=100, show=FALSE, penalization=0, cte =TRUE, Algorithm=1,
                     OptimMethod="CG"){
  
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
  K=1
  inames=rownames(X)
  ynames=colnames(Y)
  xnames=colnames(X)
  dimnames=paste("Comp.", 1:S)
  
  
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
  result$IncludeConst=cte
  
  myfit=PLSR1BinFit(Y, X, S=S, tolerance=tolerance, maxiter=maxiter, show=show, penalization=penalization, cte =cte, Algorithm=Algorithm)
  
  rownames(myfit$T)=inames
  colnames(myfit$T)=dimnames
  rownames(myfit$W)=xnames
  colnames(myfit$W)=dimnames
  C=matrix(0, K, S)
  rownames(myfit$C)=ynames
  colnames(myfit$C)=dimnames
  rownames(myfit$P)=xnames
  colnames(myfit$P)=dimnames
  
  result$XScores=myfit$T
  result$XWeights=myfit$W
  result$XLoadings=myfit$P
  result$YWeights=myfit$beta
  result$XStructure=cor(result$X,myfit$T)
  result$BinaryFit=myfit$fit
  
  class(result)="PLSR1Bin"
  return(result)
}