PLSRBin <- function(Y, X, S=2, InitTransform=5, grouping=NULL, tolerance=0.00005,
                     maxiter=100, show=FALSE, penalization=0.1, cte =TRUE,
                     OptimMethod="CG", Multiple=FALSE){

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

  I2=dim(Y)[1]
  K=1
  inames=rownames(X)
  ynames=colnames(Y)
  xnames=colnames(X)
  dimnames=paste("Comp.", 1:S)


  result$Method="PLSR for binary responses"
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

  myfit=PLSRBinFit(Y=Y, X=X, S=S, tolerance=tolerance, maxiter=maxiter, show=show, penalization=penalization, cte =cte)

  rownames(myfit$TT)=inames
  colnames(myfit$TT)=dimnames
  C=matrix(0, K, S)
  rownames(myfit$B)=xnames
  colnames(myfit$B)=ynames
  rownames(myfit$P)=xnames
  colnames(myfit$P)=dimnames

  result$XScores=myfit$TT
  result$XLoadings=myfit$P
  result$YScores=myfit$U
  result$YLoadings=myfit$Q
  rownames(result$YLoadings)=ynames
  colnames(result$YLoadings)=paste("Dim", 1:S)
  result$YWeights=myfit$B
  result$XStructure=cor(result$X,myfit$TT)
  result$BinaryFits=myfit$fit
  result$Intercepts=myfit$q0
  result$LinTerm=myfit$Linterm
  result$Expected=myfit$Expected
  result$Predictions=myfit$Predictions
  rownames(result$Predictions)=inames
  colnames(result$Predictions)=ynames
  result$PercentCorrect=myfit$PercentCorrect
  result$PercentCorrectCols=myfit$PercentCorrectCols
  
  maxima=rep(0,I2)
  for (i in 1:I2)
    maxima[i]=which(result$Expected[i,]==max(result$Expected[i,]))
  result$Maxima=maxima

  class(result)="PLSRBin"
  return(result)
}
