AddContVars2Biplot <- function(bip,  X, dims=NULL, Scaling = 5, Fit=NULL){
  n = nrow(X)
  p = ncol(X)
  if (is.null(dims)) dims=dim(bip$RowCoordinates)[2]
  # Setting the properties of data
  if (is.null(rownames(X))) 
    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "I")
  RowNames = rownames(X)
  if (is.null(colnames(X))) 
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "V")
  VarNames = colnames(X)
  
  Biplot = list()
  Biplot$Title = " Biplot"
  Biplot$Type = "External" 
  Biplot$Non_Scaled_Data = X
  Biplot$ncols=p
  Biplot$alpha=1
  Biplot$Dimension=dims
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile, 0.25)
  Biplot$P75 = apply(X, 2, quantile, 0.75)
  ContinuousDataTransform = c("Raw Data", "Substract the global mean", "Double centering", 
                              "Column centering", "Standardize columns", "Row centering", 
                              "Standardize rows", "Divide by the column means and center",
                              "Normalized residuals from independence", "Divide by the range",
                              "Within groups standardization")
  if (is.numeric(Scaling)) 
    Scaling = ContinuousDataTransform[Scaling]
  Biplot$Initial_Transformation = Scaling
  Data = InitialTransform(X, transform = Scaling)
  X = Data$X
  rownames(X) = RowNames
  colnames(X) = VarNames
  Biplot$Scaled_Data = X
  Biplot$Structure=cor(X,bip$RowCoordinates)
  
  DE=cbind(matrix(1,n,1),bip$RowCoordinates)
  b=t(solve(t(DE)%*%DE)%*%t(DE)%*%X)
  Biplot$b0=b[,1]
  Biplot$ColCoordinates=b[,2:(dims+1)]
  esp=DE %*% t(b)
  res=X-esp
  SCR=apply(res^2, 2, sum)
  SCT=apply(X^2, 2, sum)
  R2=(1-SCR/SCT)*100
  Biplot$R2=R2
  class(Biplot)="ContSupVarsBiplot"
  bip$ContSupVarsBiplot=Biplot
  return(bip)
}
