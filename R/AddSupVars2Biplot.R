# The following function adds information about supplementary variables to any already constructed biplot
# It may serve to add the variables to a  Principal Coordinates Analysis (or MDS)
# in order to build a biplot. Unless specified, numerical variables are added with linear regression, 
# factors with logistic biplots, frequencies with weighted averages and abundances with gaussian regression. 
# If the type of variable is not specified the program tries to guess.
# The type of Variable can be specified in a vector TypeVar and the type of fit in a vector TypeFit. 
# If fitting types are not supplied, default values for each type are used


AddSupVars2Biplot <- function(bip, X){
  
}

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
  Biplot$Means = apply(X, 2, mean)
  Biplot$Medians = apply(X, 2, median)
  Biplot$Deviations = apply(X, 2, sd)
  Biplot$Minima = apply(X, 2, min)
  Biplot$Maxima = apply(X, 2, max)
  Biplot$P25 = apply(X, 2, quantile)[2, ]
  Biplot$P75 = apply(X, 2, quantile)[4, ]
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


plot.ContSupVarsBiplot <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                   ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, PchVar=1, ColorVar=1, ...){
  A =x$RowCoordinates[, c(F1, F2)]
  B=x$ColCoordinates[,c(F1,F2)]
  
  n = dim(A)[1]
  b0=x$b0
  VarLabels=rownames(B)
  if (mode=="s")
    Scales = GetBiplotScales(x, TypeScale = TypeScale, ValuesScale = ValuesScale)
  
  p=dim(x$ColCoordinates)[1]
  if (is.numeric(B)) B=matrix(B, nrow=p)
  for (j in 1:p) 
    VarBiplot(B[j, 1], B[j, 2], b0=b0[j], xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, label = VarLabels[j], mode = mode, 
              ticks = Scales$Ticks[[j]], ticklabels = Scales$Labels[[j]], ts = TypeScale, PchPoint = PchVar[j], tl=0.01, ...)
  
}



AddBinVars2Biplot <- function(bip, Y, IncludeConst=TRUE, penalization=0.2, freq=NULL, tolerance = 1e-05, maxiter = 100) {
  
  n=dim(Y)[1]
  p=dim(Y)[2]
  dimens=dim(bip$RowCoordinates)[2]
  x=bip$RowCoordinates
  Pco=list()
  Pco$ColumnParameters=matrix(0,p,dimens+1)
  Res=list()
  Res$Deviances=matrix(0,p,1)
  Res$Dfs=matrix(0,p,1)
  Res$pvalues=matrix(0,p,1)
  Res$Bonferroni=matrix(0,p,1)
  Res$Nagelkerke=matrix(0,p,1)
  Res$R2=matrix(0,p,1)
  Res$PercentsCorrec=matrix(0,p,1)
  Pco$DevianceTotal=0
  Pco$p=1
  Pco$TotalPercent=0
  
  for (i in 1:p){
    y=Y[,i]
    fit=RidgeBinaryLogistic(y,x,tolerance = tolerance, maxiter = maxiter, penalization=penalization, cte=IncludeConst)
    if (IncludeConst)
      Pco$ColumnParameters[i,]=fit$beta
    else
      Pco$ColumnParameters[i,]=c(0,fit$beta)
      
    Res$Deviances[i]=fit$Dif
    Res$Dfs[i]=fit$df
    Res$pvalues[i]=fit$p
    Res$Bonferroni[i]=(fit$p * p)* ((fit$p * p)<=1) + (((fit$p * p)>1))
    Res$Nagelkerke[i]=fit$Nagelkerke
    Res$R2[i]=fit$R2
    Res$PercentsCorrec[i]=fit$PercentCorrect
    Pco$TotalPercent=Pco$TotalPercent+sum(y==fit$Prediction)
  }
  if (IncludeConst)
  rownames(Pco$ColumnParameters)=colnames(Pco$Data)
  Pco$TotalPercent=Pco$TotalPercent/(n*p)
  Pco$DevianceTotal=sum(Res$Deviances)
  Pco$TotalDf=sum(Res$Dfs)
  Pco$p=1-pchisq(Pco$DevianceTotal, df = Pco$TotalDf)
  Res=as.data.frame(Res)
  rownames(Res)=colnames(Y)
  Pco$VarInfo=Res
  class(Pco)="BinSupVarsBiplot"
  bip$BinSupVarsBiplot=Pco
  return(bip)
}

plot.Supplementary.Variables <- function(bip, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                         ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, ...){
  if (!is.null(bip$ContSupVarsBiplot))
    plot(bip$ContSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode, TypeScale=TypeScale) 
  
  if (!is.null(bip$BinSupVarsBiplot))
    plot(bip$BinSupVarsBiplot, F1=F1, F2=F2, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, mode=mode) 
  
}


plot.BinSupVarsBiplot <- function(x, F1=1, F2=2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, TypeScale = "Complete", 
                                  ValuesScale = "Original", mode="s", dp = 0, PredPoints=0, ...){
  p=dim(x$ColumnParameters)[1]
  ColLabels=rownames(x$VarInfo)
  
  for (i in 1:p)
    PlotBinaryVar(b0=x$ColumnParameters[i,1], bi1=x$ColumnParameters[i,F1+1], bi2=x$ColumnParameters[i,F2+1], xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                  mode=mode, label=ColLabels[i], Color="blue")
}
