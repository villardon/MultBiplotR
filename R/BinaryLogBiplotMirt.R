BinaryLogBiplotMirt <- function(x, dimens = 2, tolerance = 1e-04, maxiter = 30, penalization=0.2,  Rotation = "varimax", ...){
  # joint algorithm for logistic biplots
  n <- nrow(x)
  p <- ncol(x)

  print("Calculating  Row Coordinates - MIRT")
  
  mod2 <- mirt(x,dimens, ...)
  a <- fscores(mod2,method = "EAP", rotate=Rotation,  full.scores = TRUE)
  
  # # Centering the coordinates
  # med=apply(a,2,mean)
  # a = a - matrix(1,n,1) %*% matrix(med, 1, dimens)
  
  print("  ")
  print("Calculating Column Coordinates - Logistic Regression")
  
  Res=list()
  Res$Biplot="Binary Logistic"
  Res$RowCoordinates=a
  Res$ColumnParameters=matrix(0,p,dimens+1)
  Res$NullDeviances=matrix(0,p,1)
  Res$ModelDeviances=matrix(0,p,1)
  Res$Deviances=matrix(0,p,1)
  Res$Dfs=matrix(0,p,1)
  Res$pvalues=matrix(0,p,1)
  Res$Bonferroni=matrix(0,p,1)
  Res$Nagelkerke=matrix(0,p,1)
  Res$R2=matrix(0,p,1)
  Res$PercentsCorrec=matrix(0,p,1)
  Res$DevianceTotal=0
  Res$p=1
  Res$TotalPercent=0
  Res$SSRes=matrix(0,p,1)
  Res$SSTot=matrix(0,p,1)
  
  for (i in 1:p){
    cat(paste(" ",i))
    y=x[,i]
    fit=RidgeBinaryLogistic(y,a,tolerance = tolerance, maxiter = maxiter, penalization=penalization, cte=TRUE)
    Res$ColumnParameters[i,]=t(fit$beta)
    Res$ModelDeviances[i]=fit$Deviance
    Res$NullDeviances[i]=fit$NullDeviance
    Res$Deviances[i]=fit$Dif
    Res$Dfs[i]=fit$df
    Res$pvalues[i]=fit$p
    Res$Bonferroni[i]=(fit$p * p)* ((fit$p * p)<=1) + (((fit$p * p)>1))
    Res$R2[i]=fit$R2
    Res$CoxSnell[i]=fit$CoxSnell
    Res$Nagelkerke[i]=fit$Nagelkerke
    Res$MacFaden[i]=fit$MacFaden
    Res$PercentsCorrec[i]=fit$PercentCorrect
    Res$TotalPercent=Res$TotalPercent+sum(y==fit$Prediction)
    Res$SSRes[i]=fit$SSRes
    Res$SSTot[i]=fit$SSTot
  }
  
  d = sqrt(rowSums(cbind(1,Res$ColumnParameters[, 2:(dimens + 1)])^2))
  Res$Loadings = solve(diag(d)) %*% Res$ColumnParameters[, 2:(dimens + 1)]
  Res$Tresholds = Res$ColumnParameters[, 1]/d
  Res$Communalities = rowSums(Res$Loadings^2)
  
  rownames(Res$ColumnParameters)=colnames(x)
  colnames(Res$ColumnParameters)=c("Const.",paste("Dim",1:dimens, sep=""))
  
  rownames(Res$ColumnParameters)=colnames(x)
  
  Res$TotalPercent=Res$TotalPercent/(n*p)
  Res$ModelDevianceTotal=sum(Res$ModelDeviances)
  Res$NullDevianceTotal=sum(Res$NullDeviances)
  Res$DevianceTotal=sum(Res$Deviances)
  
  Res$TotalSSRes=sum(Res$SSRes)
  Res$TotalSSTot=sum(Res$SSTot)
  
  nn=length(x)
  Res$TotCoxSnell=1-exp(-1*Res$DevianceTotal/nn)
  Res$TotNagelkerke=Res$TotCoxSnell/(1-exp((Res$NullDevianceTotal/(-2)))^(2/nn))
  Res$TotMacFaden=1-(Res$ModelDevianceTotal/Res$NullDevianceTotal)
  Res$TotR2=1-(Res$TotalSSRes=sum(Res$SSRes)/Res$TotalSSTot)
  Res$TotalDf=sum(Res$Dfs)
  Res$p=1-pchisq(Res$DevianceTotal, df = Res$TotalDf)
  Res$ClusterType="us"
  Res$Clusters = as.factor(matrix(1,nrow(Res$RowCoordinates), 1))
  Res$ClusterColors="blue"
  Res$ClusterNames="ClusterTotal"
  class(Res) = "Binary.Logistic.Biplot"
  
  return(Res)
}

