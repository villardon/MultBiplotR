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
