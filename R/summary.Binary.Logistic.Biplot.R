summary.Binary.Logistic.Biplot <- function(object, Normal=TRUE, Latex=FALSE, Kable=FALSE, ...){
  
  print("BINARY LOGISTIC BIPLOT")
  print(paste("Type of Biplot : ", object$Type))
  print(paste("Initial Configuration : ", object$InitialConfig))
  print(paste("Method : ", object$Method))
  print(paste("Rotation : ", object$Rotation))
  print("-----------")
  dims=dim(object$RowCoordinates)[2]
  print("COLUMN PARAMETERS")
  print(object$ColumnParameters)
  
  print("-----------")
  print("COLUMNS FIT")
  
  RR=cbind(object$Deviances, object$Dfs, object$pvalues, object$Nagelkerke, object$CoxSnell, object$MacFaden, object$PercentsCorrec*100, object$Sensitivity*100, object$Specificity*100)
  colnames(RR)=c("Deviance", "D.F", "P-val", "Nagelkerke", "Cox-Snell", "MacFaden", "% Correct", "Sensitivity", "Specificity")
  rownames(RR)=rownames(object$ColumnParameters)
  Total=c(object$DevianceTotal, object$TotalDf, object$p, object$TotNagelkerke, object$TotCoxSnell, object$TotMacFaden, object$TotalPercent*100, object$TotalSensitivity*100, object$TotalSpecificity*100)
  RR=rbind(RR,Total)
  print(RR)
  print("------------------------")
  print("Thresholds, Loadings and Communalities")
  LO=cbind(object$Tresholds, object$Loadings, object$Communalities)
  colnames(LO)=c("Thresholds", paste("Dim",1:dims,sep=""), "Communalities")
  rownames(LO)=rownames(object$ColumnParameters)
  print(LO)
}