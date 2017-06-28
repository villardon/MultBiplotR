summary.Binary.Logistic.Biplot <- function(object){
  
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
  
  RR=cbind(object$Deviances, object$Dfs, object$pvalues, object$Nagelkerke, object$CoxSnell, object$MacFaden, object$R2, object$PercentsCorrec*100)
  colnames(RR)=c("Deviance", "D.F", "P-val", "Nagelkerke", "Cox-Snell", "MacFaden", "R2", "% Correct")
  rownames(RR)=rownames(object$ColumnParameters)
  Total=c(object$DevianceTotal, object$TotalDf, object$p, object$TotNagelkerke, object$TotCoxSnell, object$TotMacFaden, object$TotR2, object$TotalPercent*100)
  RR=rbind(RR,Total)
  print(RR)
  
  print("-----------")
  print("Thresholds, Loadings and Communalities")
  LO=cbind(object$Tresholds, object$Loadings, object$Communalities)
  colnames(LO)=c("Thresholds", paste("Dim",1:dims,sep=""), "Communalities")
  rownames(LO)=rownames(object$ColumnParameters)
  print(LO)
}