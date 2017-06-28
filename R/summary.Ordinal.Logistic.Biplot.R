summary.Ordinal.Logistic.Biplot <- function(object, ...){
  cat(object$Title)
  cat("\nType of Biplot : ", object$Type)
  cat("\nType of fit : ",object$Fit)
  cat("\nInitial (External) Configuration : ",object$InitialConfiguration)
  cat("\nPenalization for the logistic regression : ",object$Penalization)
  cat("\n\n\nCOLUMN PARAMETERS")
  cat("\nThresholds\n")
  print(object$ColumnParameters$thresholds)
  cat("\n\nSlopes (Discriminations) \n")
  print(object$ColumnParameters$coefficients)
  n=dim(object$RowCoordinates)[1]
  p=dim(object$Communalities)[1]
  percentcorrect=sum(object$ColumnParameters$fit[,5]*n)/(n*p)
  global=matrix(c(object$LogLikelihood, object$Deviance, object$df, object$pval, percentcorrect, object$CoxSnell, object$MacFaden, object$Nagelkerke ), 1,8)
  rownames(global)="Global"
  colnames(global)=colnames(object$ColumnParameters$fit[,1:8])
  cat("\n\nMeasures of fit (Global and for each separate variable) \n")
  print(rbind(object$ColumnParameters$fit[,1:8],global))
  cat("\n\nFactor Structure (Loadings and Communalities) \n")
  print(cbind(object$loadings, object$Communalities))
  
  cat("\nAIC : ", object$AIC)
  cat("\nBIC : ", object$BIC)
}
