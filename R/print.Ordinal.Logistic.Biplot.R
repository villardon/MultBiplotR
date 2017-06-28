print.Ordinal.Logistic.Biplot <- function(x, ...){
  cat(x$Title)
  cat("\nType of Biplot", x$Type)
  cat("\nType of fit : ",x$Fit)
  cat("\nInitial Configuration : ",x$InitialConfiguration)
  cat("\nLog-Likelyhood : ", x$LogLikelihood)
  cat("\nAIC : ", x$AIC)
  cat("\nBIC : ", x$BIC)
  cat("\nNagelkerke pseudo R2 : ", x$Nagelkerke)
}