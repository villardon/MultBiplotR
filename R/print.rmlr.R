print.rmlr <- function(x, ...){
  if (x$Penalization==0)
    cat("Ridge Multinomial Logistic Regression\n")
  else{
  cat("Ridge Multinomial Logistic Regression\n")
  cat("Penalization : ", x$Penalization,"\n\n")}
  cat("Coefficients :\n")
  print(x$beta)
  cat("\n")
  cat("Residual Deviance: ", x$Difference,"\n")
  cat("AIC: ", x$AIC,"\n")
  cat("BIC: ", x$BIC,"\n")
}
