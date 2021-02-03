summary.RidgeBinaryLogistic <- function(x, ...) {
  n=dim(x$x)[1]
  p=dim(x$x)[2]
  cat("\nBINARY LOGISTIC REGRESSION - with Ridge  penalization\n\n")
  cat("Penalization :", x$Penalization, "\n")
  cat("\n Coefficients: \n")
  stderr=sqrt(diag(x$Covariances))
  Z=x$beta/stderr
  Betap=(1-pnorm(abs(Z)))*2
    Coef=cbind(round(x$beta,4), round(stderr, 4), round(Z, 3), round(Betap, 4), round(exp(x$beta),4))
    colnames(Coef)=c("Beta", "Std. Err.", "Z", "Pr(>|z|)", "Exp(B)")
  print(Coef)
  
  cat("\n Classification Table\n")
  print(x$Classification)
  
  cat("\n Classification Table (percentages)\n")
  print(round(prop.table(x$Classification, margin=1)*100, digits=2))
  cat("\n % Correct :",x$PercentCorrect*100)
  
  cat("\n\nNull deviance: ", x$NullDeviance, " on", n-1, "degrees of freedom")
  cat("\nResidual deviance: ", x$Deviance, " on", n-p, "degrees of freedom")
  cat("\nDifference: ", x$Dif, " on", p-1, "degrees of freedom (p=",x$p,")")
  cat("\nNagelkerke: ", x$Nagelkerke)
  cat("\nMacFaden: ", x$MacFaden)
  cat("\nCox-Snell: ", x$CoxSnell,"\n\n")
}




print.RidgeBinaryLogistic <- function(x, ...) {
  cat("\nBINARY LOGISTIC REGRESSION - with Ridge  penalization\n\n")
  cat("Penalization", x$Penalization, "\n")
  cat("\n Coefficients: \n")
  print(round(x$beta,4))
}