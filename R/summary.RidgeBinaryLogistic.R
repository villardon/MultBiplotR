summary.RidgeBinaryLogistic <- function(object, ...) {
  n=dim(object$x)[1]
  p=dim(object$x)[2]
  cat("\nBINARY LOGISTIC REGRESSION - with Ridge  penalization\n\n")
  cat("Penalization :", object$Penalization, "\n")
  cat("\n Coefficients: \n")
  stderr=sqrt(diag(object$Covariances))
  Z=object$beta/stderr
  Betap=(1-pnorm(abs(Z)))*2
    Coef=cbind(round(object$beta,4), round(stderr, 4), round(Z, 3), round(Betap, 4), round(exp(object$beta),4))
    colnames(Coef)=c("Beta", "Std. Err.", "Z", "Pr(>|z|)", "Exp(B)")
  print(Coef)
  
  cat("\n Classification Table\n")
  print(object$Classification)
  
  cat("\n Classification Table (percentages)\n")
  print(round(prop.table(object$Classification, margin=1)*100, digits=2))
  cat("\n % Correct :",object$PercentCorrect*100)
  cat("\n\nNull deviance: ", object$NullDeviance, " on", n-1, "degrees of freedom")
  cat("\nResidual deviance: ", object$Deviance, " on", n-p, "degrees of freedom")
  cat("\nDifference: ", object$Dif, " on", p-1, "degrees of freedom (p=",object$p,")")
  cat("\nNagelkerke: ", object$Nagelkerke)
  cat("\nMacFaden: ", object$MacFaden)
  cat("\nCox-Snell: ", object$CoxSnell,"\n\n")
}


print.RidgeBinaryLogistic <- function(object, ...) {
  cat("\nBINARY LOGISTIC REGRESSION - with Ridge  penalization\n\n")
  cat("Penalization", object$Penalization, "\n")
  cat("\n Coefficients: \n")
  print(round(object$beta,4))
}