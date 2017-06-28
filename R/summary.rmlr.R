summary.rmlr <- function(object, ...){
  if (object$Penalization==0)
    cat("Multinomial Logistic Regression\n")
  else{
    cat("Ridge Multinomial Logistic Regression\n")
    cat("Penalization : ", object$Penalization,"\n\n")}
  cat("\nSummary of the reponse\n")
  print(object$itab)
  cat("\n Model Fit :\n")
  fit=matrix(c(object$Deviance, object$Difference, object$df , object$p), 1,4)
  rownames(fit)="Final Model"
  colnames(fit)=c("    -2 log Lik  ", "Difference (Null)", "    df", "P-value")
  print(fit)
  cat("\n\n")
  Ps=rbind(object$CoxSnell, object$Nagelkerke, object$MacFaden)
  rownames(Ps)=c("CoxSnell", "Nagelkerke", "MacFaden")
  colnames(Ps)="Pseudo R-Squared"
  print(Ps)
  cat("\nOther fit indices (Information criteria)\n")
  cat("AIC: ", object$AIC,"\n")
  cat("BIC: ", object$BIC,"\n")
  
  nlevels=dim(object$beta)[1]
  lablevel=rownames(object$beta)
  
  cat("\nParameter estimates :\n")
  for (i in 1:nlevels){
    cat("Response level : ", lablevel[i],"\n")
    inf=cbind(object$beta[i,], object$stderr[i,],  (object$beta[i,]/object$stderr[i,])^2, 1 - pchisq((object$beta[i,]/object$stderr[i,])^2, df = 1),
              exp(object$beta[i,]), exp(object$beta[i,]-1.96*object$stderr[i,]), exp(object$beta[i,]+1.96*object$stderr[i,]))
    colnames(inf)=c("    Beta", "  Std. Error", "   Wald", "   p-value", "exp(Beta)", "CI : Lower", "CI : Upper")
    print(round(inf, digits=4))}
  
  cat("\nClassification :\n")
  cat("\nTable :\n")
  print(object$Table)
  cat("\nPercentages :\n")
  percent=round((diag(object$Table)/object$itab[,1])*100, digits=2)
  print(percent)
  cat("\nGlobal Percentage :", object$PercentCorrect,"\n")
}