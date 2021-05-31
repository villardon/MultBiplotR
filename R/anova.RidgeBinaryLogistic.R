anova.RidgeBinaryLogistic <- function(object, object2, ...){
  if (object$Deviance > object2$Deviance){
    Difference= object$Deviance - object2$Deviance
    DifDf=object2$df - object$df}
  else{
    Difference= object2$Deviance - object$Deviance
    DifDf=object$df - object2$df}
  p=pchisq(Difference, DifDf, lower.tail = FALSE)
  res=rbind(c(object$Deviance, object$df, object$p), c(object2$Deviance, object2$df, object2$p), c(Difference, DifDf, p))
  rownames(res)=c("Model 1", "Model 2", "Difference")
  colnames(res)=c("Deviance", "df", "p-value")
  return(res)
}
