anova.RidgeBinaryLogistic <- function(m1, m2){
  if (m1$Deviance > m2$Deviance){
    Difference= m1$Deviance - m2$Deviance
    DifDf=m2$df - m1$df}
  else{
    Difference= m2$Deviance - m1$Deviance
    DifDf=m1$df - m2$df}
  p=pchisq(Difference, DifDf, lower.tail = FALSE)
  res=rbind(c(m1$Deviance, m1$df, m1$p), c(m2$Deviance, m2$df, m2$p), c(Difference, DifDf, p))
  rownames(res)=c("Model 1", "Model 2", "Difference")
  colnames(res)=c("Deviance", "df", "p-value")
  return(res)
}
