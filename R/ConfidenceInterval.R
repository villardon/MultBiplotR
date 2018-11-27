ConfidenceInterval<- function(x, Desv=NULL, df=NULL, Confidence=0.95){
  if (!is.numeric(x)) stop("You must use Numeric data to calculate a confidence interval")
  n=length(x)
  if (is.null(Desv)) Desv = sd(x)
  if (is.null(df)) df=length(x)-1
  prob=1-(1-Confidence)/2
  error <- qt(prob,df=df)*Desv/sqrt(n)
  left <- mean(x)-error
  right <- mean(x)+error
  return(c(left, right))
}
