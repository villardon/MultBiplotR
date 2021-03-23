plot.LASSOReg <- function(x, ...){
  p=dim(x$Beta)[1]
  ymin=min(x$Beta)
  ymax=max(x$Beta)
  op <- par(mfrow=c(2,1))
  plot(x$Lambdas,x$R2, type="l", xlab="lambda", ylab="R2", main="Goodness of fit")
  
  plot(x$Lambdas,x$Beta[1,], type="l", col=1, ylim=c(ymin,ymax), xlab="lambda", ylab="Betas", main="Coefficients")
  for (j in 2:p){
    points(x$Lambdas,x$Beta[j,], type="l", col=j)
  }
  par(op)
}