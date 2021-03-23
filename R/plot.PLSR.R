plot.PLSR <- function(x, ParameterBoxPlot=FALSE, ParameterCI=TRUE, Correlations=FALSE, ...){
  
  I=length(x$RegParameters[,1])
  VarLabels= rownames(x$RegParameters)
  if (x$Validation=="Cross"){
    if (ParameterBoxPlot){
      boxplot(t(x$CrossParameters), main="Box plot for the Jackkinfe pseudo-samples")
      abline(h=0, col="red")}
    
    if (ParameterCI){
      signif=x$RegParameters[,4]>0.05
      plot(1:I, x$RegParameters[,1], ylim=range(c(x$RegParameters[,5], 
                                                     x$RegParameters[,6])), pch=19, xlab="Variables", ylab="Mean +/- CI", 
           main="Scatter plot with confidence intervals error", col=signif+1)
      arrows(1:I, x$RegParameters[,5], 1:I, x$RegParameters[,6], length=0.05, angle=90, code=3, col=signif+1)
      text(1:I, x$RegParameters[,6], labels = VarLabels, pos=3, col=signif+1)
      abline(h=0, col="red", lty=3)
    }
    
  }
}
