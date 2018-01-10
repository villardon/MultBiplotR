plot.PLSR <- function(plsr, ParameterBoxPlot=FALSE, ParameterCI=TRUE, Correlations=FALSE, ...){
  
  I=length(plsr$RegParameters[,1])
  VarLabels= rownames(plsr$RegParameters)
  if (plsr$Crossvalidation){
    if (ParameterBoxPlot){
      boxplot(t(plsr$CrossParameters), main="Box plot for the Jackkinfe pseudo-samples")
      abline(h=0, col="red")}
    
    if (ParameterCI){
      signif=plsr$RegParameters[,4]>0.05
      plot(1:I, plsr$RegParameters[,1], ylim=range(c(plsr$RegParameters[,5], 
                                                     plsr$RegParameters[,6])), pch=19, xlab="Variables", ylab="Mean +/- CI", 
           main="Scatter plot with confidence intervals error", col=signif+1)
      arrows(1:I, plsr$RegParameters[,5], 1:I, plsr$RegParameters[,6], length=0.05, angle=90, code=3, col=signif+1)
      text(1:I, plsr$RegParameters[,6], labels = VarLabels, pos=3, col=signif+1)
      abline(h=0, col="red", lty=3)
    }
    
  }
}
