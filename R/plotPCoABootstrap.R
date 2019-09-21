plot.PCoABootstrap <- function(x, F1=1, F2=2, Move2Center=TRUE, BootstrapPlot="Ellipse", confidence=0.95, Colors=NULL, ...){
  if (class(x)!="PCoABootstrap") stop("Youmust provide an object of class PCoABootstrap")
  if (is.null(Colors)) Colors = matrix(1, n, 1)
  
  n=length(x$Coordinates)
  if (BootstrapPlot=="Ellipse")
    for (i in 1:n){
      ellipses=ConcEllipse(t(x$Coordinates[[i]][c(F1,F2),]), confidence=confidence)
      plot(ellipses, , col=Colors[i], ...)
    }
  else{
    for (i in 1:n){
      fraction=Fraction(t(x$Coordinates[[i]][c(F1,F2),]))
      plot(fraction, type=BootstrapPlot, col=Colors[i], ...)
    }
  }
  
}

