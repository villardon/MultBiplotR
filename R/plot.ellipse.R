plot.ellipse <- function(x, add=TRUE, labeled= FALSE , center=FALSE, centerlabel="Center", initial=FALSE,  ...){
  if (add)
    points(x$ellipse, type="l", ...)
  else
    plot(x$ellipse, type="l", ...)
  
  if (labeled) text(x$ellipse[1,1], x$ellipse[1,2], labels=x$confidence, pos=4, ... )
  if (center) {
    points(x$center[1], x$center[2], pch=16, cex=1.2, ...)
    text(x$center[1], x$center[2], labels=centerlabel, cex=1.2, pos=4, ...)
  }
  
  if (initial) points(x$data[,1], x$data[,2], pch=16, cex=0.8, ...)

}

