plot.fraction <- function(x, add=TRUE, center=FALSE, centerlabel="Center", initial=FALSE, type="ch",  ...){

  if (!add) plot(x$fraction[,1], x$fraction[,2], asp=1, cex=0)
  
  if (type == "ch"){
      hpts = chull(x$fraction)
      hpts <- c(hpts, hpts[1])
      lines(x$fraction[hpts, ], ...)
    }
  
  if (type == "st"){
    for (j in 1:(nrow(x$fraction))) {
      lines(rbind(x$fraction[j, ],x$center), ...)
  }}
  
  if (center) {
    points(x$center[1], x$center[2], pch=16, cex=1.2, ...)
    text(x$center[1], x$center[2], labels=centerlabel, cex=1.2, pos=4, ...)
  }
  if (initial) points(x$data[,1], x$data[,2], pch=16, cex=0.8, ...)
}