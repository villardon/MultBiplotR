plot.Procrustes <- function(x, F1=1, F2=2, ...){
  dev.new()
  A=cbind(x$X,x$Yrot)
  plot(A[,F1],A[,F2], cex=0, asp=1, ...)
  points(x$X[,F1],x$X[,F2], ...)
  text(x$X[,F1],x$X[,F2],rownames(x$X), ...)
  plines(x$X,x$Yrot)
  
}