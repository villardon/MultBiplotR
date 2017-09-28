ColContributionPlot <- function(bip, A1=1, A2=2, Colors=NULL, Labs=NULL, MinQuality=0, ...){
  if (is.null(bip$Structure)) stop("Correlation circle is not available for this object")
  B=bip$ColContributions[,c(A1,A2)]/100
  Qualities=apply(B,1, sum)
  p=dim(B)[1]
  C=matrix(c(0, 0, 1, 0, 1, 1, -0, 1), 4, 2, byrow = TRUE)
  plot(C[,1], C[,2], cex=0, asp=1, xaxt = "n", yaxt = "n", xlab=paste("Axis", A1), ylab=paste("Axis", A2), bty="n", 
       main=paste(bip$Title, "- Contribution Plot"))
  if (is.null(Labs)) Labs=rownames(B)
  if (is.null(Colors)) Colors=rep("black", p)
  
  for (i in 1:p){
    if (Qualities[i]>MinQuality){
      arrows(0, 0, B[i,1], B[i,2], length = 0.1, angle = 20, col = Colors[i], ...)
      ang = atan(B[i,2]/B[i,1]) * 180/pi
      if (B[i,1]>0) pos=4
      else pos=2
      text(B[i,1], B[i,2], labels=Labs[i], col=Colors[i], srt=ang, offset = 0.2, pos = pos, ...) 
      }
    #text(c1, c2, label, cex = CexPoint, pos = markerpos, , srt = angle, col = Color, ...)
  }
  # rect(-1, -1, 1, 1)
  for (i in 1:10){
    QuarterCircle((i/10)^2, lty=3, color="red")
    text((i/10)^2, 0, labels=i/10, cex=0.7, pos=1, col="red")
    text(0, (i/10)^2, labels=i/10, cex=0.7, pos=2, col="red")
  }
  
}

QuarterCircle=function (radius = 1, origin = c(0, 0), color=1, ...) 
{
  t <- seq(0, pi/2, by = 0.01)
  a <- origin[1]
  b <- origin[2]
  r <- radius
  x <- a + r * cos(t)
  y <- b + r * sin(t)
  points(x, y, type = "l",col=color, ...)
}