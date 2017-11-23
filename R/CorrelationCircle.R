CorrelationCircle <- function(bip, A1=1, A2=2, Colors=NULL, Labs=NULL, ...){
  if (is.null(bip$Structure)) stop("Correlation circle is not available for this object")
  B=bip$Structure[,c(A1,A2)]
  p=dim(B)[1]
  C=matrix(c(-1, -1, 1, -1, 1, 1, -1, 1), 4, 2, byrow = TRUE)
  plot(C[,1], C[,2], cex=0, asp=1, xaxt = "n", yaxt = "n", xlab=paste("Axis", A1), ylab=paste("Axis", A2), bty="n", 
       main=paste(bip$Title, "- Correlation Circle"))
  if (is.null(Labs)) Labs=rownames(B)
  if (is.null(Colors)) Colors=rep("black", p)
  
  if (!is.null(bip$SupStructure)){
    nvs=dim(bip$SupStructure)[1]
    p=p+nvs
    B=rbind(B, bip$SupStructure )
    Colors=c(Colors, rep("green", nvs))
    Labs=c(Labs, rownames(bip$SupStructure))
  }
  
  for (i in 1:p){
  arrows(0, 0, B[i,1], B[i,2], length = 0.1, angle = 20, col = Colors[i], ...)
    ang = atan(B[i,2]/B[i,1]) * 180/pi
    if (B[i,1]>0) pos=4
    else pos=2
    text(B[i,1], B[i,2], labels=Labs[i], col=Colors[i], srt=ang, offset = 0.2, pos = pos, ...) 
    #text(c1, c2, label, cex = CexPoint, pos = markerpos, , srt = angle, col = Color, ...)
    }
  
  for (i in 1:10){
    Circle(i/10, lty=3, col="red")
    text(i/10, 0, labels=i/10, cex=0.5, col="red")
    text(0, i/10, labels=i/10, cex=0.5, col="red")
    text(-1*i/10, 0, labels=i/10, cex=0.5, col="red")
    text(0, -1*i/10, labels=i/10, cex=0.5, col="red")
  }
  
}
