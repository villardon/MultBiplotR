ColContributionPlot <- function(bip, A1=1, A2=2, Colors=NULL, Labs=NULL, MinQuality=0, CorrelationScale=FALSE, ContributionScale=TRUE, AddSigns2Labs=TRUE,  ...){
  B=bip$Structure[,c(A1,A2)]
  Qualities=apply(B^2,1, sum)
  p=dim(B)[1]
  C=matrix(c(0, 0, 1, 0, 1, 1, -0, 1), 4, 2, byrow = TRUE)
  # plot(C[,1], C[,2], cex=0, asp=1, xaxt = "n", yaxt = "n", xlab=paste("Axis", A1), ylab=paste("Axis", A2), bty="n", main=paste(bip$Title, "- Contribution Plot"))
  
  plot(C[,1], C[,2], cex=0, asp=1, xlab=paste("Axis", A1), ylab=paste("Axis", A2), main=paste(bip$Title, "- Contribution Plot"), xlim=c(0, 1), ylim=c(0, 1), axes=F)
  #grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  #axis(1, pos=0)
  #axis(2, pos=0)
  if (is.null(Labs)) Labs=rownames(B)
  if (is.null(Colors)) Colors=rep("black", p)
  
  if (!is.null(bip$SupStructure)){
    nvs=dim(bip$SupStructure)[1]
    p=p+nvs
    B=rbind(B, bip$SupStructure )
    Colors=c(Colors, rep("green", nvs))
    Labs=c(Labs, rownames(bip$SupStructure))
    Qualities=c(Qualities, apply(bip$SupStructure^2,1, sum))
  }
  
  if (AddSigns2Labs){
    Labs= paste(Labs,"(",sign(B[,1]),",",sign(B[,2]),")", sep="")
    
    signs=sign(B[,1])*sign(B[,2])
    for (i in 1:p){
      if ((sign(B[i,1]) < 0 ) & (sign(B[i,2]) < 0)) Colors[i]="blue"
      if ((sign(B[i,1]) < 0 ) & (sign(B[i,2]) > 0)) Colors[i]="pink"
      if ((sign(B[i,1]) > 0 ) & (sign(B[i,2]) < 0)) Colors[i]="orange"
    }

    legend(x="topright", c("Q1 ++", "Q3 --", "Q2 -+", "Q4 +-") , col = c("black", "blue", "pink", "orange"), text.col = "black", pch = 16, cex=2)
    
  }
  
  B=abs(B)
  for (i in 1:10){
    if (CorrelationScale){
      segments(0, i/10, 1, i/10, col= 'gray')
      segments(i/10, 0, i/10, 1, col= 'gray')
      text((i/10), 0, labels=(i/10), cex=0.7, pos=1, col="red")
      text(0, (i/10), labels=(i/10), cex=0.7, pos=2, col="red")
      }
    
    if (ContributionScale){
      segments(0, sqrt(i/10), 1, sqrt(i/10), col= 'gray')
      segments(sqrt(i/10), 0, sqrt(i/10), 1, col= 'gray')
      text(sqrt(i/10), 1, labels=(i/10), cex=0.7, pos=3, col="gray")
      text(1, sqrt(i/10), labels=(i/10), cex=0.7, pos=4, col="gray")
      QuarterCircle(sqrt(i/10), lty=1, color="gray", lwd=1.5)
      text(sqrt(i/10)/sqrt(2), sqrt(i/10)/sqrt(2), srt=45, labels=(i/10), cex=0.7, pos=2, col="gray")
    }
  }
  
  
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