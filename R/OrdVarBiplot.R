OrdVarBiplot <- function(bi1, bi2,threshold, xmin = -3, xmax = 3, ymin = -3, ymax = 3, label = "Point", mode = "a", CexMarks=0.7,  CexPoint = 0.8, PchPoint = 1, Color = "green", 
                           tl = 0.03, textpos=1, ...){

  b1 = bi1/(bi1^2 + bi2^2)
  b2 = bi2/(bi1^2 + bi2^2)
  b = b2/b1

  ang = atan(bi2/bi1) * 180/pi
  
  x1 = xmin
  y1 = b * xmin
  if ((y1 > ymin - 0.001) & (y1 < ymax + 0.001)) 
    if ((x1 * b1 + y1 * b2) < 0) 
      ini = c(x1, y1)
  else {
    final = c(x1, y1)
  }
  x1 = xmax
  y1 = b * xmax
  if ((y1 > ymin - 0.001) & (y1 < ymax + 0.001)) 
    if ((x1 * b1 + y1 * b2) < 0) 
      ini = c(x1, y1)
  else {
    final = c(x1, y1)
  }
  x1 = ymin/b
  y1 = ymin
  if ((x1 > xmin) & (x1 < xmax)) 
    if ((x1 * b1 + y1 * b2) < 0) 
      ini = c(x1, y1)
  else {
    final = c(x1, y1)
  }
  x1 = ymax/b
  y1 = ymax
  if ((x1 > xmin) & (x1 < xmax)) 
    if ((x1 * b1 + y1 * b2) < 0) 
      ini = c(x1, y1)
  else {
    final = c(x1, y1)
  }
  if (mode == "p") {
    points(b1, b2, pch = 16, col = Color, cex = CexPoint, ...)
    c1 = b1
    c2 = b2
  }
  if (mode == "a") {
    arrows(0, 0, b1, b2, length = 0.1, angle = 20, col = Color, ...)
    c1 = b1
    c2 = b2
  }
  if (mode == "ah") {
    arrows(0, 0, b1, b2, length = 0.1, angle = 20, col = Color, ...)
    c1 = final[1]
    c2 = final[2]
  }
  if (mode == "b") {
    lines(rbind(ini, final), col = Color, ...)
    c1 = final[1]
    c2 = final[2]
  }
  if (mode == "h") {
    lines(rbind(c(0, 0), final), col = Color, ...)
    c1 = final[1]
    c2 = final[2]
  }
  
  
  if (mode == "s") {
    lines(rbind(ini, final), col = Color, lwd = 1, lty = 1, ...)
    OrCoor=OrdVarCoordinates(tr=threshold, c(bi1,bi2), plotresponse=F)
    c1 = final[1]
    c2 = final[2]
    ang = atan(bi2/bi1) * 180/pi
    M=OrCoor$points
    k=dim(M)[1]
    deltax <- tl * sin(ang * pi/180)
    deltay <- tl * cos(ang * pi/180)
    Mn <- cbind(M[, 1] + deltax, M[, 2] - deltay)
    for (i in 1:k) {
      if (InBox(M[i, 1], M[i, 2], xmin, xmax, ymin, ymax)) {
        lines(rbind(M[i, 1:2], Mn[i, 1:2]), col = Color, lwd = 1, ...)
        text(Mn[i, 1], Mn[i, 2], OrCoor$labels[i], pos = 1, offset = 0.2, cex = CexMarks, srt = ang, col = Color, ...)
      }
    }
    
  }

  
  if ((b1 < 0) & (textpos==1)) textpos=3
  if ((b1 < 0) & (textpos==2)) textpos=4  
  if ((b1 < 0) & (textpos==3)) textpos=1
  if ((b1 < 0) & (textpos==4)) textpos=2

  text(c1, c2, label, cex = CexPoint, pos = textpos , offset = 0.3, srt = ang, col = Color, ...)

}