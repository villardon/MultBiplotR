PlotBinaryVar <- function(b0, bi1, bi2, xmin = -3, xmax = 3, ymin = -3, ymax = 3, label = "Point", mode = mode, CexPoint = 0.8, PchPoint = 1, Color = "blue",
                          tl = 0.01, ts = "Probability", Position="Angle",  ...){
  b1 = bi1/(bi1^2 + bi2^2)
  b2 = bi2/(bi1^2 + bi2^2)
  b = b2/b1
  x1 = xmin
  y1 = b * xmin
  
  if ((y1 > ymin - 0.001) & (y1 < ymax + 0.001)) 
    if ((x1 * b1 + y1 * b2) < 0) 
      ini = c(x1, y1)
  else {
    final = c(x1, y1)
    markerpos = 2
    angle = 0
  }
  x1 = xmax
  y1 = b * xmax
  if ((y1 > ymin - 0.001) & (y1 < ymax + 0.001)) 
    if ((x1 * b1 + y1 * b2) < 0) 
      ini = c(x1, y1)
  else {
    final = c(x1, y1)
    markerpos = 4
    angle = 0
  }
  x1 = ymin/b
  y1 = ymin
  if ((x1 > xmin) & (x1 < xmax)) 
    if ((x1 * b1 + y1 * b2) < 0) 
      ini = c(x1, y1)
  else {
    final = c(x1, y1)
    markerpos = 4
    angle = 270
  }
  x1 = ymax/b
  y1 = ymax
  if ((x1 > xmin) & (x1 < xmax)) 
    if ((x1 * b1 + y1 * b2) < 0) 
      ini = c(x1, y1)
  else {
    final = c(x1, y1)
    markerpos = 4
    angle = 90
  }
  
  c1 = final[1]
  c2 = final[2]
  
  c=c(0,0)
  c[1] = bi1/(bi1^2+bi2^2)
  c[2] = bi2/(bi1^2+bi2^2)
  b=matrix(0,11,2)
  b[1,] = (logit(0.1) - b0) * c
  b[2,] = (logit(0.2) - b0) * c
  b[3,] = (logit(0.3) - b0) * c
  b[4,] = (logit(0.4) - b0) * c
  b[5,] = (logit(0.5) - b0) * c
  b[6,] = (logit(0.6) - b0) * c
  b[7,] = (logit(0.7) - b0) * c
  b[8,] = (logit(0.8) - b0) * c
  b[9,] = (logit(0.9) - b0) * c
  b[10,] = (logit(0.25) - b0) * c
  b[11,] = (logit(0.75) - b0) * c
  
  if (mode == "p") {
    points(b1, b2, pch = 16, col = Color, cex = CexPoint, ...)
    c1 = b1
    c2 = b2
    if (b1 < 0) 
      markerpos = 2
    else markerpos = 4
    angle = 0
  }
  
  if (mode == "a") {
    points(b[5,1], b[5,2], col=Color, pch=16)
    arrows(b[5, 1], b[5, 2], b[11, 1], b[11, 2], length = 0.1, angle = 20, col = Color)
    c1 = b[11, 1]
    c2 = b[11, 2]
    if (c1 < 0) 
      markerpos = 2
    else markerpos = 4
    angle = 0
  }
  
  if (mode == "ah") {
    arrows(b[5, 1], b[5, 2], b[11, 1], b[11, 2], length = 0.1, angle = 20, col = Color)
    c1 = final[1]
    c2 = final[2]
  }
  
  if (mode == "b") {
    lines(rbind(ini, final), col = Color, ...)
    if (InBox(b[5,1],  b[5,2], xmin, xmax, ymin, ymax))
      points(b[5,1], b[5,2], col=Color, pch=16)
    c1 = final[1]
    c2 = final[2]
  }
  
  if (mode == "s") {
    
    lines(rbind(ini, final), col = Color, lwd = 1, lty = 1, ...)
    
    if (InBox(b[5,1],  b[5,2], xmin, xmax, ymin, ymax))
      points(b[5,1], b[5,2], col=Color, pch=16)
    c1 = final[1]
    c2 = final[2]
    ang = atan(bi2/bi1) * 180/pi
    k = 9
    ticks=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    ticklabels=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
    M = b[1:9,]
    deltax <- tl * sin(ang * pi/180)
    deltay <- tl * cos(ang * pi/180)
    Mn <- cbind(M[, 1] + deltax, M[, 2] - deltay)
    for (i in 1:k) {
      if (InBox(M[i, 1], M[i, 2], xmin, xmax, ymin, ymax)) {
        lines(rbind(M[i, 1:2], Mn[i, 1:2]), col = Color, lwd = 1, ...)
        text(Mn[i, 1], Mn[i, 2], ticklabels[i], pos = 1, offset = 0.2, cex = 0.5, srt = ang, col = Color, ...)
      }
    }
    
  }
  
  if (Position == 'Angle'){
    angle=atan(bi2/bi1) * 180/pi
    if (bi1 < 0) 
      markerpos = 2
    else markerpos = 4
  }
  
  
  text(c1, c2, label, cex = CexPoint, pos = markerpos, offset = 0.2, srt = angle, col = Color, ...)
  
}