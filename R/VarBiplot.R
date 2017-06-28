VarBiplot <- function(bi1, bi2, b0=0, xmin = -3, xmax = 3, ymin = -3, ymax = 3, label = "Point", mode = "a", 
                      CexPoint = 0.8, PchPoint = 1, Color = "blue", 
                      ticks = c(-3, -2.5, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 2.5, 3), ticklabels = round(ticks, digits = 2),
                      tl = 0.03, ts = "Complete", Position="Angle", AddArrow=FALSE, ...) {
  # mode a= arrows b = both extremes h= half (just the positive part)
  #      s= scale 
  # Plots a var on a biplot

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
  if (mode == "p") {
    points(bi1, bi2, pch = 16, col = Color, cex = CexPoint, ...)
    c1 = bi1
    c2 = bi2
    if (bi1 < 0) 
      markerpos = 2
    else markerpos = 4
    angle = 0
  }
  if (mode == "a") {
    arrows(0, 0, bi1, bi2, length = 0.1, angle = 20, col = Color, ...)
    c1 = bi1
    c2 = bi2
    if (bi1 < 0) 
      markerpos = 2
    else markerpos = 4
    angle = 0
  }
  if (mode == "ah") {
    arrows(0, 0, bi1, bi2, length = 0.1, angle = 20, col = Color, ...)
    lines(rbind(c(0, 0), final), col = Color, lty=3, ...)
    c1 = final[1]
    c2 = final[2]
  }
  if (mode == "b") {
    lines(rbind(ini, final), col = Color, ...)
    arrows(0, 0, bi1, bi2, length = 0.1, angle = 20, col = Color, lwd=2 ,  ...)
    c1 = final[1]
    c2 = final[2]
  }
  if (mode == "h") {
    lines(rbind(c(0, 0), final), col = Color, ...)
    c1 = final[1]
    c2 = final[2]
  }
  if (mode == "s") {
    if (ts == "BoxPlot") {
      lty = 3
    } else {
      lty = 1
    }
    lines(rbind(ini, final), col = Color, lwd = 1, lty = lty, ...)
    c1 = final[1]
    c2 = final[2]
    ang = atan(bi2/bi1) * 180/pi
    k = length(ticks)
    M = cbind((ticks - b0) * b1, (ticks - b0) * b2)
    deltax <- tl * sin(ang * pi/180)
    deltay <- tl * cos(ang * pi/180)
    Mn <- cbind(M[, 1] + deltax, M[, 2] - deltay)
    for (i in 1:k) {
      if (InBox(M[i, 1], M[i, 2], xmin, xmax, ymin, ymax)) {
        lines(rbind(M[i, 1:2], Mn[i, 1:2]), col = Color, lwd = 1, ...)
        text(Mn[i, 1], Mn[i, 2], ticklabels[i], pos = 1, offset = 0.2, cex = 0.5, srt = ang, col = Color, ...)
      }
    }
    if (ts == "BoxPlot") {
      points(M[3, 1], M[3, 2], pch = 16, col = Color, ...)
      lines(rbind(M[2, ], M[4, ]), col = Color, lwd = 2, lty = 1, ...)
    }
  }
  
  if (AddArrow) arrows(0, 0, bi1, bi2, length = 0.1, angle = 20, col = Color, lwd = 2, ...)
  
  if (Position == 'Angle'){
    angle=atan(bi2/bi1) * 180/pi
    if (bi1 < 0) 
      markerpos = 2
    else markerpos = 4
   }

  text(c1, c2, label, cex = CexPoint, pos = markerpos, offset = 0.2, srt = angle, col = Color, ...)
}

