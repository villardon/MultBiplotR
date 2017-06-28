textsmart <- function(A, Labels, CexPoints=1, ColorPoints="black", ...) {
  n = dim(A)[1]
  
  if (length(CexPoints == 1)) 
    CexPoints = rep(CexPoints, n)
  
  if (length(ColorPoints == 1)) 
    ColorPoints = rep(ColorPoints, n)
  for (i in 1:n) {
    if (A[i, 1] > 0) 
      markerpos = 4
    else markerpos = 2
    text(A[i, 1], A[i, 2], rownames(A)[i], cex=CexPoints, col = ColorPoints[i], pos = markerpos, offset = 0.2, ...)
  }
}

