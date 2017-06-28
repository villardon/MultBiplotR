OrdVarCoordinates <- function(tr, b = c(1, 1), inf = -12, sup = 12, step = 0.01, plotresponse = FALSE, label="Item", labx= "z", laby="Probability", catnames=NULL, Legend=TRUE, LegendPos=1) {
  ncat = length(tr) + 1
  if (is.null(catnames)) catnames=as.character(1:ncat)
  b = matrix(b, 1, 2)
  tr=matrix(tr, 1, ncat-1)
  z = seq(from = inf, to = sup, by = step)
  npoints = length(z)
  eta = matrix(1, npoints, 1) %*% tr - z %*% matrix(1, 1, ncat - 1)

  pa = cbind(exp(eta)/(1 + exp(eta)), matrix(1, npoints, 1))
  p = matrix(0, npoints, ncat)
  p[, 1] = pa[, 1]
  p[, 2:ncat] = pa[, 2:ncat] - pa[, 1:ncat - 1]
  
  maxima = apply(p, 1, which.max)

  intersec = which(maxima[2:npoints] - maxima[1:npoints - 1] > 0)
  res = list()
  seno=b[1]/sum(b^2)
  coseno=b[2]/sum(b^2)
  res$z = matrix((z[intersec] + z[intersec + 1])/2, length(intersec), 1)
  x=res$z*seno
  y=res$z*coseno
  res$points = cbind(x,y)
  if (b[1]>0)
    res$labels = paste(maxima[intersec], "-", maxima[intersec + 1], sep = "")
  else
    res$labels = paste(maxima[intersec+1], "-", maxima[intersec], sep = "")
  if (length(intersec) == (ncat - 1)) {
    res$hidden = FALSE
  } else {
    res$hidden = TRUE
    res$cathidden = NULL
    for (i in 1:ncat) if (length(which(maxima == i)) == 0) 
      res$cathidden = c(res$cathidden, i)
  }
  class(res) = "OrdVarCoord"
  
  if (plotresponse == TRUE) {
    plot(z, p[, 1], type = "l", col=1, ylim=c(0, 1) , main=label, xlab=labx, ylab=laby)
    for (i in 2:(ncat)) points(z, p[, i], type = "l", col=i)
    
    if (Legend){
    LegendPositions= c("bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right", "center")
    if (is.numeric(LegendPos)) LegendPos=LegendPositions[LegendPos]
      legend(x=LegendPos, catnames , col = 1:ncat, text.col = "black", pch = 16)}
    
  }

  return(res)
}