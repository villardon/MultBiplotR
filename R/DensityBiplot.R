DensityBiplot <- function(X, y = NULL, grouplabels = NULL, ncontours = 6, groupcolors = NULL, ncolors=20, ColorType=4) {
  if (is.null(y))
    y = matrix(1, dim(X)[1], 1)
  
  if (is.factor(y)) y=as.integer(y)
  
  nlevels = max(y)
  
  if (is.null(grouplabels)) grouplabels=1:nlevels
  
  if (is.null(groupcolors)) groupcolors = 1:nlevels 
  
  xliml=min(X[,1])
  xlimu=max(X[,1])
  yliml=min(X[,2])
  ylimu=max(X[,2])
  
  switch(ColorType, "1" = {colores = rainbow(ncolors)},
         "2" = {colores = heat.colors(ncolors)},
         "3" = {colores = terrain.colors(ncolors)},
         "4" = {colores = topo.colors(ncolors)},
         "5" = {colores = cm.colors(ncolors)})
  #print(colores)
  
  for (i in 1:nlevels) {
    x = X[which(y == i), ]
    if(!is.null(nrow(x))){
      f1 <- kde2d(x[, 1], x[, 2], n = 400, lims = c(xliml, xlimu,yliml, ylimu))
      contours = round(seq(0, max(f1$z), length.out = ncontours), digits = 2)
      contour(f1, levels = contours, col = groupcolors[i], add = TRUE, lwd=2)
    }
}
}
