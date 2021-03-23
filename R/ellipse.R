ellipse <- function(center=c(0,0), radius=c(1,1), angle= NULL, orient= diag(2), add=TRUE, npoints=100, ...){
  # angle in radians
  if (!is.null(angle)) orient=matrix(c(cos(angle), sin(angle), -1* sin(angle), cos(angle)), 2, 2)
  ang=seq(from=0, to=2*pi, by =(2*pi/npoints))
  z=matrix(0,(npoints+1), 2)
  for (j in 1:(npoints+1)){
    z[j,1]=cos(ang[j])
    z[j,2]=sin(ang[j])}
  
  z=z %*% diag(radius) %*% orient
  z[,1]=z[,1]+center[1]
  z[,2]=z[,2]+center[2]
  
  if (add) 
    points(z[,1], z[,2], type="l", ...)
  else
    plot(z[,1], z[,2], type="l", ...)
}