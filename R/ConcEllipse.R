ConcEllipse <- function(data,confidence=1,npoints=100){
  center=apply(data,2,mean)
  n=dim(data)[1]
  datacen=(diag(n)-matrix(1, n, n)/n) %*% data
  DE=svd(datacen)
  K=DE$u
  d=matrix(0,n,1)
  for (i in 1:n)
    d[i]=sqrt(sum(K[i,]^2))
  
  radius = quantile(d,confidence)
  retained=which(d<radius)
  ang=seq(from=0, to=2*pi, by =(2*pi/npoints))
  z=matrix(0,(npoints+1), 2)
  for (j in 1:(npoints+1)){
    z[j,1]=cos(ang[j])
    z[j,2]=sin(ang[j])}
  z=radius*z   
  z=z %*% diag(DE$d) %*% t(DE$v)
  z[,1]=z[,1]+center[1]
  z[,2]=z[,2]+center[2]
  res=list(data=data, confidence=confidence, ellipse=z, center=center, retained=retained, distances=d, radius=radius)
  class(res)="ellipse"
  return(res)
}