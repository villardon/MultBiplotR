plot.MGC <- function(x, vars=NULL, groups=x$Classification, CexPoints=0.2, Confidence=0.95, ...){
  n=dim(x$Data)[1]
  nvars=dim(x$Data)[2]
  if (is.null(vars)) vars=1:nvars
  if (length(vars)==1){
    d=matrix(0,200,x$NG)
    dt=matrix(0,200,1)
    minimo=min(x$Data[,vars])
    maximo=max(x$Data[,vars])
    x1=matrix(seq(minimo,maximo, length.out =200), 200,1)
    for (i in 1:x$NG){
      d[,i] <-multnormdens(x1, mean=x$Centers[i,vars], sigma=x$Covariances[[i]][vars,vars])
      dt=dt+d[,i]*x$GroupProbabilities[i]
    }
    plot(x1,dt, type="l", col=1, ...)
    par(usr=c(minimo, maximo, 0, max(d)*1.5))
    for (i in 1:x$NG){
      points(x1,d[,i], type="l", col=i+1, ...)
    }
  }
  
  if (length(vars)==2){ 
    dat=x$Data[,vars]
    gg=matrix((2:(x$NG+1)), x$NG,1)
    col1=t(col2rgb(gg))/255
    COL=x$P %*% col1
    colores=rgb(COL)
    colores[which(x$Outliers==1)]="#000000"
    
    plot(dat,  col= colores, cex = CexPoints+0.2, pch=16)
    for (i in 1:x$NG)
      PlotGaussianEllipse(mean=x$Centers[i,vars], sigma=x$Covariances[[i]][vars,vars], n=n, col=(i+1), ... )
  }
  
  if (length(vars)>2){ 
    dat=as.data.frame(x$Data)
    dat=dat[,vars]
    gg=matrix((2:(x$NG+1)), x$NG,1)
    col1=t(col2rgb(gg))/255
    COL=x$P %*% col1
    colores=rgb(COL)
    colores[which(x$Outliers==1)]="#000000"
    #pairs(dat, diag.panel = panel.NPdens, col= colores, main="MGC Plot", cex = CexPoints, ...)
    #splom(dat | (x$Classification))
  }
}


panel.NPdens <- function(x, groups=NULL, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  if (is.null(groups)){
  d <- density(x)
  par(usr = c(usr[1:2], 0, max(d$y*1.4)) )
  points(d, type="l", ...)}
  
}

PlotGaussianEllipse <- function(mean=c(0,0), sigma=diag(2) , n= 1000, npoints=100, confidence=0.95, add=TRUE, distF=TRUE, ...){
  DE=svd(sigma)
  K=DE$u
  radius = sqrt((2*(n-1)/(n-2))*qf(confidence,2,(n-2)))
  ang=seq(from=0, to=2*pi, by =(2*pi/npoints))
  z=matrix(0,(npoints+1), 2)
  for (j in 1:(npoints+1)){
    z[j,1]=cos(ang[j])
    z[j,2]=sin(ang[j])}
  z=radius*z   
  z=z %*% sqrt(diag(DE$d)) %*% t(DE$v)
  z[,1]=z[,1]+mean[1]
  z[,2]=z[,2]+mean[2]
  if (add) points(z, type="l", ...)
  else  plot(z, type="l", ...)
}