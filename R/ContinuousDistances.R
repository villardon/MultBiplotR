# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013
# Gower, J. (2006) Similarity, dissimilarity and distance - Measures of. Encyclopedia of Statistical Sciences. Vol12. 2nd Edition. Wiley.
ContinuousDistances <- function(x, y=NULL,  coef = "Pythagorean", r = 1) {
	distances = c("Pythagorean", "Taxonomic", "City", "Minkowski", "Divergence", "dif_sum", "Camberra", "Bray_Curtis", "Soergel", "Ware_Hedges", "Gower")
	if (is.numeric(coef)) coef = distances[coef]
  if (is.null(y)) y=x
	n = nrow(x)
	p = ncol(x)
  s = nrow(y)
  q=ncol(y)
  NamesX=rownames(x)
  NamesY=rownames(y)
  
  print(NamesX)
  
  if (coef=="Gower") rank=apply(rbind(y,x),2,max)-apply(rbind(y,x),2,min)
  if (!(p==q)) stop("The matrices should have the same number of columns")

dis=matrix(0,s,n)
	for (i in 1:s) for (j in 1:n) {
	  switch(coef, Pythagorean = {
	    dis[i, j] = sqrt(sum((y[i, ] - x[j, ])^2))
	  },Taxonomic = {
	    dis[i,j]=sqrt(sum(((y[i,]-x[j,])^2)/r^2))
	  },City = {
	    dis[i,j]=sum(abs(y[i,]-x[j,]))
	  },Minkowski = {
	    dis[i,j]=(sum(abs(y[i,]-x[j,])^r))^(1/r)
	  },Divergence = {
	    dis[i,j]=sqrt(sum((y[i,]-x[j,])^2/(y[i,]+x[j,])^2))
	  },dif_sum = {
	    dis[i,j]=sum(abs(y[i,]-x[j,])/abs(y[i,]+x[j,]))
	  },Camberra = {
	    dis[i,j]=sum(abs(y[i,]-x[j,])/(abs(y[i,])+abs(x[j,])))
	  },Bray_Curtis = {
	    dis[i,j]=sum(abs(y[i,]-x[j,]))/sum(y[i,]+x[j,])
	  },Soergel = {
	    dis[i,j]=sum(abs(y[i,]-x[j,]))/sum(apply(rbind(y[i,],x[j,]),2,max))
	  },Ware_Hedges = {
	    dis[i,j]=sum(1-apply(rbind(y[i,],x[j,]),2,min)/apply(rbind(y[i,],x[j,]),2,max))
	  },Gower = {
	    dis[i,j]=sum(abs(y[i,]-x[j,])/rank)
	  })
	}

  rownames(dis)=NamesY
  colnames(dis)=NamesX
	class(dis) = "proximities"
	return(dis)
}
