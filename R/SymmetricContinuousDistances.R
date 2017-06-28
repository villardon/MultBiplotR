SymmetricContinuousDistances <- function(x, coef = "Pythagorean", r = 1) {
  distances = c("Pythagorean", "Taxonomic", "City", "Minkowski", "Divergence", "dif_sum", "Camberra", "Bray_Curtis", "Soergel", "Ware_Hedges")
  if (is.numeric(coef)) coef = distances[coef]
  n = nrow(x)
  p = ncol(x)
  NamesX=rownames(x)
  dis=matrix(0,n,n)
  for (i in 1:n) for (j in i:n) {
    switch(coef, Pythagorean = {
      dis[i, j] = sqrt(sum((x[i, ] - x[j, ])^2))
    },Taxonomic = {
      dis[i,j]=sqrt(sum(((x[i,]-x[j,])^2)/r^2))
    },City = {
      dis[i,j]=sum(abs(x[i,]-x[j,]))
    },Minkowski = {
      dis[i,j]=(sum(abs(x[i,]-x[j,])^r))^(1/r)
    },Divergence = {
      dis[i,j]=sqrt(sum((x[i,]-x[j,])^2/(x[i,]+x[j,])^2))
    },dif_sum = {
      dis[i,j]=sum(abs(x[i,]-x[j,])/abs(x[i,]+x[j,]))
    },Camberra = {
      dis[i,j]=sum(abs(x[i,]-x[j,])/(abs(x[i,])+abs(x[j,])))
    },Bray_Curtis = {
      dis[i,j]=sum(abs(x[i,]-x[j,]))/sum(x[i,]+x[j,])
    },Soergel = {
      dis[i,j]=sum(abs(x[i,]-x[j,]))/sum(apply(rbind(x[i,],x[j,]),2,max))
    },Ware_Hedges = {
      dis[i,j]=sum(1-apply(rbind(x[i,],x[j,]),2,min)/apply(rbind(x[i,],x[j,]),2,max))
    })
    dis[j,i] = dis[i,j]
  }
  rownames(dis)=NamesX
  colnames(dis)=NamesX
  class(dis) = "distances"
  return(dis)
}