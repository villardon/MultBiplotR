# Autor: Jose Luis Vicente Villardon
# Dpto. de Estadistica
# Universidad de Salamanca
# Revisado: Noviembre/2013

GowerSimilarities <- function(x, y=NULL, transformation="sqrt(1-S)") {
  
  if (is.data.frame(x)) {
    n = dim(x)[1]
    p = dim(x)[2]
    
    if (is.null(y)) y=x
    
    if (is.data.frame(y)) {
      n1 = dim(y)[1]
      p1 = dim(y)[2]
      
      if (!(p==p1)) stop("Number of columns of the two matrices are not the same")
      
    } else stop("Suplementary data is not organized as a data frame")
  } else stop("Main data is not organized as a data frame")
  
  sim=matrix(0,n,n1)
  
  for (i in 1:n) for (j in 1:n1) {
    similarity=matrix(0,p,1)
    peso=matrix(1,p,1)
    for (k in 1:p) {
      if (is.factor(x[[k]])){
        if (x[i, k] == y[j, k]) similarity[k] = 1
        if ((length(levels(x[[k]]))==2) & (as.integer(x[i, k])==1) & (as.integer(y[j, k])==1)) peso[k]=0}
      
      if (is.numeric(x[[k]])) similarity[k]= 1- (abs(x[i, k] - y[j, k]) / (max(c(x[, k],y[,k]))-min(c(x[, k],y[,k]))))
      
      if (is.integer(x[[k]])){
        if (x[i, k] == y[j, k]) similarity[k] = 1
        if ((x[i, k]==0) & (y[j, k]==0)) peso[k]=0}
    }
    
    sim[i, j] = sum(similarity * peso)/sum(peso)
  }
  
  switch(transformation, `Identity` = {
    dis=sim
  }, `Identity` = {
    dis=sim
  }, `1-S` = {
    dis=1-sim
  }, `sqrt(1-S)` = {
    dis = sqrt(1 - sim)
  }, `-log(s)` = {
    dis=-1*log(sim)
  }, `1/S-1` = {
    dis=1/sim -1 
  }, `sqrt(2(1-S))` = {
    dis== sqrt(2*(1 - sim))
  }, `1-(S+1)/2` = {
    dis=1-(sim+1)/2
  }, `1-abs(S)` = {
    dis=1-abs(sim)
  }, `1/(S+1)` = {
    dis=1/(sim)+1
  })
  
  return(dis)
  
}



  
  