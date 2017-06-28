SelectVarsAnova <- function(X, group, Cut="Percentile", alpha=0.01 ){
  Cuts=c("Percentile", "Bonferroni")
  X=as.matrix(X)
  g = length(levels(group))
  n = dim(X)[1]
  p = dim(X)[2]
  
  media=apply(X,2,mean)
  RowNames=rownames(X)
  ColNames=colnames(X)
  
  X2=X- matrix(1,n,1) %*% matrix(media,1,p)
  
  G = Factor2Binary(group)
  ng=diag(t(G) %*% G)
  XB= diag(1/ng) %*% t(G) %*% X2 
  B = sqrt(diag(ng)) %*% XB
  
  TSS=apply(X2^2,2,sum)
  BSS=apply(B^2,2,sum)
  WSS=TSS-BSS
  Fexp=(BSS/(g-1))/(WSS/(n-g))
  pval=1- pf(Fexp, g-1, n-g)
  
  
  switch(Cut, Percentile = {
    alpha=quantile(pval, alpha)
  },Bonferroni = {
    alpha=alpha/p
  })
  
  selecc=which(pval<alpha)
  X2=data.frame(group,X[,selecc])
  op <- par(mfrow=c(2,1))
  hist(pval)
  boxplot(pval)
  par(op)
 return(X2)
}
