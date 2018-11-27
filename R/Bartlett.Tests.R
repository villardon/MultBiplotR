Bartlett.Tests <- function(X, groups=NULL){
  n=dim(X)[1]
  p=dim(X)[2]

  if (is.null(groups)) {
    groups=as.factor(rep(1,n))
    levels(groups)="Complete Sample"}

  if (!is.factor(groups)) stop("The variable defining the groups must be a factor")

  g=length(levels(groups))
  Levels=levels(groups)
  varnames=colnames(X)

  res=apply(X,2,bartlett.test, g=groups)
  Bartlett=matrix(0, p, 3)
 for (i in 1:p){
   Bartlett[i,1]=res[[i]]$statistic
   Bartlett[i,2]=res[[i]]$parameter
   Bartlett[i,3]=res[[i]]$p.value
 }
  rownames(Bartlett)=varnames
  colnames(Bartlett)=c("Chi-squared", "d.f.", "p-value")
  return(Bartlett)
}

