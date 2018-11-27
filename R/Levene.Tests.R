Levene.Tests <- function(X, groups=NULL){
  n=dim(X)[1]
  p=dim(X)[2]

  if (is.null(groups)) {
    groups=as.factor(rep(1,n))
    levels(groups)="Complete Sample"}

  if (!is.factor(groups)) stop("The variable defining the groups must be a factor")

  g=length(levels(groups))
  Levels=levels(groups)
  varnames=colnames(X)

  res=apply(X,2,leveneTest, group=groups)
  Levene=matrix(0, p, 4)
  for (i in 1:p){
    Levene[i,1]=res[[i]]$`F value`[1]
    Levene[i,2]=res[[i]]$Df[1]
    Levene[i,3]=res[[i]]$Df[2]
    Levene[i,4]=res[[i]]$`Pr(>F)`[1]
  }
  rownames(Levene)=varnames
  colnames(Levene)=c("F-value", "d.f.1", "d.f.2", "p-value")
  return(Levene)
}
