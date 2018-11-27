NormalityTests <- function(X, groups=NULL, plot=FALSE, SortByGroups=FALSE){

  k=0
  n=dim(X)[1]
  p=dim(X)[2]

  if (is.null(groups)) {
    groups=as.factor(rep(1,n))
    levels(groups)="Complete Sample"}

  if (!is.factor(groups)) stop("The variable defining the groups must be a factor")

  g=length(levels(groups))
  Levels=levels(groups)
  varnames=colnames(X)

  Descriptives=list()

  if (SortByGroups==TRUE){
    if (plot) op=par(mfrow=c(g, p))
    for (i in 1:g){
      XX = as.matrix(X[which(groups == Levels[i]), ])
      swt=apply(XX, 2, shapiro.test)
      SW=matrix(0,p,2)
      for (j in 1:p){
        SW[j,]=c(swt[[j]]$statistic, swt[[j]]$p.value)
        if (plot) qqnorm(XX[,j] , main=paste(Levels[i],varnames[j], sep="-"))
        }
      rownames(SW)=varnames
      colnames(SW) = c("Statistic", "p-value")
      Descriptives[[i]]=SW
    }
    if (plot) par(op)

  names(Descriptives)= paste(Levels, "-", "Shapiro-Wilk normality tests")
  }

  if (SortByGroups==FALSE){
    if (plot) op=par(mfrow=c(p, g))
    for (i in 1:p){
      XX = as.numeric(X[,i])
      swt=by(XX,groups,shapiro.test)
      SW=matrix(0,g,2)
      for (j in 1:g){
       SW[j,]=c(swt[[j]]$statistic, swt[[j]]$p.value)
       if (plot) qqnorm(XX[which(groups == Levels[j])] , main=paste(varnames[i], Levels[j], sep="-"))}
      rownames(SW)=Levels
      colnames(SW) = c("Statistic", "p-value")
      Descriptives[[i]]=SW

    }
    if (plot) {
      title("Normal Q-Q Plot", outer=TRUE)
      par(op)}
    names(Descriptives)= paste(varnames, "-", "Shapiro-Wilk normality tests")
  }
  return(Descriptives)

}
