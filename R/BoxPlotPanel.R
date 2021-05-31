BoxPlotPanel <- function(X, groups=NULL, nrows=NULL, panel=TRUE, notch=FALSE, GroupsTogether=TRUE, ...){
  separated=!panel
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
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

  if (GroupsTogether){

    if (is.null(nrows))
      nrows=round(sqrt(p))
print(nrows)
    ncols=ceiling(p/nrows)

    if (separated==FALSE)
      op=par(mfrow=c(nrows, ncols))

    for (j in 1:p){
      if (separated==TRUE) dev.new()
      boxplot(X[,j]~groups,notch=notch, main=varnames[j], ...)}

    if (separated==FALSE)
      par(op)
  }
  else{
    if (is.null(nrows))
      nrows=round(sqrt(g))

    ncols=ceiling(g/nrows)

    if (separated==FALSE)
      op=par(mfrow=c(nrows, ncols))

    for (i in 1:g){
      XX = as.matrix(X[which(groups == Levels[i]), ])
      x=numeric()
      grupvar=numeric()
      for (j in 1:p){
        x=c(x, as.numeric(XX[,j]))
        grupvar=c(grupvar, rep(j, length(XX[,j])))
      }
      grupvar=as.factor(grupvar)
      levels(grupvar)=varnames
      if (separated==TRUE) dev.new()
      boxplot(x~grupvar,notch=notch, main=Levels[i], ...)}
  }

}
