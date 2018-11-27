HistogramPanel <- function(X, groups=NULL, nrows=NULL, separated=FALSE, notch=FALSE){

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

  if (is.null(nrows))
    nrows=round(sqrt(p))

  ncols=ceiling(p/nrows)

  if (separated==FALSE)
    op=par(mfrow=c(nrows, ncols))

  for (j in 1:p){
    if (separated==TRUE) dev.new()
    boxplot(X[,j]~groups,notch=notch, main=varnames[j])}

  if (separated==FALSE)
    par(op)

}
