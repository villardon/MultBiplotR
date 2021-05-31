HistogramPanel <- function(X, nrows=NULL, separated=FALSE, ...){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  k=0
  n=dim(X)[1]
  p=dim(X)[2]
  varnames=colnames(X)

  if (is.null(nrows))
    nrows=round(sqrt(p))

  ncols=ceiling(p/nrows)

  if (separated==FALSE)
    op=par(mfrow=c(nrows, ncols))

  for (j in 1:p){
    if (separated==TRUE) dev.new()
    hist(X[,j], main=varnames[j],...)}

  if (separated==FALSE)
    par(op)
}
